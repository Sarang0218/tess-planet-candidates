"""Archival-image check for a background star behind a high-proper-motion host.  Downloads Digitized Sky Survey and 2MASS cutouts from
NASA SkyView, finds the star at its old position (which dates the image), then measures the flux at the position the star occupies during
the TESS observations.  A background eclipsing binary bright enough to mimic the transit would show up there.
Usage: deep2_archive.py TIC K"""
import sys, os, io, json, urllib.request, urllib.parse, numpy as np, pandas as pd, warnings
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')
R = '/home/sarang/tess_search'; TIC, K = int(sys.argv[1]), int(sys.argv[2]); D = f'{R}/deep2/archive/{TIC}'; os.makedirs(D, exist_ok=True)
v = json.load(open(f'{R}/deep2/vet/{TIC}_{K}.json')); depth = v['trap_depth_ppm'] * 1e-6
g = pd.read_csv(f'{R}/deep2/pix/{TIC}/tic_120arcsec.csv').drop_duplicates('ID'); a = g[g.ID == TIC].iloc[0]; pmE, pmN = a.pmRA / 1000, a.pmDEC / 1000; pm = np.hypot(pmE, pmN)
src = json.load(open(f'{R}/deep2/lc/{TIC}.meta.json'))['src']; secs = sorted(int(s) for s in src); ep_tess = 2018.55 + (np.mean(secs) - 1) * 27.4 / 365.25        # mean date of the TESS data
dmax = -2.5 * np.log10(depth)                              # a fully eclipsing background star must be within this many magnitudes of the target
print(f'TIC {TIC}: G = {a.GAIAmag:.2f}, proper motion {pm*1000:.0f} mas/yr; mean TESS epoch {ep_tess:.1f}; a background star must be brighter than G ~ {a.GAIAmag + dmax:.1f} (a 100 % eclipse) to mimic a {depth*1e6:.0f} ppm dip')
def off(yr): return pmE * (yr - 2000.0), pmN * (yr - 2000.0)   # arcsec E, N from the epoch-2000 catalogue position
res = []; panels = []
for sv in ('DSS1 Red', 'DSS1 Blue', 'DSS2 Red', 'DSS2 Blue', 'DSS2 IR', '2MASS-J'):
    fn = f'{D}/{sv.replace(" ", "_")}.fits'
    if not os.path.exists(fn):
        url = 'https://skyview.gsfc.nasa.gov/current/cgi/runquery.pl?' + urllib.parse.urlencode(dict(Position=f'{a.ra},{a.dec}', Survey=sv, Pixels=180, Size=0.05, Return='FITS', Scaling='Linear'))
        try: open(fn, 'wb').write(urllib.request.urlopen(url, timeout=120).read())
        except Exception as ex: print(f'  {sv}: download failed ({ex})'); continue
    try: h = fits.open(fn)[0]; im = np.asarray(h.data, float); w = WCS(h.header)
    except Exception as ex: print(f'  {sv}: not an image'); continue
    if not np.isfinite(im).any() or np.nanstd(im) == 0: print(f'  {sv}: no coverage'); continue
    yy, xx = np.indices(im.shape); x0, y0 = [float(z) for z in w.all_world2pix(a.ra, a.dec, 0)]; sc = 0.05 * 3600 / 180                      # arcsec per pixel; East is -x
    px = lambda E, N: (x0 - E / sc, y0 + N / sc); sky = np.nanmedian(im); sig = 1.4826 * np.nanmedian(np.abs(im - sky))
    best = None                                                                          # the star's own image along its proper-motion track dates the plate
    for yr in np.arange(1945, 2006, 0.5):
        cx, cy = px(*off(yr)); m = np.hypot(xx - cx, yy - cy) < 2.0 / sc; fl = np.nansum(im[m] - sky)
        if best is None or fl > best[0]: best = (fl, yr)
    yr0 = best[1]; cx, cy = px(*off(yr0)); m = np.hypot(xx - cx, yy - cy) < 4.0 / sc; wgt = np.clip(im[m] - sky, 0, None)
    if wgt.sum() <= 0: print(f'  {sv}: star not found'); continue
    sx, sy = np.sum(xx[m] * wgt) / wgt.sum(), np.sum(yy[m] * wgt) / wgt.sum(); dE, dN = -(sx - x0) * sc, (sy - y0) * sc; yr_fit = 2000 + (dE * pmE + dN * pmN) / pm ** 2   # epoch from the measured displacement
    nx, ny = px(*off(ep_tess)); sep = np.hypot(sx - nx, sy - ny) * sc; rap = 2.5 / sc
    ap = lambda cx_, cy_: np.nansum((im - sky)[np.hypot(xx - cx_, yy - cy_) < rap])
    rng = np.random.default_rng(0); rnd = np.array([ap(rng.uniform(15, im.shape[1] - 15), rng.uniform(15, im.shape[0] - 15)) for _ in range(400)]); rs = 1.4826 * np.median(np.abs(rnd - np.median(rnd))); f_now = ap(nx, ny) - np.median(rnd); f_star = ap(sx, sy) - np.median(rnd)
    # zero point from field stars (positions moved to the plate epoch); photographic plates are non-linear, so this is approximate
    zs = []
    for _, r in g[(g.ID != TIC) & np.isfinite(g.GAIAmag) & (g.GAIAmag > 14) & (g.GAIAmag < 19)].iterrows():
        E = (r.ra - a.ra) * np.cos(np.radians(a.dec)) * 3600 + (0 if not np.isfinite(r.pmRA) else r.pmRA / 1000) * (yr_fit - 2000); N = (r.dec - a.dec) * 3600 + (0 if not np.isfinite(r.pmDEC) else r.pmDEC / 1000) * (yr_fit - 2000); qx, qy = px(E, N)
        if 12 < qx < im.shape[1] - 12 and 12 < qy < im.shape[0] - 12:
            fl = ap(qx, qy) - np.median(rnd)
            if fl > 8 * rs: zs.append(r.GAIAmag + 2.5 * np.log10(fl))
    zp = np.median(zs) if len(zs) >= 3 else (a.GAIAmag + 2.5 * np.log10(f_star) if f_star > 0 else np.nan); lim = zp - 2.5 * np.log10(3 * rs) if np.isfinite(zp) else np.nan
    usable = sep > 5.0
    print(f'  {sv}: star found {np.hypot(dE, dN):.1f}" from its 2000 position -> image date about {yr_fit:.0f}; it is {sep:.1f}" from its TESS-epoch position there. Flux at the TESS-epoch position: {f_now / rs:+.1f} sigma; 3-sigma limit G ~ {lim:.1f} ({len(zs)} field stars for the zero point)' + ('' if usable else '   [too close to the star image to be useful]'))
    res.append(dict(survey=sv, year=float(yr_fit), sep=float(sep), nsig=float(f_now / rs), lim=float(lim), usable=bool(usable))); panels.append((sv, im, sky, sig, (sx, sy), (nx, ny), yr_fit, sc, sep))
ok = [r for r in res if r['usable'] and np.isfinite(r['lim'])]
if ok:
    clean = [r for r in ok if r['nsig'] < 3]; b = max(clean or ok, key=lambda r: (r['sep'] > 8, r['lim'])); need = a.GAIAmag + dmax
    print(f'  deepest usable image: {b["survey"]} ({b["year"]:.0f}), nothing at the present position down to G ~ {b["lim"]:.1f} (flux {b["nsig"]:+.1f} sigma); a mimicking star would need G < {need:.1f}  ->  ' + ('a background star able to mimic the signal is excluded by this image' if b['lim'] > need and b['nsig'] < 3 else 'flux at the position in every usable image (often the halo of the star itself): look at the figure' if b['nsig'] >= 3 else f'background stars brighter than G ~ {b["lim"]:.1f} are excluded; fainter ones (needing an eclipse deeper than {min(100, 100 * depth * 10 ** (0.4 * (b["lim"] - a.GAIAmag))):.0f} %) are not'))
json.dump(dict(tic=TIC, G=float(a.GAIAmag), need_G=float(a.GAIAmag + dmax), ep_tess=float(ep_tess), best=(b if ok else None), images=res), open(f'{R}/deep2/vet/{TIC}_{K}.archive.json', 'w'))
if panels:
    fig, ax = plt.subplots(1, len(panels), figsize=(3.3 * len(panels), 3.7), squeeze=False)
    for axx, (sv, im, sky, sig, (sx, sy), (nx, ny), yr, sc, sep) in zip(ax[0], panels):
        hw = 30 / sc; axx.imshow(im, origin='lower', cmap='Greys', vmin=sky - 2 * sig, vmax=sky + 25 * sig); axx.set_xlim(nx - hw, nx + hw); axx.set_ylim(ny - hw, ny + hw)
        axx.add_patch(plt.Circle((nx, ny), 4 / sc, fill=False, color='#c0392b', lw=1.6)); axx.plot([sx], [sy], '+', color='#2a6fdb', ms=14, mew=1.6)
        axx.set_title(f'{sv}, about {yr:.0f}', fontsize=10, loc='left'); axx.set_xticks([]); axx.set_yticks([]); axx.set_xlabel(f'star then (+) is {sep:.1f}″ from where it is now (circle)', fontsize=8)
    fig.suptitle(f'TIC {TIC}: is anything behind the star?  1′ × 1′, North up, East left', x=0.01, ha='left', fontsize=11); plt.tight_layout(rect=(0, 0, 1, 0.93)); plt.savefig(f'{D}/archive_TIC{TIC}.png', dpi=110); print('  figure', f'{D}/archive_TIC{TIC}.png')
