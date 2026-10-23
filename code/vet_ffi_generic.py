"""Gaia neighbours, stacked on-target PSF test with TESS-SPOC FFI TPFs, and TRICERATOPS for TIC 198079261."""
import os, glob, io, json, numpy as np, pandas as pd, urllib.request, urllib.parse, warnings, lightkurve as lk
from astropy.coordinates import SkyCoord
warnings.filterwarnings('ignore')
import sys
TIC = int(sys.argv[1]); t16 = f'{TIC:016d}'; D = f'/home/sarang/tess_search/deep/{TIC}'; O = f'/home/sarang/tess_search/cand_{TIC}'; os.makedirs(O, exist_ok=True)
e = np.load(f'{D}/ephem_ffi.npz'); P, t0, dur = float(e['P']), float(e['t0']), float(e['dur']); depth = float(sys.argv[2]) * 1e-6
ra0, dec0 = float(sys.argv[3]), float(sys.argv[4])
# --- Gaia neighbours
q = f"SELECT source_id, ra, dec, phot_g_mean_mag, parallax, ruwe, non_single_star FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT(ra,dec), CIRCLE({ra0},{dec0},120./3600)) ORDER BY phot_g_mean_mag"
g = pd.read_csv(io.StringIO(urllib.request.urlopen(urllib.request.Request('https://gea.esac.esa.int/tap-server/tap/sync', data=urllib.parse.urlencode(dict(REQUEST='doQuery', LANG='ADQL', FORMAT='csv', QUERY=q)).encode()), timeout=120).read().decode()))
ti = (np.hypot((g.ra - ra0) * np.cos(np.radians(dec0)), g.dec - dec0)).argmin(); tg = g.iloc[ti]
g['sep'] = SkyCoord(g.ra, g.dec, unit='deg').separation(SkyCoord(tg.ra, tg.dec, unit='deg')).arcsec
g['dG'] = g.phot_g_mean_mag - tg.phot_g_mean_mag; g['need'] = depth * (1 + 10 ** (-.4 * g.dG)) / 10 ** (-.4 * g.dG)
cap = g[(g.sep > .5) & (g.need < 1)]
print(f'target G={tg.phot_g_mean_mag:.2f} plx={tg.parallax:.2f} RUWE={tg.ruwe:.2f} NSS={tg.non_single_star}; neighbours <120" able to mimic: {len(cap)}')
print(cap[['source_id', 'phot_g_mean_mag', 'dG', 'sep', 'need', 'parallax']].round(3).to_string(index=False)); g.to_csv(f'{O}/gaia_neighbors_120arcsec.csv', index=False)
# --- TESS-SPOC FFI TPFs: stacked fit of difference images with target PSF
aps, secs, amp = [], [], []
for s in [int(x) for x in sys.argv[5].split(',')]:
    n = f'hlsp_tess-spoc_tess_phot_{t16}-s{s:04d}_tess_v1_tp.fits'; p = f'{D}/{n}'
    for _try in range(5):
        if os.path.exists(p) and os.path.getsize(p) > 100000: break
        try: urllib.request.urlretrieve(f'https://archive.stsci.edu/hlsps/tess-spoc/s{s:04d}/target/{t16[0:4]}/{t16[4:8]}/{t16[8:12]}/{t16[12:16]}/{n}', p)
        except Exception as ex: print('tpf retry', s, ex)
    if not os.path.exists(p): continue
    tpf = lk.read(p, quality_bitmask='default'); t = tpf.time.value; fl = tpf.flux.value; ok = np.all(np.isfinite(fl.reshape(len(t), -1)), 1); t, fl = t[ok], fl[ok]
    ph = (((t - t0) / P + .5) % 1 - .5) * P; nn = np.floor((t - t0) / P + .5).astype(int)
    intr, near = np.abs(ph) < dur / 2, (np.abs(ph) > dur) & (np.abs(ph) < 4 * dur)
    diffs = [np.nanmedian(fl[near & (nn == k)], 0) - np.nanmean(fl[intr & (nn == k)], 0) for k in np.unique(nn[intr]) if (intr & (nn == k)).sum() > 3 and (near & (nn == k)).sum() > 10]
    if len(diffs) < int(os.environ.get('MIN_TR', 3)): continue
    diff = np.nanmean(diffs, 0); oot = np.nanmedian(fl[near], 0)
    err = (np.nanstd(diffs, 0) / np.sqrt(len(diffs))) if len(diffs) > 2 else np.nanstd(fl[near], 0) * np.sqrt(1 / max(intr.sum(), 1) + 1 / max(near.sum(), 1))
    x0, y0 = tpf.wcs.world_to_pixel(SkyCoord(tg.ra, tg.dec, unit='deg')); yy, xx = np.indices(diff.shape)
    m = np.hypot(xx - x0, yy - y0) < 2.5; w = 1 / err[m] ** 2
    a = np.sum(w * diff[m] * oot[m]) / np.sum(w * oot[m] ** 2); ae = np.sqrt(1 / np.sum(w * oot[m] ** 2))
    # also: best single-pixel location of the dip vs target
    sn = diff / err; dy, dx = np.unravel_index(np.nanargmax(sn), sn.shape)
    print(f'  S{s}: {len(diffs)} transits; on-target PSF amplitude {a*1e6:.0f} ± {ae*1e6:.0f} ppm; peak S/N pixel ({dx},{dy}) vs target ({x0:.1f},{y0:.1f}), peak S/N {sn[dy, dx]:.1f}')
    amp.append((a, ae)); secs.append(s); r_, c_ = np.where(tpf.pipeline_mask); aps.append(np.array([[tpf.column + cc, tpf.row + rr] for rr, cc in zip(r_, c_)]))
A = np.array(amp); wm = np.sum(A[:, 0] / A[:, 1] ** 2) / np.sum(1 / A[:, 1] ** 2); we = 1 / np.sqrt(np.sum(1 / A[:, 1] ** 2))
print(f'combined on-target amplitude: {wm*1e6:.0f} ± {we*1e6:.0f} ppm  ({wm/we:.1f} sigma); light-curve depth {depth*1e6:.0f}')
# --- TRICERATOPS
import triceratops.triceratops as tr
tg_ = tr.target(ID=TIC, sectors=np.array(secs)); tg_.calc_depths(tdepth=depth, all_ap_pixels=aps)
print(tg_.stars[['ID', 'Tmag', 'sep (arcsec)', 'fluxratio', 'tdepth']].head(8).to_string())
T, F = [], []
for fn in sorted(glob.glob(f'{D}/*_tess_v1_lc.fits')):
    lc = lk.read(fn, quality_bitmask='default').remove_nans().normalize(); cad = np.nanmedian(np.diff(lc.time.value))
    if cad > float(os.environ.get("MAXCAD_MIN", 6)) / 1440: continue                # finest-cadence sectors only for the shape
    ph = ((lc.time.value - t0) / P + .5) % 1 - .5
    fl = lc.flatten(window_length=int(round(.75 / cad)) | 1, mask=np.abs(ph * P) < 2 * dur).remove_outliers(sigma_upper=4, sigma_lower=50)
    ph = (((fl.time.value - t0) / P + .5) % 1 - .5) * P; mm = np.abs(ph) < 0.15
    T.append(ph[mm]); F.append(np.ma.filled(np.ma.asarray(fl.flux.value[mm], dtype=float), np.nan))
t = np.concatenate(T); f = np.concatenate(F); ok = np.isfinite(f); t, f = t[ok], f[ok]; o = np.argsort(t); t, f = t[o], f[o]
b = np.arange(t.min(), t.max(), 6 / 1440); idx = np.digitize(t, b); u = np.unique(idx)
tb = np.array([t[idx == k].mean() for k in u]); fb = np.array([f[idx == k].mean() for k in u]); eb = np.median([f[idx == k].std() / np.sqrt(np.sum(idx == k)) for k in u])
res = []
for run in range(5):
    tg_.calc_probs(time=tb, flux_0=fb, flux_err_0=eb, P_orb=P, N=300000, parallel=True); res.append((tg_.FPP, tg_.NFPP))
    print(f'run {run}: FPP = {tg_.FPP:.4f}  NFPP = {tg_.NFPP:.4f}', flush=True)
print(tg_.probs[['ID', 'scenario', 'prob']].sort_values('prob', ascending=False).head(10).to_string())
R = np.array(res); print(f'\nTRICERATOPS FPP = {R[:,0].mean():.3f} ± {R[:,0].std():.3f}, NFPP max = {R[:,1].max():.4f}')
json.dump(dict(FPP=list(R[:, 0]), NFPP=list(R[:, 1]), on_target_amp=[wm, we]), open(f'{O}/vet_summary.json', 'w'))
tg_.probs.to_csv(f'{O}/triceratops_probs.csv', index=False)
