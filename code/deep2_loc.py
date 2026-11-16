"""Where on the sky is the dip?  Joint fit over every pixel file already in deep2/pix/{TIC}/ (sectors with even one transit count):
each sector's out-of-transit image is fitted with all catalogue stars in the stamp (fixed TESS-band flux ratios; common shift, PSF width and
scale free); the in-transit difference images of all sectors are then fitted together with ONE point source of common brightness change
placed at a trial sky offset from the target.  Output: chi-square map over +-50", best position, and how strongly the target and each
neighbour are (dis)favoured.   Usage: deep2_loc.py TIC K [other_signal_K ...]"""
import sys, os, json, glob, numpy as np, pandas as pd, warnings, lightkurve as lk
from astropy.coordinates import SkyCoord
from scipy.ndimage import median_filter
from scipy.optimize import least_squares
from scipy.special import erf
from scipy.stats import chi2 as chi2d, norm
warnings.filterwarnings('ignore')
R = '/home/sarang/tess_search'; TIC, K = int(sys.argv[1]), int(sys.argv[2]); OTH = [int(x) for x in sys.argv[3:]]
v = json.load(open(f'{R}/deep2/vet/{TIC}_{K}.json')); P, t0, dur = v['P'], v['T0'], v['T14_h'] / 24
oth = [json.load(open(f'{R}/deep2/vet/{TIC}_{k}.json')) for k in OTH]
D = f'{R}/deep2/pix/{TIC}'; g = pd.read_csv(f'{D}/tic_120arcsec.csv').drop_duplicates('ID'); g = g[np.isfinite(g.Tmag)].reset_index(drop=True); A = g[g.ID == TIC].iloc[0]
g['dE0'] = (g.ra - A.ra) * np.cos(np.radians(A.dec)) * 3600; g['dN0'] = (g.dec - A.dec) * 3600; g['rel'] = 10 ** (-.4 * (g.Tmag - A.Tmag))
g['pmE'] = pd.to_numeric(g.pmRA, errors='coerce').fillna(0) / 1000; g['pmN'] = pd.to_numeric(g.pmDEC, errors='coerce').fillna(0) / 1000; pmA = (float(g.pmE[g.ID == TIC].iloc[0]), float(g.pmN[g.ID == TIC].iloc[0]))
print(f'TIC {TIC}: proper motion {pmA[0]*1000:+.0f}, {pmA[1]*1000:+.0f} mas/yr (catalogue positions are for epoch 2000 and are moved to each sector\'s date)')
def at_epoch(yr):                                            # offsets of every catalogue star from the TARGET at that date, and the target's own displacement since 2000
    dt = yr - 2000.0; return g.dE0.values + (g.pmE.values - pmA[0]) * dt, g.dN0.values + (g.pmN.values - pmA[1]) * dt, pmA[0] * dt, pmA[1] * dt

def psf(xx, yy, x, y, s):                                   # pixel-integrated circular Gaussian, unit flux
    q = np.sqrt(2) * s
    return .25 * (erf((xx + .5 - x) / q) - erf((xx - .5 - x) / q)) * (erf((yy + .5 - y) / q) - erf((yy - .5 - y) / q))

S = []; ntr_tot = 0; yrs = []
for p in sorted(glob.glob(f'{D}/*_tp.fits')):
    try: tpf = lk.read(p, quality_bitmask='default')
    except Exception as ex: print('  read', os.path.basename(p), ex); continue
    t = tpf.time.value; fl = np.asarray(tpf.flux.value, float); ok = np.all(np.isfinite(fl.reshape(len(t), -1)), 1); t, fl = t[ok], fl[ok]
    if len(t) < 200: continue
    cad = np.nanmedian(np.diff(t)); kk = int(round(0.75 / cad)) | 1; med = np.nanmedian(fl, 0); fl = fl - median_filter(fl, size=(kk, 1, 1), mode='nearest') + med
    keep = np.ones(len(t), bool)
    for o in oth: keep &= np.abs(((t - o['T0']) / o['P'] + .5) % 1 - .5) * o['P'] > 0.75 * o['T14_h'] / 24
    t, fl = t[keep], fl[keep]
    ph = (((t - t0) / P + .5) % 1 - .5) * P; nn = np.floor((t - t0) / P + .5).astype(int)
    intr, near = np.abs(ph) < 0.4 * dur, (np.abs(ph) > 0.7 * dur) & (np.abs(ph) < 3.5 * dur)
    good = [k for k in np.unique(nn[intr]) if (intr & (nn == k)).sum() >= max(1, int(0.4 * 0.8 * dur / cad)) and (near & (nn == k)).sum() >= 6]
    if not good: continue
    sel = np.isin(nn, good); i_, n_ = intr & sel, near & sel
    diff = np.nanmean(fl[n_], 0) - np.nanmean(fl[i_], 0); sd = 1.4826 * np.nanmedian(np.abs(fl[n_] - np.nanmedian(fl[n_], 0)), 0); err = sd * np.sqrt(1 / i_.sum() + 1 / n_.sum()) + 1e-9
    lo = 4.5 * dur / P; nulls = []                             # empirical per-pixel noise: the same difference image at other orbital phases
    for sh in (np.linspace(lo, 1 - lo, 31) if lo < 0.4 else []):
        ph2 = (((t - t0) / P + .5 + sh) % 1 - .5) * P; n2 = np.floor((t - t0) / P + .5 + sh).astype(int); i2, r2 = np.abs(ph2) < 0.4 * dur, (np.abs(ph2) > 0.7 * dur) & (np.abs(ph2) < 3.5 * dur)
        g2 = [k for k in np.unique(n2[i2]) if (i2 & (n2 == k)).sum() >= max(1, int(0.4 * 0.8 * dur / cad)) and (r2 & (n2 == k)).sum() >= 6]
        if len(g2) >= max(1, len(good) // 2): s2 = np.isin(n2, g2); nulls.append((np.nanmean(fl[r2 & s2], 0) - np.nanmean(fl[i2 & s2], 0)) * np.sqrt(len(g2) / len(good)))
    if len(nulls) >= 10: e2 = 1.4826 * np.nanmedian(np.abs(np.array(nulls) - np.nanmedian(nulls, 0)), 0); err = np.maximum(err, e2)
    oot = np.nanmedian(fl[n_], 0); yy, xx = np.indices(oot.shape)
    yr = 2000.0 + (np.nanmedian(t) + 2457000 - 2451545.0) / 365.25; dEs, dNs, mE, mN = at_epoch(yr); ra_t = A.ra + mE / 3600 / np.cos(np.radians(A.dec)); de_t = A.dec + mN / 3600; yrs.append(yr)
    xa, ya = [float(z) for z in tpf.wcs.world_to_pixel(SkyCoord(ra_t, de_t, unit='deg'))]
    cE = SkyCoord(ra_t + 10 / 3600 / np.cos(np.radians(A.dec)), de_t, unit='deg'); cN = SkyCoord(ra_t, de_t + 10 / 3600, unit='deg')
    JE = (np.array([float(z) for z in tpf.wcs.world_to_pixel(cE)]) - [xa, ya]) / 10; JN = (np.array([float(z) for z in tpf.wcs.world_to_pixel(cN)]) - [xa, ya]) / 10     # pixels per arcsec E / N
    sx = xa + dEs * JE[0] + dNs * JN[0]; sy = ya + dEs * JE[1] + dNs * JN[1]
    ins = (sx > -1) & (sx < oot.shape[1]) & (sy > -1) & (sy < oot.shape[0]); m = np.hypot(xx - xa, yy - ya) < 3.5; bg = np.nanpercentile(oot, 10)
    def model(q): return q[3] * sum(r * psf(xx, yy, x + q[0], y + q[1], q[2]) for r, x, y in zip(g.rel.values[ins], sx[ins], sy[ins])) + q[4]
    mr = np.hypot(xx - xa, yy - ya) < float(os.environ.get('RREG', 2.2))          # register the scene on the pixels around the target itself
    fo = least_squares(lambda q: (model(q) - oot)[mr], [0, 0, 0.8, (oot - bg)[mr].sum(), bg], bounds=([-1.5, -1.5, 0.4, 0, -np.inf], [1.5, 1.5, 2.0, np.inf, np.inf])).x
    nnull = len(nulls)
    S.append(dict(dimg=diff, px=sx[ins] + fo[0], py=sy[ins] + fo[1], ptm=g.Tmag.values[ins], pid=g.ID.values[ins], nnull=nnull, sec=int(tpf.sector), d=diff[m], w=1 / err[m] ** 2, xx=xx[m], yy=yy[m], x0=xa + fo[0], y0=ya + fo[1], sg=fo[2], F=fo[3], JE=JE, JN=JN, ntr=len(good))); ntr_tot += len(good)
    print(f'  S{int(tpf.sector)}: {len(good)} transit(s), {i_.sum()} in-transit cadences; image fit shift ({fo[0]*21:+.1f}",{fo[1]*21:+.1f}"), PSF sigma {fo[2]:.2f} px', flush=True)
if not S: print('no usable transits in the local pixel files'); sys.exit()

def sec_chi(s, E=0.0, N=0.0):                               # one sector alone: chi2 of the best point source at (E, N) plus a constant
    pm = s['F'] * psf(s['xx'], s['yy'], s['x0'] + E * s['JE'][0] + N * s['JN'][0], s['y0'] + E * s['JE'][1] + N * s['JN'][1], s['sg']); w = s['w']
    dd = s['d'] - np.sum(w * s['d']) / w.sum(); pp = pm - np.sum(w * pm) / w.sum(); a_ = np.sum(w * pp * dd) / np.sum(w * pp * pp); return (np.sum(w * dd * dd) - a_ * a_ * np.sum(w * pp * pp)) / (len(dd) - 2)
for s in S: s['rc'] = min(sec_chi(s, e, n) for e in (-20, 0, 20) for n in (-20, 0, 20))
medrc = np.median([s['rc'] for s in S]); bad = [s for s in S if s['rc'] > max(4.0, 4 * medrc) or s['nnull'] < 10 and s['rc'] > 4]
for s in bad: print(f"  dropping S{s['sec']}: residuals {s['rc']:.1f}x the noise estimate ({s['nnull']} null phases) - pixel data not usable")
badid = {id(s) for s in bad}; S = [s for s in S if id(s) not in badid]; ntr_tot = sum(s['ntr'] for s in S)
for s in S: s['w'] = s['w'] / max(s['rc'], 1.0)             # per-sector error inflation
if not S: print('no usable sectors'); sys.exit()

def chi(E, N):                                              # common amplitude (target-flux units), free constant per sector
    num = den = c0 = 0.0
    for s in S:
        pm = s['F'] * psf(s['xx'], s['yy'], s['x0'] + E * s['JE'][0] + N * s['JN'][0], s['y0'] + E * s['JE'][1] + N * s['JN'][1], s['sg']); w = s['w']
        dd = s['d'] - np.sum(w * s['d']) / w.sum(); pp = pm - np.sum(w * pm) / w.sum(); num += np.sum(w * pp * dd); den += np.sum(w * pp * pp); c0 += np.sum(w * dd * dd)
    a = num / den; return c0 - a * a * den, a, 1 / np.sqrt(den), c0
ax = np.arange(-50, 51, 2.0); M = np.array([[chi(e, n)[0] for e in ax] for n in ax]); j, i = np.unravel_index(np.argmin(M), M.shape)
fr = least_squares(lambda q: [np.sqrt(max(chi(q[0], q[1])[0], 0))], [ax[i], ax[j]], diff_step=0.02).x; cb, ab, eb, c0 = chi(*fr)
ndof = sum(len(s['d']) - 1 for s in S) - 3; sc = max(cb / ndof, 1.0)                      # inflate errors if the best fit is poor
SYS = 4.0                                                   # arcsec: registration systematics (0.2 pixel) seen when the alignment region is changed
def nsig(E, N):                                             # best chi2 anywhere within SYS of the star's position
    dc = (min(chi(E + r * np.cos(a_), N + r * np.sin(a_))[0] for r, a_ in [(0, 0)] + [(SYS / 2, x) for x in np.arange(0, 6.2, 1.571)] + [(SYS, x) for x in np.arange(0, 6.2, 0.785)]) - cb) / sc
    return norm.isf(chi2d.sf(max(dc, 0), 2) / 2) if dc > 0 else 0.0
reg = M - cb < 2.30 * sc; eE = (ax[reg.any(0)].max() - ax[reg.any(0)].min()) / 2 if reg.any() else np.nan; eN = (ax[reg.any(1)].max() - ax[reg.any(1)].min()) / 2 if reg.any() else np.nan
ct, at, et, _ = chi(0, 0); det = np.sqrt(max(c0 - cb, 0) / sc)
print(f'TIC {TIC} k={K}: {len(S)} sectors, {ntr_tot} transits; detection of a point-source dip in the pixels: {det:.1f} sigma (reduced chi2 {cb/ndof:.2f})')
print(f'  best position: {fr[0]:+.1f}" East, {fr[1]:+.1f}" North of the target (1-sigma region about ±{eE:.0f}" x ±{eN:.0f}");  target is {nsig(0, 0):.1f} sigma from the best position (positions closer than {SYS:.0f}" are not distinguished)')
print(f'  depth if on the target: {at*1e6:.0f} ± {et*np.sqrt(sc)*1e6:.0f} ppm of the target flux (light curve: {v["trap_depth_ppm"]:.0f} ppm)')
g['dE'], g['dN'] = at_epoch(np.mean(yrs))[:2]; g['nsig'] = [nsig(e, n) for e, n in zip(g.dE, g.dN)]; g['sep'] = np.hypot(g.dE, g.dN); g['need'] = v['trap_depth_ppm'] * 1e-6 * (1 + g.rel) / g.rel
c = g[(g.sep < 60)].sort_values('nsig'); print('  catalogue stars, by how well the dip position matches them (need = eclipse depth needed on that star):')
print(c[['ID', 'Tmag', 'sep', 'dE', 'dN', 'need', 'nsig']].head(8).round(2).to_string(index=False))
al = g[(g.ID != TIC) & (g.need < 1) & (g.sep < 60)]; best_alt = al.sort_values('nsig').iloc[0] if len(al) else None
if best_alt is not None: print(f'  least-excluded other star able to mimic the signal: TIC {int(best_alt.ID)} (T={best_alt.Tmag:.1f}, {best_alt.sep:.1f}" away) at {best_alt.nsig:.1f} sigma')
if os.environ.get('PERSEC'):                                 # each sector on its own: where does it put the dip?
    Sall = S; pos = []
    for s1 in Sall:
        S = [s1]; m1 = np.array([[chi(e, n)[0] for e in np.arange(-24, 25, 2.0)] for n in np.arange(-24, 25, 2.0)]); jj, ii = np.unravel_index(np.argmin(m1), m1.shape); c00 = chi(0, 0)
        print(f"    S{s1['sec']}: {s1['ntr']} tr; best ({np.arange(-24, 25, 2.0)[ii]:+.0f}\",{np.arange(-24, 25, 2.0)[jj]:+.0f}\"); dchi2(target-best) {c00[0]-m1.min():.1f}; depth on target {c00[1]*1e6:.0f}±{c00[2]*1e6:.0f} ppm"); pos.append((np.arange(-24, 25, 2.0)[ii], np.arange(-24, 25, 2.0)[jj]))
    S = Sall; pos = np.array(pos); print(f'    per-sector positions: mean ({pos[:,0].mean():+.1f}",{pos[:,1].mean():+.1f}"), scatter ({pos[:,0].std():.1f}",{pos[:,1].std():.1f}") -> error of mean ({pos[:,0].std()/np.sqrt(len(pos)):.1f}",{pos[:,1].std()/np.sqrt(len(pos)):.1f}")')
Sall = S; best1 = None                                      # the single sector that shows the dip best, for the figure
for s1 in Sall:
    S = [s1]; c1 = chi(fr[0], fr[1]); z1 = np.sqrt(max(c1[3] - c1[0], 0))
    if best1 is None or z1 > best1[0]: best1 = (z1, s1)
S = Sall; b1 = best1[1]; np.savez(f'{D}/locbest_{K}.npz', diff=b1['dimg'], x0=b1['x0'], y0=b1['y0'], sec=b1['sec'], px=b1['px'], py=b1['py'], ptm=b1['ptm'], pid=b1['pid'], z=best1[0])
np.save(f'{D}/locmap_{K}.npy', (M - cb) / sc)
json.dump(dict(tic=TIC, k=K, nsec=len(S), ntr=ntr_tot, det_sigma=det, best_E=fr[0], best_N=fr[1], err_E=eE, err_N=eN, target_nsig=nsig(0, 0), depth_on_target_ppm=at * 1e6, depth_on_target_err=et * np.sqrt(sc) * 1e6, red_chi2=cb / ndof,
               alt_tic=None if best_alt is None else int(best_alt.ID), alt_nsig=None if best_alt is None else float(best_alt.nsig), alt_sep=None if best_alt is None else float(best_alt.sep), alt_Tmag=None if best_alt is None else float(best_alt.Tmag),
               sectors=[s['sec'] for s in S]), open(f'{R}/deep2/vet/{TIC}_{K}.loc.json', 'w'), default=float)
