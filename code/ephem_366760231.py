import glob, json, numpy as np, lightkurve as lk, warnings
from scipy.optimize import minimize_scalar
from astropy.time import Time
warnings.filterwarnings('ignore')
D = '/home/sarang/tess_search/deep/366760231'; O = '/home/sarang/tess_search/cand_366760231'
fit = json.load(open(f'{O}/fit_results.json')); P0, T0 = fit['period_d'][0], fit['T0_BTJD'][0]; dur = fit['T14_h'][0] / 24; dep = fit['depth_ppm'][0] * 1e-6
files = sorted(glob.glob(f'{D}/*s_lc.fits')) + sorted(glob.glob(f'{D}/hlsp_tess-spoc*.fits'))
tt, te, ep = [], [], []
for fn in files:
    lc = lk.read(fn, quality_bitmask='default').remove_nans().normalize(); cad = np.nanmedian(np.diff(lc.time.value)); s = lc.meta['SECTOR']
    fl = lc.flatten(window_length=int(round(.75 / cad)) | 1).remove_outliers(sigma_upper=4, sigma_lower=50)
    t = np.asarray(fl.time.value, float); f = np.ma.filled(np.ma.asarray(fl.flux.value, dtype=float), np.nan); ok = np.isfinite(f); t, f = t[ok], f[ok]
    n0 = np.round((np.median(t) - T0) / P0); tc0 = T0 + n0 * P0; sig = 1.4826 * np.median(np.abs(f - 1))
    def chi(dt):   # box model with fixed depth/duration shifted by dt, phase-folded within the sector
        ph = ((t - tc0 - dt) / P0 + .5) % 1 - .5; m = np.abs(ph * P0) < dur / 2
        return np.sum((f - np.where(m, 1 - dep, 1)) ** 2) / sig ** 2
    grid = np.linspace(-0.08, 0.08, 321); c = np.array([chi(g) for g in grid]); g0 = grid[np.argmin(c)]
    ok1 = grid[c < c.min() + 1]; err = max((ok1.max() - ok1.min()) / 2, cad)
    tt.append(tc0 + g0); te.append(err); ep.append(n0); print(f'S{s}: mid-time {tc0+g0:.4f} ± {err*1440:.0f} min (epoch {n0:.0f})')
tt, te, ep = map(np.array, (tt, te, ep)); A = np.vstack([np.ones_like(ep), ep]).T; W = np.diag(1 / te ** 2)
cov = np.linalg.inv(A.T @ W @ A); (Tref, P) = cov @ A.T @ W @ tt; sT, sP = np.sqrt(np.diag(cov))
res = tt - (Tref + P * ep); chi2 = np.sum(res ** 2 / te ** 2)
print(f'linear ephemeris: P = {P:.6f} ± {sP:.6f} d, T0 = {Tref + 2457000:.5f} ± {sT:.5f} BJD_TDB (epoch 0 = fit reference); chi2/dof = {chi2:.1f}/{len(tt)-2}')
json.dump(dict(P=P, sP=sP, T0_BJD=Tref + 2457000, sT0=sT, cov=cov.tolist()), open(f'{O}/ephemeris.json', 'w'), indent=1)
now = Time('2026-10-22').jd - 2457000; n1 = int(np.ceil((now - Tref) / P)); lines = []
for n in range(n1, n1 + 14):
    tc = Tref + n * P; s = np.sqrt(cov[0, 0] + n ** 2 * cov[1, 1] + 2 * n * cov[0, 1])
    lines.append(f'| {n} | {tc + 2457000:.4f} | {Time(tc + 2457000, format="jd").iso[:16]} UTC | ±{s*1440:.0f} min |')
open(f'{O}/upcoming_transits.md', 'w').write('| epoch | BJD_TDB mid-transit | date (approx.) | 1σ |\n|---|---|---|---|\n' + '\n'.join(lines) + '\n'); print('\n'.join(lines[:3]))
