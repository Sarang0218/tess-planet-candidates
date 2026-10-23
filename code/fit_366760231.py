"""MCMC transit fit for TIC 366760231 on unbinned 2-min PDCSAP data (sectors 19, 59, 60).
Free: P, T0 (reference epoch near the data centroid), k = Rp/R*, b, rho* (prior from TIC M, R), ln(jitter), per-sector offsets absorbed by flattening.
Fixed quadratic limb darkening appropriate for Teff~3200 K, logg~5 in the TESS band."""
import glob, numpy as np, lightkurve as lk, batman, emcee, warnings, json
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')
OUT = '/home/sarang/tess_search/cand_366760231'
e = np.load('/home/sarang/tess_search/deep/366760231/ephem_all.npz'); P0, T00, D0 = float(e['P']), float(e['t0']), max(float(e['dur']), 0.04)
MS, RS, MS_E, RS_E, TEFF, DIST = 0.3787, 0.3887, 0.0203, 0.0116, 3456., 55.31   # TIC v8.2; errors ~TIC typical for M dwarfs
U = (0.20, 0.40)
RHO0 = 1.41 * MS / RS ** 3; RHO_E = RHO0 * np.hypot(MS_E / MS, 3 * RS_E / RS)
import os
TAG = '_freerho' if os.environ.get('FREE_RHO') else ''
if TAG: RHO_E = 1e3

T, F, E, S = [], [], [], []
for fn in sorted(glob.glob('/home/sarang/tess_search/deep/366760231/*s_lc.fits')):
    lc = lk.read(fn, quality_bitmask='default').remove_nans().normalize(); sec = lc.meta['SECTOR']
    ph = ((lc.time.value - T00) / P0 + .5) % 1 - .5
    mask = np.abs(ph * P0) < 1.5 * D0
    cad = np.nanmedian(np.diff(lc.time.value))
    fl = lc.flatten(window_length=int(round(.75 / cad)) | 1, mask=mask).remove_outliers(sigma_upper=4, sigma_lower=50)
    t = np.asarray(fl.time.value, float); f = np.ma.filled(np.ma.asarray(fl.flux.value, dtype=float), np.nan); ef = np.ma.filled(np.ma.asarray(fl.flux_err.value, dtype=float), np.nan)
    ph = ((t - T00) / P0 + .5) % 1 - .5; near = (np.abs(ph * P0) < 4 * D0) & np.isfinite(f) & np.isfinite(ef)
    T.append(t[near]); F.append(f[near]); E.append(ef[near]); S.append(np.full(near.sum(), sec))
t, f, ef, sec = map(np.concatenate, (T, F, E, S))
n_all = np.round((t - T00) / P0); nref = np.round(np.median(n_all)); Tref0 = T00 + nref * P0
print(f'{len(t)} points near transit; reference epoch index {nref:.0f} -> T0 guess {Tref0:.4f} BTJD; rho* prior {RHO0:.1f} ± {RHO_E:.1f} g/cc')

params = batman.TransitParams(); params.t0 = 0.; params.per = P0; params.rp = .05; params.a = 30.; params.inc = 89.; params.ecc = 0.; params.w = 90.; params.u = list(U); params.limb_dark = 'quadratic'
def model(theta, tt):
    P, T0, k, b, rho, lns = theta
    aR = (rho * 6.674e-8 * (P * 86400) ** 2 / (3 * np.pi)) ** (1 / 3)     # rho in g/cc
    inc = np.degrees(np.arccos(b / aR))
    params.per, params.t0, params.rp, params.a, params.inc = P, T0, k, aR, inc
    return batman.TransitModel(params, tt, supersample_factor=3, exp_time=2 / 1440).light_curve(params), aR, inc
def lnprob(theta):
    P, T0, k, b, rho, lns = theta
    if not (0 < k < .5 and 0 <= b < 1 + k and .5 < rho < 200 and -12 < lns < -4 and abs(P - P0) < .01 and abs(T0 - Tref0) < .1): return -np.inf
    m, aR, inc = model(theta, t)
    if aR < 2: return -np.inf
    s2 = ef ** 2 + np.exp(2 * lns)
    return -.5 * np.sum((f - m) ** 2 / s2 + np.log(s2)) - .5 * ((rho - RHO0) / RHO_E) ** 2

ndim, nw = 6, 40
p0 = np.array([P0, Tref0, np.sqrt(1.5e-3), .9, RHO0, -8.])
pos = p0 + np.array([1e-5, 1e-3, 3e-3, .1, 1., .2]) * np.random.default_rng(0).standard_normal((nw, ndim))
pos[:, 3] = np.abs(pos[:, 3])
sampler = emcee.EnsembleSampler(nw, ndim, lnprob)
sampler.run_mcmc(pos, 12000, progress=False)
try: tau = sampler.get_autocorr_time(tol=0); print('autocorr times:', np.round(tau, 0))
except Exception as ex: tau = None; print('autocorr', ex)
ch = sampler.get_chain(discard=4000, thin=10, flat=True)
P, T0, k, b, rho, lns = ch.T
aR = (rho * 6.674e-8 * (P * 86400) ** 2 / (3 * np.pi)) ** (1 / 3); inc = np.degrees(np.arccos(b / aR))
T14 = P / np.pi * np.arcsin(np.sqrt(np.clip((1 + k) ** 2 - b ** 2, 0, None)) / (aR * np.sin(np.radians(inc)))) * 24
dep = np.array([1 - model([P[i], T0[i], k[i], b[i], rho[i], 0], np.array([T0[i]]))[0][0] for i in np.random.default_rng(1).choice(len(P), 3000)]) * 1e6
Rp = k * np.random.default_rng(2).normal(RS, RS_E, len(k)) * 109.08
Teq = TEFF * np.sqrt(1 / (2 * aR))                                          # A_B = 0, full redistribution
S = (RS ** 2 * (TEFF / 5772) ** 4) / ((aR * RS * 0.00465047) ** 2)
def q(x, d=4): lo, me, hi = np.percentile(x, [16, 50, 84]); return me, hi - me, me - lo
res = {name: q(x) for name, x in [('period_d', P), ('T0_BTJD', T0), ('Rp_Rs', k), ('b', b), ('rho_star_gcc', rho), ('a_Rs', aR), ('inc_deg', inc),
                                  ('T14_h', T14), ('depth_ppm', dep), ('Rp_Rearth', Rp), ('Teq_K', Teq), ('insolation_Searth', S), ('jitter_ppm', np.exp(lns) * 1e6)]}
res['T0_BJD'] = (res['T0_BTJD'][0] + 2457000, res['T0_BTJD'][1], res['T0_BTJD'][2])
for kk, (me, up, lo) in res.items(): print(f'{kk:18s} {me:.6f} +{up:.6f} -{lo:.6f}')
print('acceptance', np.mean(sampler.acceptance_fraction).round(3))
json.dump({k_: list(map(float, v)) for k_, v in res.items()}, open(f'{OUT}/fit_results{TAG}.json', 'w'), indent=1)
np.save(f'{OUT}/chain{TAG}.npy', ch)

# figure: phase-folded data with best model, residuals
best = np.median(ch, 0); ph_h = (((t - best[1]) / best[0] + .5) % 1 - .5) * best[0] * 24
o = np.argsort(ph_h); tt = best[1] + np.linspace(-4, 4, 800) / 24; mm = model(best, tt)[0]
fig, ax = plt.subplots(2, 1, figsize=(8, 7), sharex=True, gridspec_kw=dict(height_ratios=[3, 1]))
ax[0].plot(ph_h, (f - 1) * 1e3, '.', color='.75', ms=2, label='TESS 2-min (S17-S52)')
eb = np.arange(-4, 4.01, 10 / 60); idx = np.digitize(ph_h, eb); bx = (eb[:-1] + eb[1:]) / 2
by = np.array([np.mean(f[idx == i]) if np.sum(idx == i) > 5 else np.nan for i in range(1, len(eb))]); be = np.array([np.std(f[idx == i]) / np.sqrt(np.sum(idx == i)) if np.sum(idx == i) > 5 else np.nan for i in range(1, len(eb))])
ax[0].errorbar(bx, (by - 1) * 1e3, be * 1e3, fmt='ok', ms=4, label='10-min bins'); ax[0].plot((tt - best[1]) * 24, (mm - 1) * 1e3, 'r-', lw=2, label='median-posterior model')
ax[0].set_ylabel('relative flux - 1 (ppt)'); ax[0].set_ylim(-6, 3); ax[0].legend(fontsize=8); ax[0].set_title(f'TIC 366760231: P = {res["period_d"][0]:.5f} d, Rp = {res["Rp_Rearth"][0]:.2f} R$_\\oplus$')
resid = f - model(best, t)[0]; br = np.array([np.mean(resid[idx == i]) if np.sum(idx == i) > 5 else np.nan for i in range(1, len(eb))])
ax[1].plot(ph_h, resid * 1e3, '.', color='.75', ms=2); ax[1].plot(bx, br * 1e3, 'ok', ms=4); ax[1].axhline(0, c='r'); ax[1].set_ylim(-4, 4)
ax[1].set_xlabel('hours from mid-transit'); ax[1].set_ylabel('resid (ppt)'); ax[1].set_xlim(-4, 4)
plt.tight_layout(); plt.savefig(f'{OUT}/fig1_transit_fit{TAG}.png', dpi=130)
