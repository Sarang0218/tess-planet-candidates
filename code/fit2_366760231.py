"""Two-planet analysis of TIC 366760231: the published 3.137-d candidate (.01) and the new 1.143-d signal (.02).
1) (P, T0) likelihood grid for each signal over all 9 sectors -> linear ephemerides.
2) Joint MCMC on phase-folded, 3-min-binned transits with ONE shared stellar density (TIC prior, or free with FREE_RHO=1):
   tests whether both transit shapes are consistent with the same star, and gives k, b, T14, depth for each.
Light curves: SPOC 2-min S17,18,19,24,25,52 and TESS-SPOC FFI (200 s) S58,59,85; 0.75-d detrend with both transits masked."""
import glob, os, json, numpy as np, warnings, lightkurve as lk, batman, emcee, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from astropy.time import Time
warnings.filterwarnings('ignore')
D = '/home/sarang/tess_search/deep/366760231'; O = '/home/sarang/tess_search/cand_366760231/two_planet'
MS, RS, MS_E, RS_E, TEFF = 0.3787, 0.3887, 0.0203, 0.0116, 3456.; U = (0.20, 0.40)
RHO0 = 1.41 * MS / RS ** 3; RHO_E = RHO0 * np.hypot(MS_E / MS, 3 * RS_E / RS)
TAG = '_freerho' if os.environ.get('FREE_RHO') else ''
if TAG: RHO_E = 1e4
G0 = {'b': (3.1373207, 1962.04155, 0.6 / 24), 'c': (1.143413, 1999.7835, 1.08 / 24)}
T, F, E = [], [], []
for fn in sorted(glob.glob(f'{D}/*s_lc.fits')) + sorted(glob.glob(f'{D}/hlsp_tess-spoc*_lc.fits')):
    lc = lk.read(fn, quality_bitmask='default', flux_column='pdcsap_flux').remove_nans().normalize(); lc = lc[np.isfinite(lc.flux.value)]
    t = lc.time.value; cad = np.nanmedian(np.diff(t)); m = np.zeros(len(t), bool)
    for P, T0, d in G0.values(): m |= np.abs(((t - T0) / P + .5) % 1 - .5) * P < 1.0 * d
    fl = lc.flatten(window_length=int(round(.75 / cad)) | 1, mask=m).remove_outliers(sigma_upper=4, sigma_lower=50)
    t = np.asarray(fl.time.value, float); f = np.ma.filled(np.ma.asarray(fl.flux.value, dtype=float), np.nan); ok = np.isfinite(f); t, f = t[ok], f[ok]
    T.append(t); F.append(f); E.append(np.full(len(t), 1.4826 * np.median(np.abs(f - np.median(f)))))
t, f, e = map(np.concatenate, (T, F, E)); w = 1 / e ** 2
print(f'{len(t)} points, {len(T)} sectors, BTJD {t.min():.1f}-{t.max():.1f}; rho* prior {RHO0:.2f} ± {RHO_E:.2f} g/cc')

def other_mask(name):                     # drop points inside the other signal's transit
    m = np.ones(len(t), bool)
    for k, (P, T0, d) in G0.items():
        if k != name: m &= ~(np.abs(((t - T0) / P + .5) % 1 - .5) * P < 0.75 * d)
    return m
eph = {}
for name, (P0, T00, d0) in G0.items():
    m = other_mask(name); tt, ff, ww = t[m], f[m], w[m]
    n0 = np.round((np.average(tt, weights=ww) - T00) / P0); Tc0 = T00 + n0 * P0; B = np.ptp(tt)
    Ps = P0 + np.linspace(-1, 1, 161) * 2.0 * d0 * P0 / B; dTs = np.linspace(-1, 1, 121) * 0.8 * d0
    chi = np.zeros((len(Ps), len(dTs)))
    for i, P in enumerate(Ps):
        ph = (((tt - Tc0) / P + .5) % 1 - .5) * P
        for j, dT in enumerate(dTs):
            a = np.abs(ph - dT) < d0 / 2; b_ = (np.abs(ph - dT) > 0.75 * d0) & (np.abs(ph - dT) < 3 * d0)
            dep = np.average(ff[b_], weights=ww[b_]) - np.average(ff[a], weights=ww[a]); chi[i, j] = -dep ** 2 / (1 / ww[a].sum() + 1 / ww[b_].sum())
    i, j = np.unravel_index(np.argmin(chi), chi.shape); ok1 = chi < chi.min() + 1
    Pb, sP = Ps[i], np.ptp(Ps[ok1.any(1)]) / 2; Tb, sT = Tc0 + dTs[j], np.ptp(dTs[ok1.any(0)]) / 2
    # correlation is small because the pivot epoch is the weighted centre of the data
    eph[name] = dict(P=float(Pb), sP=float(max(sP, np.diff(Ps)[0])), T0_BTJD=float(Tb), T0_BJD=float(Tb + 2457000), sT0=float(max(sT, np.diff(dTs)[0])), box_snr=float(np.sqrt(-chi.min())))
    print(f'{name}: P = {Pb:.6f} ± {eph[name]["sP"]:.6f} d, T0 = {Tb + 2457000:.5f} ± {eph[name]["sT0"]:.5f} BJD_TDB, box S/N {np.sqrt(-chi.min()):.1f}')
json.dump(eph, open(f'{O}/ephemeris2.json', 'w'), indent=1)

# phase-folded, binned transits
fold = {}
for name in G0:
    P, T0 = eph[name]['P'], eph[name]['T0_BTJD']; m = other_mask(name); ph = (((t[m] - T0) / P + .5) % 1 - .5) * P
    edges = np.arange(-3.0, 3.0001, 3 / 60) / 24; idx = np.digitize(ph, edges); x, y, ye = [], [], []
    for q in range(1, len(edges)):
        s = idx == q
        if s.sum() > 5: x.append(ph[s].mean()); y.append(np.average(f[m][s], weights=w[m][s])); ye.append(1 / np.sqrt(w[m][s].sum()))
    fold[name] = (np.array(x), np.array(y), np.array(ye), P)
pm = batman.TransitParams(); pm.t0 = 0.; pm.per = 1.; pm.rp = .05; pm.a = 10.; pm.inc = 89.; pm.ecc = 0.; pm.w = 90.; pm.u = list(U); pm.limb_dark = 'quadratic'
def model(P, dt, k, b, rho, x):
    aR = (rho * 6.674e-8 * (P * 86400) ** 2 / (3 * np.pi)) ** (1 / 3)
    if aR <= b + 0.01: return None, aR
    pm.per, pm.t0, pm.rp, pm.a, pm.inc = P, dt, k, aR, np.degrees(np.arccos(b / aR))
    return batman.TransitModel(pm, x, supersample_factor=5, exp_time=3 / 1440).light_curve(pm), aR
def lnprob(th):
    rho, kb, bb, dtb, kc, bc, dtc, lns = th
    if not (0.5 < rho < 400 and 0 < kb < .3 and 0 <= bb < 1 + kb and 0 < kc < .3 and 0 <= bc < 1 + kc and abs(dtb) < .02 and abs(dtc) < .03 and -12 < lns < -5): return -np.inf
    lp = -.5 * ((rho - RHO0) / RHO_E) ** 2
    for name, k, b, dt in (('b', kb, bb, dtb), ('c', kc, bc, dtc)):
        x, y, ye, P = fold[name]; mm, aR = model(P, dt, k, b, rho, x)
        if mm is None: return -np.inf
        s2 = ye ** 2 + np.exp(2 * lns); lp += -.5 * np.sum((y - mm) ** 2 / s2 + np.log(s2))
    return lp
p0 = np.array([RHO0 if not TAG else 12., .045, .9, 0., .023, .4, 0., -9.]); nw = 48
pos = p0 + np.array([.5, 3e-3, .05, 1e-3, 2e-3, .15, 1e-3, .2]) * np.random.default_rng(0).standard_normal((nw, 8)); pos[:, 2] = np.abs(pos[:, 2]); pos[:, 5] = np.abs(pos[:, 5])
sm = emcee.EnsembleSampler(nw, 8, lnprob); sm.run_mcmc(pos, 14000, progress=False)
ch = sm.get_chain(discard=5000, thin=10, flat=True); np.save(f'{O}/chain2{TAG}.npy', ch)
try: print('autocorr', np.round(sm.get_autocorr_time(tol=0), 0))
except Exception as ex: print(ex)
def q(x): lo, me, hi = np.percentile(x, [16, 50, 84]); return [float(me), float(hi - me), float(me - lo)]
res = dict(rho_star_gcc=q(ch[:, 0]), rho_prior=[RHO0, RHO_E]); rs = np.random.default_rng(2).normal(RS, RS_E, len(ch))
for name, ik, ib, idt in (('b', 1, 2, 3), ('c', 4, 5, 6)):
    P = eph[name]['P']; k, b, rho = ch[:, ik], ch[:, ib], ch[:, 0]; aR = (rho * 6.674e-8 * (P * 86400) ** 2 / (3 * np.pi)) ** (1 / 3); inc = np.degrees(np.arccos(b / aR))
    T14 = P / np.pi * np.arcsin(np.clip(np.sqrt(np.clip((1 + k) ** 2 - b ** 2, 0, None)) / (aR * np.sin(np.radians(inc))), 0, 1)) * 24
    sel = np.random.default_rng(1).choice(len(ch), 2000); dep = np.array([1 - model(P, 0., k[i], b[i], rho[i], np.array([0.]))[0][0] for i in sel]) * 1e6
    res[name] = dict(P=P, Rp_Rs=q(k), b=q(b), a_Rs=q(aR), inc_deg=q(inc), T14_h=q(T14), depth_ppm=q(dep), Rp_Rearth=q(k * rs * 109.08), Teq_K=q(TEFF * np.sqrt(1 / (2 * aR))),
                     insolation=q(RS ** 2 * (TEFF / 5772) ** 4 / (aR * RS * 0.00465047) ** 2), dt_min=q(ch[:, idt] * 1440))
res['jitter_ppm'] = q(np.exp(ch[:, 7]) * 1e6); res['acceptance'] = float(np.mean(sm.acceptance_fraction))
json.dump(res, open(f'{O}/fit2_results{TAG}.json', 'w'), indent=1)
for name in 'bc': print(name, {kk: (np.round(v, 4).tolist() if isinstance(v, list) else round(v, 6)) for kk, v in res[name].items()})
print('rho*', np.round(res['rho_star_gcc'], 2), 'jitter', np.round(res['jitter_ppm'], 0), 'acc', round(res['acceptance'], 3))
best = np.median(ch, 0); fig, ax = plt.subplots(1, 2, figsize=(12, 4.5), sharey=True)
for a, name, ik, ib, idt in ((ax[0], 'b', 1, 2, 3), (ax[1], 'c', 4, 5, 6)):
    x, y, ye, P = fold[name]; xb = np.arange(-3, 3.001, .25) / 24; idx = np.digitize(x, xb)
    a.errorbar(x * 24, (y - 1) * 1e6, ye * 1e6, fmt='.', color='.8', ms=2, lw=.5, zorder=1)
    a.errorbar([x[idx == q_].mean() * 24 for q_ in range(1, len(xb)) if (idx == q_).any()], [(np.average(y[idx == q_], weights=1 / ye[idx == q_] ** 2) - 1) * 1e6 for q_ in range(1, len(xb)) if (idx == q_).any()],
               [1e6 / np.sqrt(np.sum(1 / ye[idx == q_] ** 2)) for q_ in range(1, len(xb)) if (idx == q_).any()], fmt='ok', ms=4, zorder=3)
    xx = np.linspace(-3, 3, 600) / 24; a.plot(xx * 24, (model(P, best[idt], best[ik], best[ib], best[0], xx)[0] - 1) * 1e6, 'r-', lw=2, zorder=4)
    a.set_title(f'TIC 366760231 {"candidate .01" if name == "b" else "new signal .02"}: P = {P:.5f} d, Rp = {res[name]["Rp_Rearth"][0]:.2f} R$_\\oplus$, T14 = {res[name]["T14_h"][0]:.2f} h'); a.set_xlabel('hours from mid-transit'); a.set_ylim(-3000, 2000)
ax[0].set_ylabel('relative flux − 1 (ppm)'); plt.tight_layout(); plt.savefig(f'{O}/fig_two_planets{TAG}.png', dpi=110)
now = Time('2026-11-16').jd - 2457000
for name in 'bc':
    P, T0, sP, sT = eph[name]['P'], eph[name]['T0_BTJD'], eph[name]['sP'], eph[name]['sT0']; n1 = int(np.ceil((now - T0) / P)); L = ['| epoch | BJD_TDB | UTC | KST | 1σ |', '|---|---|---|---|---|']
    for n in range(n1, n1 + (20 if name == 'b' else 45)):
        tc = T0 + n * P + 2457000; s = np.hypot(sT, n * sP) * 1440; L.append(f'| {n} | {tc:.4f} | {Time(tc, format="jd").iso[:16]} | {Time(tc + 9 / 24, format="jd").iso[:16]} | ±{s:.0f} min |')
    open(f'{O}/upcoming_{name}{TAG}.md', 'w').write('\n'.join(L) + '\n')
