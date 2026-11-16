"""MCMC transit fit of a deep-pass signal on the raw PDCSAP light curves of every sector (re-detrended with the transits masked).
Free: P, T0 (epoch near the data centroid), k = Rp/R*, b, rho* (prior from TIC mass and radius; FREE_RHO=1 removes it), ln(jitter).
Each cadence (1800/600/200/120 s) is modelled with its own exposure-time integration.  Usage: deep2_fit.py TIC K"""
import sys, os, json, glob, numpy as np, pandas as pd, urllib.request, warnings, lightkurve as lk, batman, emcee
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')
R = '/home/sarang/tess_search'; TIC, K = int(sys.argv[1]), int(sys.argv[2]); TAG = '_freerho' if os.environ.get('FREE_RHO') else ''
v = json.load(open(f'{R}/deep2/vet/{TIC}_{K}.json')); P0, T00, D0 = v['P'], v['T0'], max(v['T14_h'] / 24, 0.03)
st = pd.concat([pd.read_csv(f'{R}/data/{n}.csv') for n in ('targets', 'targets_allsky', 'targets_ffi_faint')]).drop_duplicates('ID').set_index('ID').loc[TIC]
MS, RS, TEFF = float(st['mass']), float(st.rad), float(st.Teff); MS_E, RS_E = 0.20 * MS, 0.05 * RS            # generous TIC-like errors for M dwarfs
U = (0.25, 0.35) if TEFF < 3600 else (0.35, 0.28); RHO0 = 1.41 * MS / RS ** 3; RHO_E = 1e3 if TAG else RHO0 * np.hypot(MS_E / MS, 3 * RS_E / RS)
D = f'{R}/deep2/pix/{TIC}'; os.makedirs(D, exist_ok=True); t16 = f'{TIC:016d}'
plan = json.loads(pd.read_csv(f'{R}/deep2/targets.csv', usecols=['ID', 'plan']).set_index('ID').loc[TIC, 'plan'])
sigs = pd.read_csv(f'{R}/deep2_signals.csv', low_memory=False); oth = sigs[(sigs.tic == TIC) & (sigs.k != K) & (sigs.snr >= 7.5) & (sigs.ntr >= 3)]
oth = oth[[min(abs(p / P0 - q) for q in (1, .5, 2)) > 0.01 for p in oth.P]]
T, F, E, S, X = [], [], [], [], []
for s, src in sorted((int(k), x) for k, x in plan.items()):
    if src == 'ffi': fn = f'hlsp_tess-spoc_tess_phot_{t16}-s{s:04d}_tess_v1_lc.fits'; url = f'https://archive.stsci.edu/hlsps/tess-spoc/s{s:04d}/target/{t16[0:4]}/{t16[4:8]}/{t16[8:12]}/{t16[12:16]}/{fn}'
    else: fn = src[3:]; x = fn.split('-')[2]; url = f'https://stpubdata.s3.amazonaws.com/tess/public/tid/s{s:04d}/{x[0:4]}/{x[4:8]}/{x[8:12]}/{x[12:16]}/{fn}'
    p = f'{D}/{fn}'
    for _ in range(4):
        if os.path.exists(p) and os.path.getsize(p) > 20000: break
        try: urllib.request.urlretrieve(url, p)
        except Exception as ex: print('  lc', s, ex)
    try: lc = lk.read(p, quality_bitmask='default', flux_column='pdcsap_flux').remove_nans().normalize(); lc = lc[np.isfinite(lc.flux.value)]
    except Exception as ex: print('  skip', s, ex); continue
    tt = lc.time.value; cad = np.nanmedian(np.diff(tt)); mask = np.abs(((tt - T00) / P0 + .5) % 1 - .5) * P0 < 1.5 * D0
    for _, o in oth.iterrows(): mask |= np.abs(((tt - o.t0) / o.P + .5) % 1 - .5) * o.P < 1.0 * o.dur
    fl = lc.flatten(window_length=int(round(.75 / cad)) | 1, mask=mask).remove_outliers(sigma_upper=4, sigma_lower=50)
    t = np.asarray(fl.time.value, float); f = np.ma.filled(np.ma.asarray(fl.flux.value, dtype=float), np.nan); ef = np.ma.filled(np.ma.asarray(fl.flux_err.value, dtype=float), np.nan)
    ok = (np.abs(((t - T00) / P0 + .5) % 1 - .5) * P0 < 4 * D0) & np.isfinite(f) & np.isfinite(ef)
    for _, o in oth.iterrows(): ok &= ~(np.abs(((t - o.t0) / o.P + .5) % 1 - .5) * o.P < 1.0 * o.dur)
    # trim orbit edges (ramps)
    g = np.flatnonzero(np.diff(tt) > 0.5)
    for ts in np.r_[tt[0], tt[g + 1]]: ok &= ~((t >= ts) & (t < ts + 0.2))
    oot = np.abs(((t - T00) / P0 + .5) % 1 - .5) * P0 > 1.2 * D0; sc = 1.4826 * np.median(np.abs(f[ok & oot] - np.median(f[ok & oot]))) if (ok & oot).sum() > 20 else np.nanmedian(ef)
    T.append(t[ok]); F.append(f[ok]); E.append(np.full(ok.sum(), sc)); S.append(np.full(ok.sum(), s)); X.append(np.full(ok.sum(), round(cad * 86400 / 10) * 10))
t, f, ef, sec, exp = map(np.concatenate, (T, F, E, S, X)); nref = np.round(np.median(np.round((t - T00) / P0))); Tref0 = T00 + nref * P0
print(f'TIC {TIC} k={K}: {len(t)} points near transit in sectors {sorted(set(sec.astype(int)))}; exposure groups {sorted(set(exp.astype(int)))} s; T0 guess {Tref0:.4f}; rho* prior {RHO0:.1f} ± {RHO_E:.1f} g/cc', flush=True)
groups = [(x, exp == x) for x in np.unique(exp)]
params = batman.TransitParams(); params.t0 = 0.; params.per = P0; params.rp = .05; params.a = 30.; params.inc = 89.; params.ecc = 0.; params.w = 90.; params.u = list(U); params.limb_dark = 'quadratic'
def aR_of(P, rho): return (rho * 6.674e-8 * (P * 86400) ** 2 / (3 * np.pi)) ** (1 / 3)
def model(theta, tt, x):
    P, T0, k, b, rho, lns = theta; aR = aR_of(P, rho); params.per, params.t0, params.rp, params.a, params.inc = P, T0, k, aR, np.degrees(np.arccos(min(b / aR, 1)))
    return batman.TransitModel(params, tt, supersample_factor=max(1, int(round(x / 150))) | 1, exp_time=x / 86400).light_curve(params)
def lnprob(theta):
    P, T0, k, b, rho, lns = theta
    if not (0 < k < .6 and 0 <= b < 1 + k and .3 < rho < 300 and -12 < lns < -2 and abs(P - P0) < .02 and abs(T0 - Tref0) < .1) or aR_of(P, rho) < 1.5: return -np.inf
    m = np.empty(len(t))
    for x, g in groups: m[g] = model(theta, t[g], x)
    s2 = ef ** 2 + np.exp(2 * lns); return -.5 * np.sum((f - m) ** 2 / s2 + np.log(s2)) - .5 * ((rho - RHO0) / RHO_E) ** 2
ndim, nw = 6, 40; p0 = np.array([P0, Tref0, np.sqrt(max(v['trap_depth_ppm'], 50) * 1e-6), .3, RHO0, -7.])
pos = p0 + np.array([P0 ** 2 / np.ptp(t) * 0.02, 1e-3, 0.05 * p0[2], .1, 0.05 * RHO0, .2]) * np.random.default_rng(0).standard_normal((nw, ndim)); pos[:, 3] = np.abs(pos[:, 3])
sm = emcee.EnsembleSampler(nw, ndim, lnprob); sm.run_mcmc(pos, int(os.environ.get('NSTEP', 9000)), progress=False)
try: print('autocorr times:', np.round(sm.get_autocorr_time(tol=0), 0))
except Exception as ex: print('autocorr', ex)
ch = sm.get_chain(discard=3000, thin=10, flat=True); P, T0, k, b, rho, lns = ch.T; aR = aR_of(P, rho); inc = np.degrees(np.arccos(np.clip(b / aR, 0, 1)))
T14 = P / np.pi * np.arcsin(np.clip(np.sqrt(np.clip((1 + k) ** 2 - b ** 2, 0, None)) / (aR * np.sin(np.radians(inc))), 0, 1)) * 24
ii = np.random.default_rng(1).choice(len(P), 2000); dep = np.array([1 - model([P[i], T0[i], k[i], b[i], rho[i], 0], np.array([T0[i]]), 20.)[0] for i in ii]) * 1e6
Rp = k * np.random.default_rng(2).normal(RS, RS_E, len(k)) * 109.08; Teq = TEFF * np.sqrt(1 / (2 * aR)); Sin = (RS ** 2 * (TEFF / 5772) ** 4) / ((aR * RS * 0.00465047) ** 2)
def q(x): lo, me, hi = np.percentile(x, [16, 50, 84]); return float(me), float(hi - me), float(me - lo)
res = {n: q(x) for n, x in [('period_d', P), ('T0_BTJD', T0), ('Rp_Rs', k), ('b', b), ('rho_star_gcc', rho), ('a_Rs', aR), ('inc_deg', inc), ('T14_h', T14), ('depth_ppm', dep), ('Rp_Rearth', Rp), ('Teq_K', Teq), ('insolation_Searth', Sin), ('jitter_ppm', np.exp(lns) * 1e6)]}
res['T0_BJD'] = (res['T0_BTJD'][0] + 2457000, res['T0_BTJD'][1], res['T0_BTJD'][2]); res['star'] = dict(Rs=RS, Ms=MS, Teff=TEFF, rho_prior=[RHO0, RHO_E]); res['npts'] = int(len(t)); res['sectors'] = sorted(set(int(x) for x in sec))
for kk, val in res.items():
    if isinstance(val, tuple): print(f'{kk:18s} {val[0]:.6f} +{val[1]:.6f} -{val[2]:.6f}')
print('acceptance', np.mean(sm.acceptance_fraction).round(3)); json.dump(res, open(f'{R}/deep2/vet/{TIC}_{K}.fit{TAG}.json', 'w'), indent=1); np.save(f'{R}/deep2/vet/{TIC}_{K}.chain{TAG}.npy', ch.astype(np.float64))
best = np.median(ch, 0); ph = (((t - best[1]) / best[0] + .5) % 1 - .5) * best[0] * 24; lim = 4 * D0 * 24; tt = np.linspace(-lim, lim, 600) / 24
fig, ax = plt.subplots(2, 1, figsize=(8, 7), sharex=True, gridspec_kw=dict(height_ratios=[3, 1])); ax[0].plot(ph, (f - 1) * 1e3, '.', color='.8', ms=2)
eb = np.linspace(-lim, lim, 41); idx = np.digitize(ph, eb); bx = (eb[:-1] + eb[1:]) / 2; mfull = np.empty(len(t))
for x, g in groups: mfull[g] = model(best, t[g], x)
w = 1 / ef ** 2; by = np.array([np.average(f[idx == i], weights=w[idx == i]) if (idx == i).sum() > 3 else np.nan for i in range(1, len(eb))]); be = np.array([1 / np.sqrt(w[idx == i].sum()) if (idx == i).sum() > 3 else np.nan for i in range(1, len(eb))])
br = np.array([np.average((f - mfull)[idx == i], weights=w[idx == i]) if (idx == i).sum() > 3 else np.nan for i in range(1, len(eb))])
ax[0].errorbar(bx, (by - 1) * 1e3, be * 1e3, fmt='ok', ms=4, label='binned, all sectors'); ax[0].plot(tt * 24, (model(best, best[1] + tt, 20.) - 1) * 1e3, 'r-', lw=2, label='median model (unsmeared)')
ax[0].set_ylim((np.nanmin(by) - 1) * 1e3 * 1.8, -(np.nanmin(by) - 1) * 1e3 * 0.9); ax[0].set_ylabel('relative flux - 1 (ppt)'); ax[0].legend(fontsize=8)
ax[0].set_title(f'TIC {TIC}: P = {res["period_d"][0]:.5f} d, Rp = {res["Rp_Rearth"][0]:.2f} R$_\\oplus$, depth {res["depth_ppm"][0]:.0f} ppm')
ax[1].errorbar(bx, br * 1e3, be * 1e3, fmt='ok', ms=4); ax[1].axhline(0, c='r'); ax[1].set_xlabel('hours from mid-transit'); ax[1].set_ylabel('resid (ppt)'); plt.tight_layout(); plt.savefig(f'{R}/deep2/vet/{TIC}_{K}.fit{TAG}.png', dpi=110)
