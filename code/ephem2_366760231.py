"""Global (P, T0) likelihood for both TIC 366760231 signals using the best-fit transit-model shape (not a box) on unbinned data.
Also per-sector mid-times for a TTV/consistency check. Writes two_planet/ephemeris_final.json and upcoming tables (UTC and KST)."""
import glob, json, numpy as np, warnings, lightkurve as lk, batman
from astropy.time import Time
warnings.filterwarnings('ignore')
D = '/home/sarang/tess_search/deep/366760231'; O = '/home/sarang/tess_search/cand_366760231/two_planet'
eph = json.load(open(f'{O}/ephemeris2.json')); fit = json.load(open(f'{O}/fit2_results.json')); rho = fit['rho_star_gcc'][0]
G = {n: (eph[n]['P'], eph[n]['T0_BTJD'] + fit[n]['dt_min'][0] / 1440, fit[n]['T14_h'][0] / 24, fit[n]['Rp_Rs'][0], fit[n]['b'][0]) for n in 'bc'}
T, F, E, S, C = [], [], [], [], []
for fn in sorted(glob.glob(f'{D}/*s_lc.fits')) + sorted(glob.glob(f'{D}/hlsp_tess-spoc*_lc.fits')):
    lc = lk.read(fn, quality_bitmask='default', flux_column='pdcsap_flux').remove_nans().normalize(); lc = lc[np.isfinite(lc.flux.value)]
    t = lc.time.value; cad = np.nanmedian(np.diff(t)); m = np.zeros(len(t), bool)
    for P, T0, d, k, b in G.values(): m |= np.abs(((t - T0) / P + .5) % 1 - .5) * P < 1.0 * d
    fl = lc.flatten(window_length=int(round(.75 / cad)) | 1, mask=m).remove_outliers(sigma_upper=4, sigma_lower=50)
    t = np.asarray(fl.time.value, float); f = np.ma.filled(np.ma.asarray(fl.flux.value, dtype=float), np.nan); ok = np.isfinite(f); t, f = t[ok], f[ok]
    T.append(t); F.append(f); E.append(np.full(len(t), 1.4826 * np.median(np.abs(f - np.median(f))))); S.append(np.full(len(t), lc.meta['SECTOR'])); C.append(np.full(len(t), cad))
t, f, e, sec, cad = map(np.concatenate, (T, F, E, S, C))
pm = batman.TransitParams(); pm.ecc = 0.; pm.w = 90.; pm.u = [0.20, 0.40]; pm.limb_dark = 'quadratic'
def shape(name, x, exp):                       # model vs time-from-mid-transit x (days)
    P, T0, d, k, b = G[name]; aR = (rho * 6.674e-8 * (P * 86400) ** 2 / (3 * np.pi)) ** (1 / 3)
    pm.per, pm.t0, pm.rp, pm.a, pm.inc = P, 0., k, aR, np.degrees(np.arccos(b / aR))
    return batman.TransitModel(pm, x, supersample_factor=5, exp_time=exp).light_curve(pm)
out = {}
for name, other in (('b', 'c'), ('c', 'b')):
    P0, T00, d0, k0, b0 = G[name]; Po, To, do = G[other][:3]
    keep = ~(np.abs(((t - To) / Po + .5) % 1 - .5) * Po < 0.75 * do); tt, ff, ee, ss, cc = t[keep], f[keep], e[keep], sec[keep], cad[keep]
    xg = np.linspace(-0.2, 0.2, 4001); tab = {c_: shape(name, xg, c_) for c_ in np.unique(np.round(cc, 5))}     # tabulated shape per cadence
    cr = np.round(cc, 5); n0 = np.round((np.average(tt, weights=1 / ee ** 2) - T00) / P0); Tc0 = T00 + n0 * P0; B = np.ptp(tt)
    def chi2(P, Tc, sel=slice(None)):
        ph = (((tt[sel] - Tc) / P + .5) % 1 - .5) * P; near = np.abs(ph) < 0.2; m = np.ones(len(ph))
        for c_ in tab: s_ = near & (cr[sel] == c_); m[s_] = np.interp(ph[s_], xg, tab[c_])
        return np.sum(((ff[sel] - m) / ee[sel]) ** 2)
    Ps = P0 + np.linspace(-1, 1, 201) * 2.5 * d0 * P0 / B; dTs = np.linspace(-1, 1, 161) * 0.5 * d0
    chi = np.array([[chi2(P, Tc0 + dT) for dT in dTs] for P in Ps]); i, j = np.unravel_index(np.argmin(chi), chi.shape)
    L = np.exp(-.5 * (chi - chi.min())); pP = L.sum(1) / L.sum(); pT = L.sum(0) / L.sum()
    mP = np.sum(Ps * pP); sP = np.sqrt(np.sum((Ps - mP) ** 2 * pP)); mT = np.sum(dTs * pT); sT = np.sqrt(np.sum((dTs - mT) ** 2 * pT))
    cov = np.sum(np.outer(Ps - mP, dTs - mT) * L) / L.sum(); null = np.sum(((ff - 1) / ee) ** 2)
    out[name] = dict(P=float(mP), sP=float(sP), T0_BTJD=float(Tc0 + mT), T0_BJD=float(Tc0 + mT + 2457000), sT0=float(sT), cov_PT=float(cov), dchi2_vs_flat=float(null - chi.min()), epoch_ref=int(n0))
    print(f'{name}: P = {mP:.7f} ± {sP:.7f} d; T0 = {Tc0 + mT + 2457000:.5f} ± {sT:.5f} BJD_TDB; corr {cov/(sP*sT):+.2f}; delta-chi2 vs no transit {null - chi.min():.1f}; edge check P[{pP[0]:.1e},{pP[-1]:.1e}] T[{pT[0]:.1e},{pT[-1]:.1e}]')
    for s_ in np.unique(ss):                                    # per-sector mid-time offsets at the global period
        sel = ss == s_; c1 = np.array([chi2(mP, Tc0 + mT + dT, sel) for dT in dTs]); p1 = np.exp(-.5 * (c1 - c1.min())); p1 /= p1.sum(); m1 = np.sum(dTs * p1)
        print(f'    S{s_}: O-C {m1*1440:+6.1f} ± {np.sqrt(np.sum((dTs - m1) ** 2 * p1))*1440:4.1f} min')
json.dump(out, open(f'{O}/ephemeris_final.json', 'w'), indent=1)
old = json.load(open('/home/sarang/tess_search/cand_366760231/ephemeris.json'))
now = Time('2026-11-16').jd - 2457000
for name in 'bc':
    P, T0, sP, sT, cv = (out[name][k] for k in ('P', 'T0_BTJD', 'sP', 'sT0', 'cov_PT')); n1 = int(np.ceil((now - T0) / P)); L_ = ['| epoch | BJD_TDB | UTC | KST | 1σ |', '|---|---|---|---|---|']
    for n in range(n1, n1 + (24 if name == 'b' else 60)):
        tc = T0 + n * P + 2457000; s = np.sqrt(sT ** 2 + (n * sP) ** 2 + 2 * n * cv) * 1440
        L_.append(f'| {n} | {tc:.4f} | {Time(tc, format="jd").iso[:16]} | {Time(tc + 9 / 24, format="jd").iso[:16]} | ±{s:.0f} min |')
    open(f'{O}/upcoming_{name}.md', 'w').write('\n'.join(L_) + '\n'); print(name, 'next:', L_[2], L_[3])
    if name == 'b':
        n_old = np.round((T0 + n1 * P + 2457000 - old['T0_BJD']) / old['P']); told = old['T0_BJD'] + n_old * old['P']
        print(f'   published-ephemeris prediction for the same transit: {told:.4f} -> new minus old = {(T0 + n1 * P + 2457000 - told)*1440:+.1f} min')
