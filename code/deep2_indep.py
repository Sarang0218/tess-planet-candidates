"""Independent re-reduction check of a candidate.  Starts from the RAW aperture flux (SAP_FLUX, no pipeline systematics correction) in the
light-curve files in deep2/pix/{TIC}/, detrends each sector with a Savitzky-Golay filter (lightkurve.flatten, no knowledge of the transit),
then with the fitted transit shape held fixed measures: (1) the depth relative to the published fit and its significance against 300 other
orbital phases; (2) where the best period lies in a +-0.3 % window; (3) the mid-transit time in each observing season against the linear
ephemeris; (4) the depth in the first and second half of the data and by cadence type.   Usage: deep2_indep.py TIC K"""
import sys, glob, json, numpy as np, warnings, lightkurve as lk, batman
from astropy.io import fits
warnings.filterwarnings('ignore')
R = '/home/sarang/tess_search'; TIC, K = int(sys.argv[1]), int(sys.argv[2]); V = f'{R}/deep2/vet'
f = json.load(open(f'{V}/{TIC}_{K}.fit.json')); P, T0 = f['period_d'][0], f['T0_BJD'][0] - 2457000; T14 = f['T14_h'][0] / 24
T, F, S, C = [], [], [], []; crowd = []
import os; QLP = os.environ.get('SRC') == 'qlp'              # SRC=qlp: MIT Quick-Look Pipeline light curves (independent photometry and detrending) from deep2/qlp/{TIC}/
for p in sorted(glob.glob(f'{R}/deep2/qlp/{TIC}/*_llc.fits') if QLP else glob.glob(f'{R}/deep2/pix/{TIC}/*_lc.fits')):
    try:
        h = fits.open(p); d = h[1].data; sec = h[0].header.get('SECTOR'); crowd.append(np.nan if QLP else h[1].header.get('CROWDSAP', np.nan)); col = ('SAP_FLUX' if not QLP else 'KSPSAP_FLUX' if 'KSPSAP_FLUX' in d.names else 'DET_FLUX')
        ok = (d['QUALITY'] == 0) & np.isfinite(d[col]) & np.isfinite(d['TIME']); t, y = np.asarray(d['TIME'][ok], float), np.asarray(d[col][ok], float)
        if len(t) < 200: continue
        gap = np.r_[True, np.diff(t) > 0.5]; start = np.maximum.accumulate(np.where(gap, t, -np.inf)); keep = (t - start > 0.3)               # drop the scattered-light ramp after each gap
        t, y = t[keep], y[keep]; cad = np.nanmedian(np.diff(t)); wl = int(round(1.0 / cad)) | 1
        fl = lk.LightCurve(time=t, flux=y / np.nanmedian(y)).flatten(window_length=wl, polyorder=2, break_tolerance=int(0.5 / cad), niters=3, sigma=3)
        yy = np.asarray(fl.flux.value, float); m = np.isfinite(yy) & (yy < 1 + 5 * 1.4826 * np.nanmedian(np.abs(yy - np.nanmedian(yy))))
        T.append(t[m]); F.append(yy[m]); S.append(np.full(m.sum(), sec)); C.append(np.full(m.sum(), cad))
    except Exception as ex: print('  skip', p.split('/')[-1], ex)
t, y, s, c = [np.concatenate(x) for x in (T, F, S, C)]; sig = np.array([1.4826 * np.nanmedian(np.abs(y[s == u] - np.nanmedian(y[s == u]))) for u in s[np.r_[True, np.diff(s) != 0]]]); e = np.concatenate([np.full((s == u).sum(), v) for u, v in zip(s[np.r_[True, np.diff(s) != 0]], sig)])
pr = batman.TransitParams(); pr.t0 = 0; pr.per = P; pr.rp = f['Rp_Rs'][0]; pr.a = f['a_Rs'][0]; pr.inc = f['inc_deg'][0]; pr.ecc = 0; pr.w = 90; pr.u = [0.25, 0.35]; pr.limb_dark = 'quadratic'
def shape(tt, cc, per=P, t0=T0):
    ph = ((tt - t0) / per + .5) % 1 - .5; out = np.zeros(len(tt)); pr.per = per
    for cv in np.unique(np.round(cc, 4)):
        mm = np.round(cc, 4) == cv; out[mm] = 1 - batman.TransitModel(pr, ph[mm] * per, supersample_factor=7 if cv > 0.004 else 1, exp_time=cv).light_curve(pr)
    return out
def amp(tt, yy, ee, cc, per=P, t0=T0):                     # least-squares scale of the published transit shape (1 = same depth), local baseline removed
    ph = np.abs(((tt - t0) / per + .5) % 1 - .5) * per; near = ph < 3 * T14
    if near.sum() < 10: return np.nan, np.nan, np.inf
    q = shape(tt[near], cc[near], per, t0); w = 1 / ee[near] ** 2; d = 1 - yy[near]; d0 = d - np.sum(w * d) / w.sum(); q0 = q - np.sum(w * q) / w.sum()
    a = np.sum(w * q0 * d0) / np.sum(w * q0 * q0); return a, 1 / np.sqrt(np.sum(w * q0 * q0)), np.sum(w * (d0 - a * q0) ** 2)
a0, ea, _ = amp(t, y, e, c); rng = np.random.default_rng(0); lo = 4 * T14 / P
nulls = np.array([amp(t, y, e, c, t0=T0 + sh * P)[0] for sh in rng.uniform(lo, 1 - lo, 300)]); nulls = nulls[np.isfinite(nulls)]
print(f'TIC {TIC} k={K}: ' + ('QLP light curves' if QLP else 'raw aperture flux') + f', {len(np.unique(s))} sectors, {len(t)} points; pipeline crowding metric (fraction of aperture flux from the target) {np.nanmean(crowd):.3f}')
print(f'  depth relative to the published fit: {a0:.2f} ± {ea:.2f} (white noise);  same measurement at 300 other phases: mean {nulls.mean():+.2f}, scatter {nulls.std():.2f}  ->  {(a0 - nulls.mean()) / nulls.std():.1f} sigma;  expected {np.nanmean(crowd):.2f} if the dip is on the target (raw flux is diluted by neighbours)')
pg = P * (1 + np.linspace(-0.003, 0.003, 241)); ch = np.array([amp(t, y, e, c, per=x, t0=T0 + round((np.median(t) - T0) / P) * (x - P) * 0)[2] for x in pg])      # T0 fixed at the published epoch: a wrong period drifts out of phase
az = np.array([amp(t, y, e, c, per=x)[0] / ea for x in pg]); best = pg[np.argmax(az)]
print(f'  period scan ±0.3 %: strongest at P = {best:.6f} d (published {P:.6f}, difference {abs(best - P) / P * 1e6:.0f} ppm of P);  SNR there {az.max():.1f}, median elsewhere in the window {np.median(az):.1f}')
seas = np.cumsum(np.r_[0, np.diff(t) > 100]); rows = []
for g in np.unique(seas):
    m = seas == g; n = round((np.median(t[m]) - T0) / P); dts = np.linspace(-0.06, 0.06, 121) * min(1, P / 0.5)
    rr = [amp(t[m], y[m], e[m], c[m], t0=T0 + dt) for dt in dts]; aa = np.array([r[0] for r in rr]); x2 = -np.where(aa > 0, (aa / np.array([r[1] for r in rr])) ** 2, 0.)      # detection chi2 gain of the template at each shift (the fit window moves with the shift, so raw chi2 values are not comparable)
    if not np.isfinite(x2).all() or aa[np.argmin(x2)] <= 0: rows.append((g, m.sum(), np.nan, np.nan, np.nan)); continue
    sc = 1.0; inr = dts[x2 - x2.min() < sc]; zz = amp(t[m], y[m], e[m], c[m]); zz = zz[0] / zz[1]
    if zz < 4: rows.append((g, m.sum(), np.nan, np.nan, zz)); continue                 # too weak in this season to time
    rows.append((g, m.sum(), dts[np.argmin(x2)] * 1440, max((inr.max() - inr.min()) / 2, dts[1] - dts[0], 0.4 * T14 / zz) * 1440, zz))
print('  mid-transit time per observing season, observed minus the linear ephemeris (minutes):')
for g, n, oc, eo, z in rows: print(f'    season {int(g)+1} (sectors {int(s[seas==g].min())}-{int(s[seas==g].max())}): ' + (f'{oc:+.1f} ± {eo:.1f}   (transit SNR in this season {z:.1f})' if np.isfinite(oc) else f'too weak to time (SNR {z:.1f})' if np.isfinite(z) else 'not measurable'))
ok = [r for r in rows if np.isfinite(r[2])]
if len(ok) > 1: x2 = sum((r[2] / r[3]) ** 2 for r in ok); print(f'    chi2 of the timing residuals against a strictly periodic signal: {x2:.1f} for {len(ok)} seasons')
us = np.unique(s); h1 = np.isin(s, us[:len(us) // 2]); A1 = amp(t[h1], y[h1], e[h1], c[h1]); A2 = amp(t[~h1], y[~h1], e[~h1], c[~h1])
print(f'  first half of the sectors: {A1[0]:.2f} ± {A1[1]:.2f};  second half: {A2[0]:.2f} ± {A2[1]:.2f}  (difference {(A1[0]-A2[0])/np.hypot(A1[1],A2[1]):+.1f} sigma, white noise)')
for lab, mm in (('30-min', c > 0.015), ('10-min', (c > 0.004) & (c <= 0.015)), ('200-s', (c > 0.0018) & (c <= 0.004)), ('2-min', c <= 0.0018)):
    if mm.sum() > 500: A = amp(t[mm], y[mm], e[mm], c[mm]); print(f'    {lab} cadence data: {A[0]:.2f} ± {A[1]:.2f}')
json.dump(dict(tic=TIC, k=K, crowdsap=float(np.nanmean(crowd)), amp=a0, amp_err=ea, null_mean=float(nulls.mean()), null_std=float(nulls.std()), best_P=float(best), snr_scan=float(az.max()), seasons=[[float(x) for x in r] for r in rows], half1=A1[:2], half2=A2[:2]), open(f'{V}/{TIC}_{K}.indep' + ('_qlp' if QLP else '') + '.json', 'w'), default=float)
