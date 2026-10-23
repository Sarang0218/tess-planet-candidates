import glob, numpy as np, warnings, json, lightkurve as lk
warnings.filterwarnings('ignore')
import triceratops.triceratops as tr
TIC = 366760231; D = '/home/sarang/tess_search/deep/366760231'
e = np.load(f'{D}/ephem_all.npz'); P, t0 = float(e['P']), float(e['t0'])
fit = json.load(open('/home/sarang/tess_search/cand_366760231/fit_results_freerho.json')); P = fit['period_d'][0]; t0 = fit['T0_BTJD'][0]
tpfs = sorted(glob.glob(f'{D}/*s_tp.fits')); sectors, aps = [], []
for fn in tpfs:
    tpf = lk.read(fn); sectors.append(tpf.sector); r, c = np.where(tpf.pipeline_mask)
    aps.append(np.array([[tpf.column + cc, tpf.row + rr] for rr, cc in zip(r, c)]))
print('sectors', sectors, 'aperture sizes', [len(a) for a in aps])
tg = tr.target(ID=TIC, sectors=np.array(sectors))
depth = fit['depth_ppm'][0] * 1e-6
tg.calc_depths(tdepth=depth, all_ap_pixels=aps)
print(tg.stars[['ID', 'Tmag', 'sep (arcsec)', 'fluxratio', 'tdepth']].head(12).to_string())
# phase-folded 2-min light curve near transit
T, F = [], []
for fn in sorted(glob.glob(f'{D}/*s_lc.fits')):
    lc = lk.read(fn, quality_bitmask='default').remove_nans().normalize(); cad = np.nanmedian(np.diff(lc.time.value))
    ph = ((lc.time.value - t0) / P + .5) % 1 - .5
    fl = lc.flatten(window_length=int(round(.75 / cad)) | 1, mask=np.abs(ph * P) < .05).remove_outliers(sigma_upper=4, sigma_lower=50)
    ph = (((fl.time.value - t0) / P + .5) % 1 - .5) * P; m = np.abs(ph) < 0.12
    T.append(ph[m]); F.append(np.ma.filled(np.ma.asarray(fl.flux.value[m], dtype=float), np.nan))
t = np.concatenate(T); f = np.concatenate(F); ok = np.isfinite(f); t, f = t[ok], f[ok]
o = np.argsort(t); t, f = t[o], f[o]
# bin to 4 minutes for speed
b = np.arange(t.min(), t.max(), 4 / 1440); idx = np.digitize(t, b)
tb = np.array([t[idx == k].mean() for k in np.unique(idx)]); fb = np.array([f[idx == k].mean() for k in np.unique(idx)])
eb = np.array([f[idx == k].std() / np.sqrt(np.sum(idx == k)) for k in np.unique(idx)])
print(f'folded LC: {len(tb)} bins, median err {np.median(eb)*1e6:.0f} ppm')
tg.calc_probs(time=tb, flux_0=fb, flux_err_0=np.median(eb), P_orb=P, N=300000, parallel=True)
df = tg.probs; print(df[['ID', 'scenario', 'prob']].sort_values('prob', ascending=False).head(12).to_string())
print(f'\nFPP = {tg.FPP:.4f}   NFPP = {tg.NFPP:.4f}')
import os
run = os.environ.get('RUN', '0'); json.dump(dict(FPP=float(tg.FPP), NFPP=float(tg.NFPP)), open(f'/home/sarang/tess_search/cand_366760231/triceratops_{run}.json', 'w'))
df.to_csv('/home/sarang/tess_search/cand_366760231/triceratops_probs.csv', index=False)
