"""TRICERATOPS for the 1.1434-d signal on TIC 366760231 (candidate-.01 transits removed). Env RUN = run index."""
import glob, os, numpy as np, warnings, json, lightkurve as lk
warnings.filterwarnings('ignore')
import triceratops.triceratops as tr
TIC = 366760231; D = '/home/sarang/tess_search/deep/366760231'; O = '/home/sarang/tess_search/cand_366760231/two_planet'
eph = json.load(open(f'{O}/ephemeris2.json')); fit = json.load(open(f'{O}/fit2_results.json'))
P, t0 = eph['c']['P'], eph['c']['T0_BTJD'] + fit['c']['dt_min'][0] / 1440; Pb, Tb = eph['b']['P'], eph['b']['T0_BTJD']
tpfs = sorted(glob.glob(f'{D}/*s_tp.fits')); sectors, aps = [], []
for fn in tpfs:
    tpf = lk.read(fn); sectors.append(tpf.sector); r, c = np.where(tpf.pipeline_mask)
    aps.append(np.array([[tpf.column + cc, tpf.row + rr] for rr, cc in zip(r, c)]))
tg = tr.target(ID=TIC, sectors=np.array(sectors)); depth = fit['c']['depth_ppm'][0] * 1e-6
tg.calc_depths(tdepth=depth, all_ap_pixels=aps)
print(tg.stars[['ID', 'Tmag', 'sep (arcsec)', 'fluxratio', 'tdepth']].head(10).to_string())
T, F = [], []
for fn in sorted(glob.glob(f'{D}/*s_lc.fits')):
    lc = lk.read(fn, quality_bitmask='default').remove_nans().normalize(); t = lc.time.value; cad = np.nanmedian(np.diff(t))
    m = (np.abs(((t - t0) / P + .5) % 1 - .5) * P < .04) | (np.abs(((t - Tb) / Pb + .5) % 1 - .5) * Pb < .03)
    fl = lc.flatten(window_length=int(round(.75 / cad)) | 1, mask=m).remove_outliers(sigma_upper=4, sigma_lower=50)
    tt = np.asarray(fl.time.value, float); ff = np.ma.filled(np.ma.asarray(fl.flux.value, dtype=float), np.nan)
    keep = ~(np.abs(((tt - Tb) / Pb + .5) % 1 - .5) * Pb < .025); ph = (((tt - t0) / P + .5) % 1 - .5) * P; m = (np.abs(ph) < 0.16) & keep
    T.append(ph[m]); F.append(ff[m])
t = np.concatenate(T); f = np.concatenate(F); ok = np.isfinite(f); t, f = t[ok], f[ok]
b = np.arange(t.min(), t.max(), 4 / 1440); idx = np.digitize(t, b); u = np.unique(idx)
tb = np.array([t[idx == k].mean() for k in u]); fb = np.array([f[idx == k].mean() for k in u]); eb = np.array([f[idx == k].std() / np.sqrt(np.sum(idx == k)) for k in u])
print(f'folded LC: {len(tb)} bins, median err {np.median(eb)*1e6:.0f} ppm')
tg.calc_probs(time=tb, flux_0=fb, flux_err_0=np.median(eb), P_orb=P, N=300000, parallel=True)
df = tg.probs; print(df[['ID', 'scenario', 'prob']].sort_values('prob', ascending=False).head(10).to_string())
print(f'\nFPP = {tg.FPP:.4f}   NFPP = {tg.NFPP:.4f}')
run = os.environ.get('RUN', '0'); json.dump(dict(FPP=float(tg.FPP), NFPP=float(tg.NFPP)), open(f'{O}/triceratops_c_{run}.json', 'w'))
df.to_csv(f'{O}/triceratops_c_probs_{run}.csv', index=False)
