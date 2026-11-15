"""Vet the 1.1434-d signal on TIC 366760231 from the raw light curves, with the 3.1373-d transits masked and the
detrending window varied (with/without masking the signal itself)."""
import glob, numpy as np, warnings, lightkurve as lk, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')
D = '/home/sarang/tess_search/deep/366760231'
Pb, Tb, Db = 3.1373207, 1962.04155, 0.6 / 24          # known candidate (BTJD)
P2, T2, D2 = 1.143413, 1999.7835, 1.08 / 24
files = sorted(glob.glob(f'{D}/*s_lc.fits')) + sorted(glob.glob(f'{D}/hlsp_tess-spoc*_lc.fits'))
def load(win, mask2):
    out = {}
    for fn in files:
        lc = lk.read(fn, quality_bitmask='default', flux_column='pdcsap_flux').remove_nans().normalize(); lc = lc[np.isfinite(lc.flux.value)]
        t = lc.time.value; cad = np.nanmedian(np.diff(t)); sec = lc.meta['SECTOR']
        mb = np.abs(((t - Tb) / Pb + .5) % 1 - .5) * Pb < 1.0 * Db; m2 = np.abs(((t - T2) / P2 + .5) % 1 - .5) * P2 < 0.8 * D2
        fl = lc.flatten(window_length=int(round(win / cad)) | 1, mask=(mb | m2) if mask2 else mb).remove_outliers(sigma_upper=3, sigma_lower=20)
        tt = np.asarray(fl.time.value, float); keep = ~(np.abs(((tt - Tb) / Pb + .5) % 1 - .5) * Pb < 1.0 * Db); fl = fl[keep]     # remove candidate-b transits
        if cad < 9 / 1440: fl = fl.bin(time_bin_size=10 / 1440).remove_nans()
        tt = np.asarray(fl.time.value, float); ff = np.ma.filled(np.ma.asarray(fl.flux.value, dtype=float), np.nan); ok = np.isfinite(ff)
        out[sec] = (tt[ok], ff[ok])
    return out
def depth(t, f, P, T0, d, shift=0.0):
    ph = (((t - T0) / P + .5 + shift) % 1 - .5) * P; a = np.abs(ph) < d / 2; b = (np.abs(ph) > 0.75 * d) & (np.abs(ph) < 3 * d)
    if a.sum() < 3 or b.sum() < 5: return np.nan, np.nan
    sg = 1.4826 * np.median(np.abs(f - np.median(f))); return (f[b].mean() - f[a].mean()) * 1e6, sg * np.sqrt(1 / a.sum() + 1 / b.sum()) * 1e6
for win in (0.5, 0.75, 1.5):
    for mask2 in (False, True):
        data = load(win, mask2); ds = {s: depth(*data[s], P2, T2, D2) for s in data}
        v = np.array([ds[s] for s in sorted(ds)]); w = 1 / v[:, 1] ** 2; dm = np.sum(v[:, 0] * w) / w.sum(); de = 1 / np.sqrt(w.sum())
        print(f'window {win} d, signal masked in detrend={mask2}: depth {dm:.0f} ± {de:.0f} ppm ({dm/de:.1f} sigma); chi2 {np.sum(w*(v[:,0]-dm)**2):.1f}/{len(v)-1}; per sector: ' + ' '.join(f'S{s}:{ds[s][0]:.0f}±{ds[s][1]:.0f}' for s in sorted(ds)))
data = load(0.75, True); t = np.concatenate([data[s][0] for s in sorted(data)]); f = np.concatenate([data[s][1] for s in sorted(data)])
n = np.floor((t - T2) / P2 + .5).astype(int)
for lab, m in (('odd', n % 2 == 1), ('even', n % 2 == 0)): print(lab, '%.0f ± %.0f' % depth(t[m], f[m], P2, T2, D2))
print('secondary (phase 0.5): %.0f ± %.0f' % depth(t, f, P2, T2, D2, 0.5))
print('at 2P, phases 0 / 0.5: %.0f ± %.0f' % depth(t, f, 2 * P2, T2, D2), '/ %.0f ± %.0f' % depth(t, f, 2 * P2, T2 + P2, D2))
# phase scan of the same box: how often does a random phase give this depth?
sh = np.linspace(0.05, 0.95, 181); dd = np.array([depth(t, f, P2, T2, D2, s)[0] for s in sh]); d0, e0 = depth(t, f, P2, T2, D2)
print(f'depth at T0 {d0:.0f} ± {e0:.0f}; other phases: mean {np.nanmean(dd):.0f}, std {np.nanstd(dd):.0f}, max {np.nanmax(dd):.0f}, min {np.nanmin(dd):.0f}')
fig, ax = plt.subplots(1, 3, figsize=(16, 4.2))
ph = (((t - T2) / P2 + .5) % 1 - .5) * P2 * 24
for a, lim, nb in ((ax[0], 13.7, 60), (ax[1], 4, 33)):
    m = np.abs(ph) < lim; e = np.linspace(-lim, lim, nb + 1); idx = np.digitize(ph[m], e)
    a.plot(ph[m], (f[m] - 1) * 1e6, '.', color='.85', ms=1)
    a.errorbar((e[:-1] + e[1:]) / 2, [(f[m][idx == q].mean() - 1) * 1e6 for q in range(1, nb + 1)], [f[m][idx == q].std() / np.sqrt((idx == q).sum()) * 1e6 for q in range(1, nb + 1)], fmt='ok', ms=4)
    a.axvspan(-D2 * 12, D2 * 12, color='r', alpha=.12); a.set_ylim(-1500, 1000); a.set_xlabel('hours from mid-transit'); a.set_ylabel('ppm')
ax[0].set_title('TIC 366760231, P = 1.1434 d fold (b transits removed)'); ax[1].set_title('zoom')
ss = sorted(data); v = np.array([depth(*data[s], P2, T2, D2) for s in ss]); ax[2].errorbar(range(len(ss)), v[:, 0], v[:, 1], fmt='ok'); ax[2].axhline(0, color='k', lw=.5); ax[2].axhline(d0, color='r', ls='--')
ax[2].set_xticks(range(len(ss))); ax[2].set_xticklabels([f'S{s}' for s in ss]); ax[2].set_title('per-sector depth (ppm)')
plt.tight_layout(); plt.savefig('/home/sarang/tess_search/plots/second_366760231.png', dpi=90)
