"""Independent test of the TIC 366760231 signal in FFI sectors not used for detection (TESS-SPOC 58,59,85; QLP 57,77,78,79,84,86)."""
import os, glob, numpy as np, urllib.request, warnings, lightkurve as lk, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from astropy.timeseries import BoxLeastSquares
warnings.filterwarnings('ignore')
TIC = 366760231; t16 = f'{TIC:016d}'; D = f'/home/sarang/tess_search/deep/{TIC}'; os.makedirs(D, exist_ok=True)
def get(url, name):
    p = f'{D}/{name}'
    if not os.path.exists(p):
        try: urllib.request.urlretrieve(url, p)
        except Exception as e: print('fail', name, e); return None
    return p
def prep(lc):
    lc = lc.remove_nans().normalize(); lc = lc[np.isfinite(lc.flux.value)]; cad = np.nanmedian(np.diff(lc.time.value))
    fl = lc.flatten(window_length=int(round(.75 / cad)) | 1).remove_outliers(sigma_upper=3, sigma_lower=20)
    b = fl.bin(time_bin_size=10 / 1440).remove_nans() if cad < 9 / 1440 else fl
    return np.asarray(b.time.value, float), np.ma.filled(np.ma.asarray(b.flux.value, dtype=float), np.nan)
data = {}
for fn in sorted(glob.glob(f'{D}/*s_lc.fits')):
    lc = lk.read(fn, quality_bitmask='default'); data[f'2min-S{lc.meta["SECTOR"]}'] = prep(lc)
for s in (58, 59, 85):
    n = f'hlsp_tess-spoc_tess_phot_{t16}-s{s:04d}_tess_v1_lc.fits'
    p = get(f'https://archive.stsci.edu/hlsps/tess-spoc/s{s:04d}/target/{t16[0:4]}/{t16[4:8]}/{t16[8:12]}/{t16[12:16]}/{n}', n)
    if p: data[f'SPOC-FFI-S{s}'] = prep(lk.read(p, quality_bitmask='default'))
for s in (57, 77, 78, 79, 84, 86):
    n = f'hlsp_qlp_tess_ffi_s{s:04d}-{t16}_tess_v01_llc.fits'
    p = get(f'https://archive.stsci.edu/hlsps/qlp/s{s:04d}/{t16[0:4]}/{t16[4:8]}/{t16[8:12]}/{t16[12:16]}/{n}', n)
    if p:
        try: data[f'QLP-S{s}'] = prep(lk.read(p, quality_bitmask='default', flux_column='kspsap_flux'))
        except Exception as e: print('qlp read fail', s, e)
for k, (t, f) in data.items(): print(f'{k:14s} {len(t):5d} pts  {t.min():.1f}-{t.max():.1f}  rms {np.nanstd(f)*1e6:.0f} ppm')
old = [k for k in data if k.startswith('2min') and not k.endswith('S52')]
new = [k for k in data if not k.startswith('2min')]
def cat(keys):
    t = np.concatenate([data[k][0] for k in keys]); f = np.concatenate([data[k][1] for k in keys]); ok = np.isfinite(f); return t[ok], f[ok]
# 1) ephemeris from the detection data only (2-min S17-25)
t, f = cat(old); sig = 1.4826 * np.median(np.abs(f - np.median(f)))
bls = BoxLeastSquares(t, f, np.full(len(f), sig)); per = np.linspace(3.1360, 3.1387, 3001)
r = bls.power(per, [0.015, 0.02, 0.03, 0.045], objective='snr'); j = np.argmax(r.depth_snr)
P0, T00, dur = r.period[j], r.transit_time[j], r.duration[j]
# period uncertainty from the SNR peak width (where SNR drops by 1)
w = per[r.depth_snr > r.depth_snr[j] - 1]; sP = (w.max() - w.min()) / 2
print(f'\ndetection data only: P={P0:.6f} ± {sP:.6f}  T0={T00:.4f}  dur={dur*24:.2f} h  SNR={r.depth_snr[j]:.1f}')
# 2) test in NEW data only: scan the predicted transit time within its uncertainty
tn, fn_ = cat(new); sgn = 1.4826 * np.median(np.abs(fn_ - np.median(fn_)))
ncyc = np.median((tn - T00) / P0); st = np.hypot(r.transit_time[j] * 0 + 0.002, ncyc * sP)
def depth(tt, ff, P, T0, d):
    ph = ((tt - T0) / P + .5) % 1 - .5; a = np.abs(ph * P) < d / 2; b = (np.abs(ph * P) > d) & (np.abs(ph * P) < 4 * d)
    return (np.median(ff[b]) - ff[a].mean()) * 1e6, sgn * np.sqrt(1 / max(a.sum(), 1)) * 1e6, a.sum()
print(f'new data: {len(tn)} pts; propagated {ncyc:.0f} cycles -> timing uncertainty ±{st*24:.2f} h')
offs = np.linspace(-max(3 * st, 0.05), max(3 * st, 0.05), 121); Ps = [P0 + k * sP for k in (-1, -.5, 0, .5, 1)]
grid = np.array([[depth(tn, fn_, P, T00 + o, dur)[0] for o in offs] for P in Ps])
i, k = np.unravel_index(np.nanargmax(grid), grid.shape); dbest, ebest, nin = depth(tn, fn_, Ps[i], T00 + offs[k], dur)
print(f'best in new data: depth {dbest:.0f} ± {ebest:.0f} ppm ({nin} in-transit pts) at offset {offs[k]*24:+.2f} h, dP {Ps[i]-P0:+.6f}')
# null: same scan at random phases (look-elsewhere included)
rng = np.random.default_rng(0); null = []
for _ in range(300):
    sh = rng.uniform(.15, .85) * P0
    null.append(np.nanmax([[depth(tn, fn_, P, T00 + sh + o, dur)[0] for o in offs[::4]] for P in Ps[::2]]))
null = np.array(null); print(f'null (same scan at random phases): mean {null.mean():.0f}, std {null.std():.0f}; fraction >= observed {np.mean(null >= dbest):.3f}')
# 3) per-sector depths at the combined ephemeris
ta, fa = cat(list(data)); blsA = BoxLeastSquares(ta, fa, np.full(len(fa), sig)); perA = np.linspace(P0 - 3 * sP, P0 + 3 * sP, 4001)
rA = blsA.power(perA, [0.015, 0.02, 0.03, 0.045], objective='snr'); jA = np.argmax(rA.depth_snr); PA, TA, dA = rA.period[jA], rA.transit_time[jA], rA.duration[jA]
print(f'\nALL data: P={PA:.6f} T0={TA:.4f} dur={dA*24:.2f} h depth={rA.depth[jA]*1e6:.0f} SNR={rA.depth_snr[jA]:.1f}')
for kk, (tt, ff) in data.items():
    ok = np.isfinite(ff); dd, ee, nn = depth(tt[ok], ff[ok], PA, TA, dA); print(f'   {kk:14s} depth {dd:6.0f} ± {ee:4.0f} ppm  (n_in={nn})')
np.savez(f'{D}/ephem_all.npz', P=PA, t0=TA, dur=dA)
fig, ax = plt.subplots(1, 2, figsize=(13, 4))
for a, keys, lab in ((ax[0], old + ['2min-S52'], '2-min data'), (ax[1], new, 'NEW FFI sectors (not used for detection)')):
    tt, ff = cat([k for k in keys if k in data]); ph = (((tt - TA) / PA + .5) % 1 - .5) * PA * 24; m = np.abs(ph) < 4
    e = np.linspace(-4, 4, 41); idx = np.digitize(ph[m], e)
    a.plot(ph[m], ff[m], '.', color='.75', ms=2); a.errorbar((e[:-1] + e[1:]) / 2, [np.mean(ff[m][idx == q]) if np.sum(idx == q) > 5 else np.nan for q in range(1, len(e))],
               [np.std(ff[m][idx == q]) / np.sqrt(np.sum(idx == q)) if np.sum(idx == q) > 5 else np.nan for q in range(1, len(e))], fmt='ok')
    a.set_ylim(.996, 1.002); a.axvspan(-dA * 12, dA * 12, color='r', alpha=.15); a.set_title(lab); a.set_xlabel('hours from mid-transit')
plt.tight_layout(); plt.savefig('/home/sarang/tess_search/plots/366760231_newsectors.png', dpi=80)
