"""Deep vetting for FFI-only candidates (TESS-SPOC FFI light curves): independent test in sectors not used for detection,
per-sector depths, odd/even, secondary. Usage: ffi_deep.py TIC P t0 dur_days"""
import sys, os, json, numpy as np, pandas as pd, urllib.request, warnings, lightkurve as lk
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from astropy.timeseries import BoxLeastSquares
warnings.filterwarnings('ignore')
TIC, P0, T00, DUR = int(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])
D = f'/home/sarang/tess_search/deep/{TIC}'; os.makedirs(D, exist_ok=True); t16 = f'{TIC:016d}'
meta = json.load(open(f'/home/sarang/tess_search/out3/{TIC}.meta.json')); used = set(meta['used']); allsec = meta['all_sectors']
data = {}
for s in allsec:
    n = f'hlsp_tess-spoc_tess_phot_{t16}-s{s:04d}_tess_v1_lc.fits'; p = f'{D}/{n}'
    if not os.path.exists(p):
        try: urllib.request.urlretrieve(f'https://archive.stsci.edu/hlsps/tess-spoc/s{s:04d}/target/{t16[0:4]}/{t16[4:8]}/{t16[8:12]}/{t16[12:16]}/{n}', p)
        except Exception as e: print('fail', s, e); continue
    try: lc = lk.read(p, quality_bitmask='default').remove_nans().normalize()
    except Exception as e: print('read fail', s, e); continue
    lc = lc[np.isfinite(lc.flux.value)]; cad = np.nanmedian(np.diff(lc.time.value))
    ph = ((lc.time.value - T00) / P0 + .5) % 1 - .5
    fl = lc.flatten(window_length=int(round(.75 / cad)) | 1, mask=np.abs(ph * P0) < 1.5 * DUR).remove_outliers(sigma_upper=3, sigma_lower=20)
    data[s] = (np.asarray(fl.time.value, float), np.ma.filled(np.ma.asarray(fl.flux.value, dtype=float), np.nan), cad)
print(f'TIC {TIC}: {len(data)} sectors loaded; detection sectors {sorted(used)}')
def cat(keys):
    t = np.concatenate([data[k][0] for k in keys]); f = np.concatenate([data[k][1] for k in keys]); ok = np.isfinite(f); return t[ok], f[ok]
det = [s for s in data if s in used]; new = [s for s in data if s not in used]
t, f = cat(det); sig = 1.4826 * np.median(np.abs(f - np.median(f)))
bls = BoxLeastSquares(t, f, np.full(len(f), sig)); dP = 3 * P0 ** 2 / np.ptp(t); per = np.linspace(P0 - dP, P0 + dP, 4001)
durs = [DUR * .7, DUR, DUR * 1.4]; r = bls.power(per, durs, objective='snr'); j = np.argmax(r.depth_snr)
P1, T1, d1 = r.period[j], r.transit_time[j], r.duration[j]; w = per[r.depth_snr > r.depth_snr[j] - 1]; sP = max((w.max() - w.min()) / 2, 1e-5)
print(f'detection-only ephemeris: P={P1:.5f} ± {sP:.5f}, T0={T1:.4f}, dur={d1*24:.2f} h, depth={r.depth[j]*1e6:.0f} ppm, SNR={r.depth_snr[j]:.1f}')
tn, fn_ = cat(new); sgn = 1.4826 * np.median(np.abs(fn_ - np.median(fn_)))
def depth(tt, ff, P, T0, d, s_=None):
    ph = ((tt - T0) / P + .5) % 1 - .5; a = np.abs(ph * P) < d / 2; b = (np.abs(ph * P) > d) & (np.abs(ph * P) < 4 * d)
    if a.sum() < 2 or b.sum() < 5: return np.nan, np.nan, 0
    return (np.median(ff[b]) - ff[a].mean()) * 1e6, (s_ or sgn) * np.sqrt(1 / a.sum()) * 1e6, a.sum()
ncyc = np.abs(np.median((tn - T1) / P1)); st = max(ncyc * sP, 0.01)
offs = np.linspace(-3 * st, 3 * st, 61); Ps = [P1 + k * sP for k in np.linspace(-2, 2, 9)]
def scan(T0base):
    best = (-1e9, 0, 0)
    for P in Ps:
        for o in offs:
            dd = depth(tn, fn_, P, T0base + o, d1)[0]
            if np.isfinite(dd) and dd > best[0]: best = (dd, P, o)
    return best
dbest, Pb, ob = scan(T1); dd, ee, nin = depth(tn, fn_, Pb, T1 + ob, d1)
print(f'NEW sectors {sorted(new)}: {len(tn)} pts; timing window ±{3*st*24:.1f} h; best depth {dd:.0f} ± {ee:.0f} ppm (n_in={nin}) at offset {ob*24:+.2f} h')
rng = np.random.default_rng(0); null = np.array([scan(T1 + rng.uniform(.15, .85) * P1)[0] for _ in range(150)])
print(f'null (identical scan at random phases): mean {np.nanmean(null):.0f} ± {np.nanstd(null):.0f} ppm; fraction >= observed: {np.mean(null >= dbest):.3f}; '
      f'significance ~{(dbest - np.nanmean(null)) / np.nanstd(null):.1f} sigma')
ta, fa = cat(list(data)); blsA = BoxLeastSquares(ta, fa, np.full(len(fa), sig))
perA = np.linspace(Pb - 4 * sP, Pb + 4 * sP, 4001); rA = blsA.power(perA, durs, objective='snr'); jA = np.argmax(rA.depth_snr)
PA, TA, dA = rA.period[jA], rA.transit_time[jA], rA.duration[jA]
print(f'\nALL sectors: P={PA:.6f} T0={TA:.4f} dur={dA*24:.2f} h depth={rA.depth[jA]*1e6:.0f} ppm SNR={rA.depth_snr[jA]:.1f}')
st_ = blsA.compute_stats(PA, dA, TA)
print(f'   odd {st_["depth_odd"][0]*1e6:.0f} ± {st_["depth_odd"][1]*1e6:.0f}   even {st_["depth_even"][0]*1e6:.0f} ± {st_["depth_even"][1]*1e6:.0f}   half-phase {st_["depth_half"][0]*1e6:.0f} ± {st_["depth_half"][1]*1e6:.0f}')
for s in sorted(data):
    tt, ff, cad = data[s]; ok = np.isfinite(ff); s_ = 1.4826 * np.median(np.abs(ff[ok] - np.median(ff[ok])))
    dd, ee, nn = depth(tt[ok], ff[ok], PA, TA, dA, s_); print(f'   S{s:3d} {"(det)" if s in used else "     "} cad {cad*1440:4.0f} min  depth {dd:7.0f} ± {ee:5.0f} ppm  n_in={nn}')
np.savez(f'{D}/ephem_ffi.npz', P=PA, t0=TA, dur=dA)
fig, ax = plt.subplots(1, 2, figsize=(13, 4))
for a, keys, lab in ((ax[0], det, f'detection sectors {sorted(det)}'), (ax[1], new, 'NEW sectors (not used for detection)')):
    tt, ff = cat(keys); ph = (((tt - TA) / PA + .5) % 1 - .5) * PA * 24; m = np.abs(ph) < 5
    e = np.linspace(-5, 5, 41); idx = np.digitize(ph[m], e)
    a.plot(ph[m], ff[m], '.', color='.8', ms=2)
    a.errorbar((e[:-1] + e[1:]) / 2, [np.mean(ff[m][idx == q]) if np.sum(idx == q) > 5 else np.nan for q in range(1, len(e))],
               [np.std(ff[m][idx == q]) / np.sqrt(np.sum(idx == q)) if np.sum(idx == q) > 5 else np.nan for q in range(1, len(e))], fmt='ok')
    a.axvspan(-dA * 12, dA * 12, color='r', alpha=.15); a.set_title(lab, fontsize=9); a.set_xlabel('hours from mid-transit')
    a.set_ylim(1 - 2.5 * rA.depth[jA], 1 + rA.depth[jA])
plt.tight_layout(); plt.savefig(f'/home/sarang/tess_search/plots/ffideep_{TIC}.png', dpi=80)
