"""Deep vetting of one candidate: all 2-min sectors, ephemeris refinement with alias scan, per-sector depths
(transit-to-transit errors) at phase 0 and 0.5, and difference-image centroids from target pixel files.
Usage: deep_vet.py TIC P t0 dur_days [max_tpf_sectors]"""
import sys, os, numpy as np, pandas as pd, urllib.request, warnings, lightkurve as lk
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from astropy.timeseries import BoxLeastSquares
warnings.filterwarnings('ignore')
TIC, P0, T00, DUR0 = int(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])
NTPF = int(sys.argv[5]) if len(sys.argv) > 5 else 3
D = f'/home/sarang/tess_search/deep/{TIC}'; os.makedirs(D, exist_ok=True)
S3 = 'https://stpubdata.s3.amazonaws.com/tess/public/tid'
av = sorted(pd.read_pickle('/home/sarang/tess_search/data/avail2min.pkl')[TIC])

def fetch(s, fn):
    p = f'{D}/{fn}'
    if not os.path.exists(p):
        t16 = fn.split('-')[2]; urllib.request.urlretrieve(f'{S3}/s{s:04d}/{t16[0:4]}/{t16[4:8]}/{t16[8:12]}/{t16[12:16]}/{fn}', p)
    return p

parts = {}
for s, fn in av:
    try: lc = lk.read(fetch(s, fn), quality_bitmask='default').remove_nans().normalize()
    except Exception as e: print('skip', s, e); continue
    cad = np.nanmedian(np.diff(lc.time.value))
    fl = lc.flatten(window_length=int(round(0.75 / cad)) | 1).remove_outliers(sigma_upper=3, sigma_lower=20).bin(time_bin_size=10 / 1440).remove_nans()
    parts[s] = (fl.time.value, fl.flux.value)
t = np.concatenate([v[0] for v in parts.values()]); f = np.concatenate([v[1] for v in parts.values()])
print(f'TIC {TIC}: sectors {list(parts)}  baseline {np.ptp(t):.0f} d')

# --- ephemeris: scan a window covering cycle-count aliases across the long gaps
sig = 1.4826 * np.median(np.abs(f - np.median(f)))
bls = BoxLeastSquares(t, f, np.full(len(f), sig))
dP = max(0.003, 3 * P0 ** 2 / np.ptp(t))
periods = np.linspace(P0 - dP, P0 + dP, 20001)
r = bls.power(periods, [DUR0 * .7, DUR0, DUR0 * 1.4], objective='snr')
order = np.argsort(r.depth_snr)[::-1]; peaks = []
for k in order:
    if all(abs(r.period[k] - x[0]) > P0 ** 2 / np.ptp(t) / 3 for x in peaks): peaks.append((r.period[k], r.depth_snr[k], r.transit_time[k], r.duration[k], r.depth[k]))
    if len(peaks) == 5: break
print('alias scan (period, SNR):', [(round(p, 5), round(s, 1)) for p, s, *_ in peaks])
P, snr, t0, dur, dep = peaks[0]
print(f'best: P={P:.6f} t0={t0:.4f} dur={dur*24:.2f} h depth={dep*1e6:.0f} ppm SNR={snr:.1f}  (runner-up alias SNR {peaks[1][1]:.1f})')

def per_transit(ts, fs, off):
    ph = ((ts - t0) / P - off + .5) % 1 - .5; n = np.floor((ts - t0) / P - off + .5).astype(int); out = []
    for k in np.unique(n):
        a = (n == k) & (np.abs(ph * P) < dur / 2); b = (n == k) & (np.abs(ph * P) > dur) & (np.abs(ph * P) < 3 * dur)
        if a.sum() >= 3 and b.sum() >= 8: out.append(np.median(fs[b]) - fs[a].mean())
    return np.array(out)

print('\nper-sector depths (ppm; error = transit-to-transit scatter / sqrt(n)):')
allp, alls = [], []
for s, (ts, fs) in parts.items():
    a0, a5 = per_transit(ts, fs, 0), per_transit(ts, fs, .5); allp += list(a0); alls += list(a5)
    e = lambda a: a.std(ddof=1) / np.sqrt(len(a)) * 1e6 if len(a) > 1 else np.nan
    print(f'  S{s:3d}: primary {a0.mean()*1e6 if len(a0) else np.nan:6.0f} ± {e(a0):4.0f} (n={len(a0):2d})   half-phase {a5.mean()*1e6 if len(a5) else np.nan:6.0f} ± {e(a5):4.0f} (n={len(a5):2d})')
allp, alls = np.array(allp), np.array(alls)
chi = np.sum((allp - allp.mean()) ** 2) / (allp.std(ddof=1) ** 2) if len(allp) > 2 else np.nan
print(f'ALL: primary {allp.mean()*1e6:.0f} ± {allp.std(ddof=1)/np.sqrt(len(allp))*1e6:.0f} ppm (n={len(allp)}), '
      f'half-phase {alls.mean()*1e6:.0f} ± {alls.std(ddof=1)/np.sqrt(len(alls))*1e6:.0f} ppm (n={len(alls)})')

# --- difference images in up to NTPF sectors (prefer those with most transits)
cnt = {s: len(per_transit(*parts[s], 0)) for s in parts}
use = sorted(cnt, key=lambda s: -cnt[s])[:NTPF]
fig, ax = plt.subplots(2, len(use), figsize=(5 * len(use), 9), squeeze=False)
for col, s in enumerate(use):
    fn = dict(av)[s].replace('s_lc.fits', 's_tp.fits')
    try: tpf = lk.read(fetch(s, fn), quality_bitmask='default')
    except Exception as e: print('tpf skip', s, e); continue
    tt = tpf.time.value; fl = tpf.flux.value; ok = np.all(np.isfinite(fl.reshape(len(tt), -1)), 1) & np.isfinite(tt); tt, fl = tt[ok], fl[ok]
    ph = (((tt - t0) / P + .5) % 1 - .5) * P; n = np.floor((tt - t0) / P + .5).astype(int)
    intr, near = np.abs(ph) < dur / 3, (np.abs(ph) > .75 * dur) & (np.abs(ph) < 2 * dur)
    diffs = [np.nanmedian(fl[near & (n == k)], 0) - np.nanmean(fl[intr & (n == k)], 0) for k in np.unique(n[intr]) if (intr & (n == k)).sum() > 3 and (near & (n == k)).sum() > 8]
    diff = np.nanmean(diffs, 0); oot = np.nanmedian(fl[near], 0)
    yy, xx = np.indices(diff.shape)
    # star position: flux-weighted centroid of the 3x3 around the brightest OOT pixel (robust to proper motion)
    by, bx = np.unravel_index(np.nanargmax(oot), oot.shape); m3 = (abs(yy - by) <= 1) & (abs(xx - bx) <= 1)
    sx, sy = (oot * xx)[m3].sum() / oot[m3].sum(), (oot * yy)[m3].sum() / oot[m3].sum()
    dy, dx = np.unravel_index(np.nanargmax(diff), diff.shape); m3d = (abs(yy - dy) <= 1) & (abs(xx - dx) <= 1); w = np.clip(diff, 0, None)
    cx, cy = (w * xx)[m3d].sum() / w[m3d].sum(), (w * yy)[m3d].sum() / w[m3d].sum()
    snr_img = np.nanmax(diff) / (np.nanstd(np.array(diffs), 0)[dy, dx] / np.sqrt(len(diffs)))
    print(f'  S{s}: {len(diffs)} transits; star centroid ({sx:.2f},{sy:.2f}); difference centroid ({cx:.2f},{cy:.2f}); '
          f'offset {np.hypot(cx - sx, cy - sy):.2f} px = {np.hypot(cx - sx, cy - sy) * 21:.0f}"; diff-peak S/N {snr_img:.1f}')
    ax[0, col].imshow(oot, origin='lower'); ax[0, col].plot(sx, sy, 'rx', ms=12); ax[0, col].set_title(f'S{s} out-of-transit')
    ax[1, col].imshow(diff, origin='lower', cmap='RdBu_r', vmin=-np.nanmax(abs(diff)), vmax=np.nanmax(abs(diff)))
    ax[1, col].plot(sx, sy, 'kx', ms=12); ax[1, col].plot(cx, cy, 'g+', ms=14, mew=2); ax[1, col].contour(tpf.pipeline_mask, levels=[.5], colors='y')
    ax[1, col].set_title(f'S{s} difference (out - in), n={len(diffs)}', fontsize=9)
plt.tight_layout(); plt.savefig(f'/home/sarang/tess_search/plots/deep_{TIC}_diffimg.png', dpi=60)

# --- folded plots per sector
fig, ax = plt.subplots(len(parts), 2, figsize=(12, 2.6 * len(parts)), squeeze=False)
for i, (s, (ts, fs)) in enumerate(parts.items()):
    for j, off in enumerate((0, .5)):
        a = ax[i, j]; ph = (((ts - t0) / P - off + .5) % 1 - .5) * P * 24; m = np.abs(ph) < 4 * dur * 24
        a.plot(ph[m], fs[m], '.', color='.7', ms=2); e = np.linspace(-4 * dur * 24, 4 * dur * 24, 33); idx = np.digitize(ph[m], e)
        a.plot((e[:-1] + e[1:]) / 2, [np.median(fs[m][idx == k]) if np.sum(idx == k) else np.nan for k in range(1, len(e))], 'o-r', ms=3)
        a.axhline(1 - dep, ls='--'); a.set_title(f'S{s} phase {off}', fontsize=8)
plt.tight_layout(); plt.savefig(f'/home/sarang/tess_search/plots/deep_{TIC}_sectors.png', dpi=60)
np.savez(f'{D}/ephem.npz', P=P, t0=t0, dur=dur, dep=dep)
