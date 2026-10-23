"""Find 10-min time bins where many stars dip together (common-mode systematics)."""
import glob, numpy as np
fs = sorted(glob.glob('/home/sarang/tess_search/out/*.npz'))
W = 10 / 1440; tot = {}; low = {}
for fn in fs:
    d = np.load(fn); t, f = d['t'], d['f']
    z = (f - np.median(f)) / (1.4826 * np.median(np.abs(f - np.median(f))))
    b = np.floor(t / W).astype(np.int64)
    for arr, store in ((b, tot), (b[z < -3], low)):
        u, c = np.unique(arr, return_counts=True)
        for k, v in zip(u, c): store[k] = store.get(k, 0) + 1
keys = np.array(sorted(tot)); n = np.array([tot[k] for k in keys]); l = np.array([low.get(k, 0) for k in keys])
frac = l / n
base = np.median(frac[n >= 20])
bad = keys[(n >= 20) & (l >= 4) & (frac > max(5 * base, 0.03))]
bad = np.unique(np.concatenate([bad - 1, bad, bad + 1]))
np.save('/home/sarang/tess_search/data/bad_bins.npy', bad)
print(f'{len(fs)} stars, median low-frac {base:.4f}, bad bins {len(bad)} = {len(bad)*10/1440:.1f} d of {len(keys)*10/1440:.0f} d')
t = keys[np.isin(keys, bad)] * W
print('bad-time clusters (BTJD):', np.round(t[np.r_[True, np.diff(t) > 0.5]], 1)[:40])
