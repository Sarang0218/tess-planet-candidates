"""Check individual transit-like events: SAP vs PDCSAP, background, quality flags/momentum dumps, centroids, difference image.
Usage: single_event_check.py TIC t1 [t2 ...] dur_h"""
import sys, os, numpy as np, pandas as pd, urllib.request, warnings, lightkurve as lk
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')
TIC = int(sys.argv[1]); times = [float(x) for x in sys.argv[2:-1]]; dur = float(sys.argv[-1]) / 24
D = f'/home/sarang/tess_search/deep/{TIC}'; os.makedirs(D, exist_ok=True); S3 = 'https://stpubdata.s3.amazonaws.com/tess/public/tid'
av = sorted(pd.read_pickle('/home/sarang/tess_search/data/avail2min.pkl')[TIC]); print('2-min sectors:', [s for s, _ in av])
def fetch(s, fn):
    p = f'{D}/{fn}'
    if not os.path.exists(p):
        t16 = fn.split('-')[2]; urllib.request.urlretrieve(f'{S3}/s{s:04d}/{t16[0:4]}/{t16[4:8]}/{t16[8:12]}/{t16[12:16]}/{fn}', p)
    return p
lcs = {}
for s, fn in av:
    try: lcs[s] = lk.read(fetch(s, fn), quality_bitmask='none')
    except Exception as e: print('skip', s, e)
fig, ax = plt.subplots(len(times), 4, figsize=(18, 3.2 * len(times)), squeeze=False)
for i, tc in enumerate(times):
    s = [k for k, lc in lcs.items() if lc.time.value.min() - .1 < tc < lc.time.value.max() + .1]
    if not s: print(f'event {tc}: no 2-min sector'); continue
    lc = lcs[s[0]]; t = lc.time.value; m = np.abs(t - tc) < 1.0
    q = lc.quality.value[m]; tt = t[m]
    pd_ = lc.pdcsap_flux.value[m] / np.nanmedian(lc.pdcsap_flux.value[m]); sap = lc.sap_flux.value[m] / np.nanmedian(lc.sap_flux.value[m])
    bkg = lc.sap_bkg.value[m]; cx, cy = lc.centroid_col.value[m], lc.centroid_row.value[m]
    intr = np.abs(tt - tc) < dur / 2; oot = (np.abs(tt - tc) > dur) & (np.abs(tt - tc) < 3 * dur)
    flagged = np.unique(q[intr & (q > 0)]); dumps = tt[(q & 32) > 0]
    near_dump = np.min(np.abs(dumps - tc)) * 24 if len(dumps) else np.inf
    def dd(x): return (np.nanmedian(x[oot]) - np.nanmean(x[intr])) * 1e6
    print(f'event {tc:.3f} (S{s[0]}): PDCSAP depth {dd(pd_):.0f} ppm, SAP depth {dd(sap):.0f} ppm; bkg in/out {np.nanmean(bkg[intr]):.1f}/{np.nanmedian(bkg[oot]):.1f}; '
          f'centroid shift col {(np.nanmean(cx[intr])-np.nanmedian(cx[oot]))*1e3:.2f} mpx, row {(np.nanmean(cy[intr])-np.nanmedian(cy[oot]))*1e3:.2f} mpx; '
          f'quality flags in transit: {list(flagged)}; nearest momentum dump {near_dump:.1f} h; in-transit points {intr.sum()}')
    for j, (y, lab) in enumerate(((pd_, 'PDCSAP'), (sap, 'SAP'), (bkg, 'background'))):
        a = ax[i, j]; a.plot((tt - tc) * 24, y, '.', ms=1.5, color='.6'); k = np.isfinite(y)
        b = np.linspace(-24, 24, 97); idx = np.digitize((tt[k] - tc) * 24, b)
        a.plot((b[:-1] + b[1:]) / 2, [np.nanmedian(y[k][idx == u]) if np.sum(idx == u) > 2 else np.nan for u in range(1, len(b))], 'r-')
        a.axvspan(-dur * 12, dur * 12, color='y', alpha=.3); a.set_title(f'S{s[0]} event {tc:.2f} {lab}', fontsize=8)
        for x in dumps: a.axvline((x - tc) * 24, c='b', ls=':', lw=.8)
    ax[i, 3].plot((tt - tc) * 24, cx - np.nanmedian(cx), '.', ms=1, label='col'); ax[i, 3].plot((tt - tc) * 24, cy - np.nanmedian(cy), '.', ms=1, label='row')
    ax[i, 3].axvspan(-dur * 12, dur * 12, color='y', alpha=.3); ax[i, 3].set_ylim(-.05, .05); ax[i, 3].legend(fontsize=6); ax[i, 3].set_title('flux-weighted centroid (px)', fontsize=8)
plt.tight_layout(); plt.savefig(f'/home/sarang/tess_search/plots/events_{TIC}.png', dpi=60)
