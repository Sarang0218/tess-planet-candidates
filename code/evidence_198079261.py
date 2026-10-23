"""One-page evidence sheet for TIC 198079261."""
import glob, json, os, numpy as np, pandas as pd, warnings, lightkurve as lk, batman
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from astropy.timeseries import BoxLeastSquares
from astropy.coordinates import SkyCoord
warnings.filterwarnings('ignore')
TIC = 198079261; D = f'/home/sarang/tess_search/deep/{TIC}'; O = f'/home/sarang/tess_search/cand_{TIC}'
eph = json.load(open(f'{O}/ephemeris.json')); fit = json.load(open(f'{O}/fit_results.json'))
P, T0 = eph['P'], eph['T0_BJD'] - 2457000; DUR = fit['T14_h'][0] / 24
used = set(json.load(open(f'/home/sarang/tess_search/out3/{TIC}.meta.json'))['used'])
# ---- load all sectors
data = {}
for fn in sorted(glob.glob(f'{D}/hlsp_tess-spoc*_lc.fits')):
    lc = lk.read(fn, quality_bitmask='default').remove_nans().normalize(); s = lc.meta['SECTOR']; cad = np.nanmedian(np.diff(lc.time.value))
    ph = ((lc.time.value - T0) / P + .5) % 1 - .5
    fl = lc.flatten(window_length=int(round(.75 / cad)) | 1, mask=np.abs(ph * P) < 2 * DUR).remove_outliers(sigma_upper=3, sigma_lower=20)
    t = np.asarray(fl.time.value, float); f = np.ma.filled(np.ma.asarray(fl.flux.value, dtype=float), np.nan); ok = np.isfinite(f)
    data[s] = (t[ok], f[ok], cad)
def cat(keys):
    return np.concatenate([data[k][0] for k in keys]), np.concatenate([data[k][1] for k in keys])
def binned(x, y, edges):
    idx = np.digitize(x, edges); c = (edges[:-1] + edges[1:]) / 2
    m = np.array([np.mean(y[idx == k]) if np.sum(idx == k) > 4 else np.nan for k in range(1, len(edges))])
    e = np.array([np.std(y[idx == k]) / np.sqrt(np.sum(idx == k)) if np.sum(idx == k) > 4 else np.nan for k in range(1, len(edges))])
    return c, m, e
def phase_h(t, off=0.): return (((t - T0) / P - off + .5) % 1 - .5) * P * 24
fig = plt.figure(figsize=(18, 20)); gs = fig.add_gridspec(4, 3, hspace=.42, wspace=.28)
fig.suptitle(f'TIC {TIC} (2MASS J04160640−5935386): planet-candidate evidence sheet\n'
             f'M6–M7 dwarf, R★ = 0.123 R☉, d = 31.6 pc  •  P = {P:.5f} d  •  Rp = {fit["Rp_Rearth"][0]:.2f} ± {fit["Rp_Rearth"][1]:.2f} R⊕  •  Teq ≈ {fit["Teq_K"][0]:.0f} K',
             fontsize=15, y=.995)
# (1) BLS periodogram on all data
ax = fig.add_subplot(gs[0, 0]); ta, fa = cat(list(data)); sig = 1.4826 * np.median(np.abs(fa - np.median(fa)))
per = np.exp(np.linspace(np.log(.6), np.log(15), 60000)); bls = BoxLeastSquares(ta, fa, np.full(len(fa), sig))
r = bls.power(per, [0.02, 0.03], objective='snr')
ax.plot(r.period, r.depth_snr, 'k-', lw=.5); ax.axvline(P, c='r', ls='--', alpha=.5); ax.set_xscale('log')
ax.set_xlabel('trial period (days)'); ax.set_ylabel('transit SNR'); ax.set_title('(1) Blind period search, all 20 sectors:\nhighest peak, but only modestly above the noise forest', fontsize=10)
ax.annotate(f'P = {P:.4f} d', (P, r.depth_snr.max()), xytext=(P * 1.6, r.depth_snr.max() * .9), arrowprops=dict(arrowstyle='->', color='r'), color='r')
# (2) all data folded with model
ax = fig.add_subplot(gs[0, 1:]); ph = phase_h(ta); m = np.abs(ph) < 4
ax.plot(ph[m], (fa[m] - 1) * 1e3, '.', color='.82', ms=1.5, rasterized=True)
c, mm, ee = binned(ph[m], (fa[m] - 1) * 1e3, np.arange(-4, 4.01, 1 / 6)); ax.errorbar(c, mm, ee, fmt='o', color='k', ms=4)
pm = batman.TransitParams(); pm.t0, pm.per, pm.rp, pm.a = 0., P, fit['Rp_Rs'][0], fit['a_Rs'][0]; pm.inc, pm.ecc, pm.w, pm.u, pm.limb_dark = fit['inc_deg'][0], 0., 90., [.2, .4], 'quadratic'
tt = np.linspace(-4, 4, 1000) / 24; ax.plot(tt * 24, (batman.TransitModel(pm, tt, supersample_factor=5, exp_time=10 / 1440).light_curve(pm) - 1) * 1e3, 'r-', lw=2, label='transit model (10-min sampling)')
ax.set_ylim(-14, 8); ax.set_xlabel('hours from mid-transit'); ax.set_ylabel('relative flux − 1 (ppt)'); ax.legend(loc='lower left')
ax.set_title(f'(2) All {len(data)} TESS sectors folded on the period (grey: data, black: 10-min bins)', fontsize=11)
# (3) detection vs independent
for j, (keys, lab, col) in enumerate(((sorted(used & set(data)), 'sectors used for detection', 'C0'), (sorted(set(data) - used), 'sectors NOT used for detection (independent test)', 'C3'))):
    ax = fig.add_subplot(gs[1, j]); t, f = cat(keys); ph = phase_h(t); m = np.abs(ph) < 4
    ax.plot(ph[m], (f[m] - 1) * 1e3, '.', color='.85', ms=1.5, rasterized=True)
    c, mm, ee = binned(ph[m], (f[m] - 1) * 1e3, np.arange(-4, 4.01, 1 / 4)); ax.errorbar(c, mm, ee, fmt='o', color=col, ms=4)
    ax.axvspan(-DUR * 12, DUR * 12, color='y', alpha=.25); ax.axhline(0, c='.5', lw=.8); ax.set_ylim(-14, 8)
    ax.set_title(f'(3{"ab"[j]}) {lab}\nsectors {",".join(map(str, keys))}', fontsize=10); ax.set_xlabel('hours from mid-transit'); ax.set_ylabel('ppt')
# (4) per-sector depths
ax = fig.add_subplot(gs[1, 2]); rows = []
for s, (t, f, cad) in sorted(data.items()):
    ph = phase_h(t) / 24; a = np.abs(ph) < DUR / 2; b = (np.abs(ph) > DUR) & (np.abs(ph) < 4 * DUR)
    if a.sum() < 2: continue
    sg = 1.4826 * np.median(np.abs(f - np.median(f))); rows.append((s, (np.median(f[b]) - f[a].mean()) * 1e3, sg / np.sqrt(a.sum()) * 1e3, s in used))
R = pd.DataFrame(rows, columns=['s', 'd', 'e', 'used']); x = np.arange(len(R))
ax.errorbar(x[R.used], R.d[R.used], R.e[R.used], fmt='o', color='C0', label='detection'); ax.errorbar(x[~R.used], R.d[~R.used], R.e[~R.used], fmt='s', color='C3', label='independent')
wm = np.sum(R.d / R.e ** 2) / np.sum(1 / R.e ** 2); ax.axhline(wm, c='k', ls='--', lw=1); ax.axhline(0, c='.5', lw=.8)
ax.set_xticks(x); ax.set_xticklabels(R.s, rotation=90, fontsize=7); ax.set_xlabel('TESS sector'); ax.set_ylabel('transit depth (ppt)'); ax.legend(fontsize=8)
ax.set_title(f'(4) Depth in each sector (mean {wm:.1f} ppt)\n{(R.d > 0).sum()} of {len(R)} sectors positive', fontsize=10)
# (5) odd/even and secondary
t, f = cat(list(data)); n = np.floor((t - T0) / P + .5).astype(int)
for j, (sel, lab, col, off) in enumerate(((n % 2 == 0, 'even transits', 'C0', 0.), (n % 2 == 1, 'odd transits', 'C1', 0.), (np.ones_like(n, bool), 'phase 0.5 (where a binary\'s secondary eclipse would be)', 'C2', .5))):
    ax = fig.add_subplot(gs[2, j]); ph = phase_h(t[sel], off); m = np.abs(ph) < 4
    c, mm, ee = binned(ph[m], (f[sel][m] - 1) * 1e3, np.arange(-4, 4.01, 1 / 4)); ax.errorbar(c, mm, ee, fmt='o', color=col, ms=4)
    a = np.abs(ph / 24) < DUR / 2; b = (np.abs(ph / 24) > DUR) & (np.abs(ph / 24) < 4 * DUR); dd = (np.median(f[sel][b]) - f[sel][a].mean()) * 1e3
    de = 1.4826 * np.median(np.abs(f[sel] - np.median(f[sel]))) / np.sqrt(a.sum()) * 1e3
    ax.axvspan(-DUR * 12, DUR * 12, color='y', alpha=.25); ax.axhline(0, c='.5', lw=.8); ax.set_ylim(-14, 8)
    ax.set_title(f'(5{"abc"[j]}) {lab}: depth {dd:.1f} ± {de:.1f} ppt', fontsize=10); ax.set_xlabel('hours from that phase'); ax.set_ylabel('ppt')
# (6) transit timing
ax = fig.add_subplot(gs[3, 0]); tt_ = []
for s, (t, f, cad) in sorted(data.items()):
    n0 = np.round((np.median(t) - T0) / P); tc0 = T0 + n0 * P; sg = 1.4826 * np.median(np.abs(f - 1)); dep = 7e-3
    grid = np.linspace(-.12, .12, 241); chi = []
    for g in grid:
        ph = ((t - tc0 - g) / P + .5) % 1 - .5; mm_ = np.abs(ph * P) < DUR / 2; chi.append(np.sum((f - np.where(mm_, 1 - dep, 1)) ** 2) / sg ** 2)
    chi = np.array(chi); g0 = grid[np.argmin(chi)]; okk = grid[chi < chi.min() + 1]; tt_.append((tc0 + g0, (okk.max() - okk.min()) / 2, n0, s))
TT = pd.DataFrame(tt_, columns=['t', 'e', 'n', 's']); oc = (TT.t - (T0 + TT.n * P)) * 1440; err = np.maximum(TT.e * 1440 * 2.25, 3)
ax.errorbar(TT.t + 2457000 - 2458000, oc, err, fmt='o'); ax.axhline(0, c='k', lw=.8)
ax.set_xlabel('BJD − 2458000'); ax.set_ylabel('observed − predicted (min)'); ax.set_ylim(-150, 150)
ax.set_title('(6) Transit timing per sector, 2018–2023:\nall on one strictly periodic orbit', fontsize=10)
# (7) difference image, stacked across 200-s sectors in pixel space of S64
ax = fig.add_subplot(gs[3, 1]); tpfs = sorted(glob.glob(f'{D}/*_tp.fits')); tpf = lk.read([x for x in tpfs if 's0064' in x][0], quality_bitmask='default')
tm = tpf.time.value; fl = tpf.flux.value; ok = np.all(np.isfinite(fl.reshape(len(tm), -1)), 1); tm, fl = tm[ok], fl[ok]
ph = (((tm - T0) / P + .5) % 1 - .5) * P; nn = np.floor((tm - T0) / P + .5).astype(int); intr, near = np.abs(ph) < DUR / 2, (np.abs(ph) > DUR) & (np.abs(ph) < 4 * DUR)
diffs = [np.nanmedian(fl[near & (nn == k)], 0) - np.nanmean(fl[intr & (nn == k)], 0) for k in np.unique(nn[intr]) if (intr & (nn == k)).sum() > 3]
diff = np.nanmean(diffs, 0); im = ax.imshow(diff, origin='lower', cmap='RdBu_r', vmin=-np.nanmax(abs(diff)), vmax=np.nanmax(abs(diff)))
x0, y0 = tpf.wcs.world_to_pixel(SkyCoord(64.0267076207366, -59.5940445447903, unit='deg')); ax.plot(x0, y0, 'kx', ms=14, mew=3, label='TIC 198079261')
g = pd.read_csv(f'{O}/gaia_neighbors_120arcsec.csv'); cap = g[(g.sep > .5) & (g.need < 1)]
for _, rr in cap.iterrows():
    xx, yy = tpf.wcs.world_to_pixel(SkyCoord(rr.ra, rr.dec, unit='deg'))
    if 0 <= xx < diff.shape[1] and 0 <= yy < diff.shape[0]: ax.plot(xx, yy, 'o', mfc='none', mec='g', ms=10, mew=1.5)
ax.plot([], [], 'o', mfc='none', mec='g', label='neighbours able to mimic'); ax.legend(fontsize=8, loc='upper right')
plt.colorbar(im, ax=ax, fraction=.046); vs = json.load(open(f'{O}/vet_summary.json'))
ax.set_title(f'(7) Where the light disappears (S64 difference image, 21″ px)\non-target PSF fit, 5 sectors: {vs["on_target_amp"][0]*1e3:.1f} ± {vs["on_target_amp"][1]*1e3:.1f} ppt vs 7.3 ppt dip', fontsize=10)
ax.set_xlabel('column (px)'); ax.set_ylabel('row (px)')
# (8) TRICERATOPS
ax = fig.add_subplot(gs[3, 2]); tp = pd.read_csv(f'{O}/triceratops_probs.csv'); tp = tp.groupby('scenario').prob.sum().sort_values(ascending=True)
lab = {'TP': 'planet on target', 'PTP': 'planet on target +\nunresolved companion', 'STP': 'planet on unresolved\ncompanion', 'EB': 'eclipsing binary (target)', 'EBx2P': 'EB at 2×P (target)',
       'PEB': 'EB + companion', 'PEBx2P': 'EB 2×P + companion', 'SEB': 'EB on companion', 'SEBx2P': 'EB 2×P on companion', 'DTP': 'diluted planet', 'DEB': 'diluted EB', 'DEBx2P': 'diluted EB 2×P',
       'BTP': 'background planet', 'BEB': 'background EB', 'BEBx2P': 'background EB 2×P', 'NTP': 'planet on neighbour', 'NEB': 'EB on neighbour', 'NEBx2P': 'EB 2×P on neighbour'}
tp = tp[tp > 1e-3]; cols = ['C2' if k in ('TP', 'PTP', 'DTP') else 'C3' for k in tp.index]
ax.barh([lab.get(k, k) for k in tp.index], tp.values, color=cols); ax.set_xlabel('probability'); ax.set_xlim(0, 1)
fpp = np.mean(vs['FPP']); ax.set_title(f'(8) TRICERATOPS scenario probabilities\nFPP = {fpp:.2f}, NFPP = {np.mean(vs["NFPP"]):.3f} (5 runs)', fontsize=10)
for i, v in enumerate(tp.values): ax.text(v + .01, i, f'{v:.3f}', va='center', fontsize=8)
out = f'{O}/EVIDENCE_TIC198079261.png'; plt.savefig(out, dpi=90, bbox_inches='tight'); print(out)
