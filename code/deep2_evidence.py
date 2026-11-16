"""One-page evidence figure for a deep-pass candidate from the vetting products in deep2/vet and deep2/pix.  Usage: deep2_evidence.py TIC K OUTDIR"""
import sys, os, json, glob, re, numpy as np, pandas as pd, warnings, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt, batman
warnings.filterwarnings('ignore')
R = '/home/sarang/tess_search'; TIC, K, OUT = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]; V = f'{R}/deep2/vet'; os.makedirs(OUT, exist_ok=True)
v = json.load(open(f'{V}/{TIC}_{K}.json')); fit = json.load(open(f'{V}/{TIC}_{K}.fit.json')); pix = json.load(open(f'{V}/{TIC}_{K}.pix.json')); loc = json.load(open(f'{V}/{TIC}_{K}.loc.json')) if os.path.exists(f'{V}/{TIC}_{K}.loc.json') else None
fr = json.load(open(f'{V}/{TIC}_{K}.fit_freerho.json')) if os.path.exists(f'{V}/{TIC}_{K}.fit_freerho.json') else None
tr = json.load(open(f'{V}/{TIC}_{K}.tric.json')) if os.path.exists(f'{V}/{TIC}_{K}.tric.json') else None
nl = open(f'{V}/null_{TIC}.log').read() if os.path.exists(f'{V}/null_{TIC}.log') else ''; mnl = re.search(r'scrambled max SNR ([\d.]+), max SDE ([\d.]+); inverted ([\d.]+); real ([\d.]+) \(SDE ([\d.]+)\)', nl)
z = np.load(f'{V}/{TIC}.prep.npz'); t, f, e, s = z['t'], z['f'], z['e'], z['s']; P, T0 = fit['period_d'][0], fit['T0_BTJD'][0]; T14 = fit['T14_h'][0] / 24
ph = (((t - T0) / P + .5) % 1 - .5) * P * 24; n = np.floor((t - T0) / P + .5).astype(int); INK, MUTED, C1, C2, ACC = '#1f2430', '#6b7280', '#2a6fdb', '#d9822b', '#c0392b'
plt.rcParams.update({'font.size': 9, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK, 'xtick.color': MUTED, 'ytick.color': MUTED, 'axes.spines.top': False, 'axes.spines.right': False})
fig, ax = plt.subplots(2, 3, figsize=(15, 8.2)); ax = ax.ravel()
def binned(x, y, w, ed):
    idx = np.digitize(x, ed); bx = (ed[:-1] + ed[1:]) / 2
    by = np.array([np.average(y[idx == i], weights=w[idx == i]) if (idx == i).sum() > 2 else np.nan for i in range(1, len(ed))]); be = np.array([1 / np.sqrt(w[idx == i].sum()) if (idx == i).sum() > 2 else np.nan for i in range(1, len(ed))]); return bx, by, be
w = 1 / e ** 2; bx, by, be = binned(ph, f, w, np.linspace(-P * 12, P * 12, 81)); ax[0].plot(ph, (f - 1) * 1e6, '.', color='.88', ms=1, rasterized=True); ax[0].errorbar(bx, (by - 1) * 1e6, be * 1e6, fmt='o', color=INK, ms=3, lw=1)
dep = fit['depth_ppm'][0]; ax[0].set_ylim(-2.2 * dep, 1.4 * dep); ax[0].set_title('Whole orbit, all sectors folded', color=INK, loc='left'); ax[0].set_xlabel('hours from mid-transit'); ax[0].set_ylabel('ppm')
lim = 4 * T14 * 24; m = np.abs(ph) < lim; bx, by, be = binned(ph[m], f[m], w[m], np.linspace(-lim, lim, 37)); ax[1].plot(ph[m], (f[m] - 1) * 1e6, '.', color='.88', ms=1.5, rasterized=True); ax[1].errorbar(bx, (by - 1) * 1e6, be * 1e6, fmt='o', color=INK, ms=4, lw=1)
pr = batman.TransitParams(); pr.t0 = 0; pr.per = P; pr.rp = fit['Rp_Rs'][0]; pr.a = fit['a_Rs'][0]; pr.inc = fit['inc_deg'][0]; pr.ecc = 0; pr.w = 90; pr.u = [0.25, 0.35]; pr.limb_dark = 'quadratic'
xx = np.linspace(-lim, lim, 500) / 24; ax[1].plot(xx * 24, (batman.TransitModel(pr, xx).light_curve(pr) - 1) * 1e6, color=ACC, lw=2); ax[1].set_ylim(-2.2 * dep, 1.4 * dep)
ax[1].set_title(f'Transit ({v["ntr"]} transits folded)', color=INK, loc='left'); ax[1].set_xlabel('hours from mid-transit'); ax[1].set_ylabel('ppm')
def dp(mm):
    a = mm & (np.abs(ph) < T14 * 12); b = mm & (np.abs(ph) > T14 * 18) & (np.abs(ph) < T14 * 72)
    return ((f[b].mean() - f[a].mean()) * 1e6, np.sqrt(np.mean(e[a] ** 2) / a.sum() + np.mean(e[b] ** 2) / b.sum()) * 1e6) if a.sum() > 2 and b.sum() > 5 else (np.nan, np.nan)
ss = sorted(np.unique(s)); d = np.array([dp(s == u) for u in ss]); ax[2].errorbar(range(len(ss)), d[:, 0], d[:, 1], fmt='o', color=INK, ms=5, lw=1); ax[2].axhline(0, color=MUTED, lw=.6); wm = np.nansum(d[:, 0] / d[:, 1] ** 2) / np.nansum(1 / d[:, 1] ** 2); ax[2].axhline(wm, color=ACC, lw=1, ls='--')
st = max(1, len(ss) // 12); ax[2].set_xticks(range(0, len(ss), st)); ax[2].set_xticklabels([f'S{int(u)}' for u in ss[::st]], rotation=90 if len(ss) > 8 else 0); ax[2].set_ylim(min(0, np.nanmin(d[:, 0] - d[:, 1])) - 0.15 * dep, np.nanmax(d[:, 0] + d[:, 1]) + 0.15 * dep); ax[2].set_title(f'Depth per sector: {int((d[:,0]>0).sum())} of {len(ss)} positive', color=INK, loc='left'); ax[2].set_ylabel('ppm')
for lab, mk, c in (('odd transits', n % 2 == 1, C1), ('even transits', n % 2 == 0, C2)):
    mm = mk & (np.abs(ph) < lim); bx, by, be = binned(ph[mm], f[mm], w[mm], np.linspace(-lim, lim, 25)); ax[3].errorbar(bx, (by - 1) * 1e6, be * 1e6, fmt='o-', color=c, ms=3, lw=1.2, label=lab)
ax[3].legend(frameon=False, labelcolor=INK); ax[3].set_ylim(-2.2 * dep, 1.4 * dep); ax[3].set_xlabel('hours from mid-transit'); ax[3].set_ylabel('ppm')
ax[3].set_title(f'Odd vs even transits ({v["oe_sig"]:.1f}σ apart)', color=INK, loc='left')
lb = f'{R}/deep2/pix/{TIC}/locbest_{K}.npz'
if os.path.exists(lb):
    z = np.load(lb); ax[4].imshow(z['diff'], origin='lower', cmap='Greys'); ax[4].plot(z['x0'], z['y0'], '+', color=ACC, ms=18, mew=2)
    ok = (z['pid'] != TIC) & (z['ptm'] < v['Tmag'] + 5); ax[4].plot(z['px'][ok], z['py'][ok], 'o', mfc='none', mec=C1, ms=9, mew=1.2)
    ax[4].set_xlim(-.5, z['diff'].shape[1] - .5); ax[4].set_ylim(-.5, z['diff'].shape[0] - .5); ax[4].set_title(f'Where the light drops (sector {int(z["sec"])})', color=INK, loc='left'); ax[4].set_xlabel('pixel (21″ each);  + the star,  o other stars'); ax[4].set_ylabel('pixel')
fs = [] if os.path.exists(lb) else sorted(glob.glob(f'{R}/deep2/pix/{TIC}/diff_s*.npy')); best = None
for fn in fs:
    a = np.load(fn); sn = np.nanmax(a[0]) / (np.nanstd(a[0]) + 1e-9)
    if best is None or sn > best[0]: best = (sn, a, int(re.search(r'diff_s(\d+)', fn).group(1)))
if best:
    a = best[1]; im = ax[4].imshow(a[0], origin='lower', cmap='Greys'); yy, xx = np.indices(a[1].shape); cen = np.hypot(xx - (a[1].shape[1] - 1) / 2, yy - (a[1].shape[0] - 1) / 2) < 2.0; cy, cx = np.unravel_index(np.nanargmax(np.where(cen, a[1], -np.inf)), a[1].shape); ax[4].plot(cx, cy, '+', color=ACC, ms=16, mew=2)
    ax[4].set_title(f'Where the light drops (sector {best[2]}); + = the star', color=INK, loc='left'); ax[4].set_xlabel('pixel (21″ each)'); ax[4].set_ylabel('pixel')
ax[5].axis('off'); L = [f'TIC {TIC}   (TESS mag {v["Tmag"]:.2f}, R★ = {v["Rs"]:.2f} R☉)', '',
     f'Period        {fit["period_d"][0]:.6f} ± {max(fit["period_d"][1:]):.6f} d', f'Mid-transit   BJD {fit["T0_BJD"][0]:.5f} ± {max(fit["T0_BJD"][1:]):.5f}',
     f'Depth         {fit["depth_ppm"][0]:.0f} (+{fit["depth_ppm"][1]:.0f}/−{fit["depth_ppm"][2]:.0f}) ppm', f'Duration      {fit["T14_h"][0]:.2f} ± {max(fit["T14_h"][1:]):.2f} h',
     f'Radius        {fit["Rp_Rearth"][0]:.2f} ± {max(fit["Rp_Rearth"][1:]):.2f} R⊕   (Rp/R★ {fit["Rp_Rs"][0]:.4f}, b {fit["b"][0]:.2f})', f'Equilibrium T {fit["Teq_K"][0]:.0f} K', '',
     f'Secondary eclipse     {v["sec_ppm"][0]:.0f} ± {v["sec_ppm"][1]:.0f} ppm', (f'Depth in pixels       {loc["depth_on_target_ppm"]:.0f} ± {loc["depth_on_target_err"]:.0f} ppm ({loc["det_sigma"]:.1f}σ)' if loc else f'On-target amplitude   {pix["amp_ppm"]:.0f} ± {pix["amp_err"]:.0f} ppm ({pix["emp_sigma"]:.1f}σ vs other phases)'),
     (f'Dip position          on the star within 4″ ({loc["nsec"]} sectors)\n' + (f'Nearest other star    disfavoured at {loc["alt_nsig"]:.1f}σ' if loc.get("alt_tic") else 'Nearest other star    none can mimic it')) if loc else f'Dip offset from star  {pix["off_arcsec"]:.1f}″ ({pix["off_sigma"]:.1f}σ)']
if fr: L.append(f'Stellar density       {fr["rho_star_gcc"][0]:.1f} (+{fr["rho_star_gcc"][1]:.1f}/−{fr["rho_star_gcc"][2]:.1f}) g/cc (catalogue {fit["star"]["rho_prior"][0]:.1f})')
if mnl: L.append(f'Search SNR            {float(mnl.group(4)):.1f} (scrambled max {float(mnl.group(1)):.1f}, inverted {float(mnl.group(3)):.1f})')
if tr: L.append(f'TRICERATOPS           FPP {np.mean(tr["FPP"]):.3f} ± {np.std(tr["FPP"]):.3f}, NFPP {np.mean(tr["NFPP"]):.4f}')
ax[5].text(0, 1, '\n'.join(L), va='top', ha='left', family='monospace', fontsize=9, color=INK, transform=ax[5].transAxes)
fig.suptitle(f'TIC {TIC}: transit-like signal every {P:.4f} d  (planet candidate, not confirmed)', color=INK, fontsize=13, x=0.02, ha='left'); plt.tight_layout(rect=(0, 0, 1, 0.96))
plt.savefig(f'{OUT}/EVIDENCE_TIC{TIC}.png', dpi=110); print('wrote', f'{OUT}/EVIDENCE_TIC{TIC}.png')
