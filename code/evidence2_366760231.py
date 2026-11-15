"""Evidence sheet for the 1.1434-d signal on TIC 366760231 (three panels: folded transit + model, per-sector depths, noise-floor test)."""
import glob, json, re, numpy as np, warnings, lightkurve as lk, batman, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')
D = '/home/sarang/tess_search/deep/366760231'; O = '/home/sarang/tess_search/cand_366760231/two_planet'
E = json.load(open(f'{O}/ephemeris_final.json')); fit = json.load(open(f'{O}/fit2_results.json'))
P, T0 = E['c']['P'], E['c']['T0_BTJD']; Pb, Tb = E['b']['P'], E['b']['T0_BTJD']; d2 = fit['c']['T14_h'][0] / 24
INK, MUTE, BLUE, GRID, SURF = '#0b0b0b', '#52514e', '#2a78d6', '#e4e3df', '#fcfcfb'
plt.rcParams.update({'font.size': 10, 'axes.edgecolor': MUTE, 'axes.labelcolor': MUTE, 'xtick.color': MUTE, 'ytick.color': MUTE, 'axes.spines.top': False, 'axes.spines.right': False,
                     'figure.facecolor': SURF, 'axes.facecolor': SURF, 'axes.titlecolor': INK, 'axes.titlesize': 11, 'axes.titleweight': 'bold', 'axes.titlelocation': 'left'})
sec = {}
for fn in sorted(glob.glob(f'{D}/*s_lc.fits')) + sorted(glob.glob(f'{D}/hlsp_tess-spoc*_lc.fits')):
    lc = lk.read(fn, quality_bitmask='default', flux_column='pdcsap_flux').remove_nans().normalize(); lc = lc[np.isfinite(lc.flux.value)]; t = lc.time.value; cad = np.nanmedian(np.diff(t))
    m = (np.abs(((t - T0) / P + .5) % 1 - .5) * P < d2) | (np.abs(((t - Tb) / Pb + .5) % 1 - .5) * Pb < 0.03)
    fl = lc.flatten(window_length=int(round(.75 / cad)) | 1, mask=m).remove_outliers(sigma_upper=4, sigma_lower=50)
    t = np.asarray(fl.time.value, float); f = np.ma.filled(np.ma.asarray(fl.flux.value, dtype=float), np.nan); ok = np.isfinite(f) & ~(np.abs(((t - Tb) / Pb + .5) % 1 - .5) * Pb < 0.03)
    sec[lc.meta['SECTOR']] = (t[ok], f[ok], 1.4826 * np.median(np.abs(f[ok] - np.median(f[ok]))))
def depth(t, f, sg):
    ph = (((t - T0) / P + .5) % 1 - .5) * P; a = np.abs(ph) < 0.45 * d2; b = (np.abs(ph) > 0.7 * d2) & (np.abs(ph) < 3 * d2)
    return (f[b].mean() - f[a].mean()) * 1e6, sg * np.sqrt(1 / a.sum() + 1 / b.sum()) * 1e6
fig, ax = plt.subplots(1, 3, figsize=(15, 4.4), gridspec_kw=dict(width_ratios=[1.5, 1, 1]))
# (a) fold
ph = np.concatenate([(((t - T0) / P + .5) % 1 - .5) * P * 24 for t, f, s in sec.values()]); fa = np.concatenate([f for t, f, s in sec.values()]); wa = np.concatenate([np.full(len(t), 1 / s ** 2) for t, f, s in sec.values()])
e = np.arange(-3.5, 3.5001, 1 / 6); idx = np.digitize(ph, e); x = (e[:-1] + e[1:]) / 2
y = np.array([np.average(fa[idx == q], weights=wa[idx == q]) for q in range(1, len(e))]); ye = np.array([1 / np.sqrt(wa[idx == q].sum()) for q in range(1, len(e))])
ax[0].errorbar(x, (y - 1) * 1e6, ye * 1e6, fmt='o', color=BLUE, ms=4.5, lw=1, zorder=3)
pm = batman.TransitParams(); pm.ecc = 0.; pm.w = 90.; pm.u = [0.2, 0.4]; pm.limb_dark = 'quadratic'; rho = fit['rho_star_gcc'][0]; aR = (rho * 6.674e-8 * (P * 86400) ** 2 / (3 * np.pi)) ** (1 / 3)
pm.per, pm.t0, pm.rp, pm.a, pm.inc = P, 0., fit['c']['Rp_Rs'][0], aR, np.degrees(np.arccos(fit['c']['b'][0] / aR)); xx = np.linspace(-3.5, 3.5, 700)
ax[0].plot(xx, (batman.TransitModel(pm, xx / 24, supersample_factor=5, exp_time=10 / 1440).light_curve(pm) - 1) * 1e6, color=INK, lw=2, zorder=2)
ax[0].axhline(0, color=GRID, lw=1, zorder=1); ax[0].set_xlabel('hours from mid-transit'); ax[0].set_ylabel('relative flux (ppm)'); ax[0].set_ylim(-1150, 600)
ax[0].set_title(f'Folded at {P:.4f} d: a {fit["c"]["depth_ppm"][0]:.0f} ppm, {fit["c"]["T14_h"][0]:.1f} h dip'); ax[0].text(3.4, -1080, '10-min bins, 9 sectors, model in black', ha='right', color=MUTE, fontsize=9)
# (b) per sector
ss = sorted(sec); v = np.array([depth(*sec[s]) for s in ss]); w = 1 / v[:, 1] ** 2; dm = np.sum(v[:, 0] * w) / w.sum()
ax[1].axhline(0, color=GRID, lw=1); ax[1].axhline(dm, color=INK, lw=1, ls=(0, (4, 3))); ax[1].errorbar(range(len(ss)), v[:, 0], v[:, 1], fmt='o', color=BLUE, ms=6, lw=1.2)
ax[1].set_xticks(range(len(ss))); ax[1].set_xticklabels([str(s) for s in ss]); ax[1].set_xlabel('TESS sector (2019–2025)'); ax[1].set_ylabel('depth (ppm)')
ax[1].set_title('Present in all 9 sectors'); ax[1].set_ylim(-80, 1250); ax[1].text(-0.3, 1150, f'dashed: mean {dm:.0f} ppm', ha='left', color=MUTE, fontsize=9)
# (c) noise floor
txt = open('/home/sarang/tess_search/deep2/null_366760231.log').read(); sn = np.array([float(a) for a, b in re.findall(r'\(([\d.]+),([\d.]+)\)', txt.split('scrambled top')[1])])
ax[2].hist(sn, bins=np.arange(5.5, 7.3, 0.1), color='#9ec5f4', edgecolor=SURF, lw=1); ax[2].axvline(8.5, color=BLUE, lw=2)
ax[2].text(8.45, ax[2].get_ylim()[1] * 0.92, 'signal\nSNR 8.5', ha='right', va='top', color=INK, fontsize=9); ax[2].text(5.5, ax[2].get_ylim()[1] * 0.92, f'{len(sn)} scrambled\nsearches', va='top', color=MUTE, fontsize=9)
ax[2].set_xlim(5.4, 9); ax[2].set_xlabel('strongest signal found (SNR)'); ax[2].set_ylabel('searches'); ax[2].set_title('Above every scrambled search')
fig.suptitle('TIC 366760231: a second transit signal at 1.1434 d (candidate, not confirmed)', x=0.01, ha='left', fontsize=13, fontweight='bold', color=INK)
plt.tight_layout(rect=(0, 0, 1, 0.94)); plt.savefig(f'{O}/EVIDENCE_TIC366760231_02.png', dpi=130)
print('per-sector', np.round(v[:, 0]), 'mean', round(dm), 'null max', sn.max(), 'n', len(sn))
