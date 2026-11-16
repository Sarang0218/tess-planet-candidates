"""Where on the sky is the dip?  Gaia DR3 neighbours that could mimic the depth, then per-transit difference images (out - in)
from TESS-SPOC FFI target-pixel files: on-target amplitude, difference-image centroid vs out-of-transit centroid (bootstrap over
transits), and the same statistic at 40 other orbital phases as the null.   Usage: deep2_pix.py TIC K [max_sectors]"""
import sys, os, io, json, glob, numpy as np, pandas as pd, urllib.request, urllib.parse, warnings, lightkurve as lk
from astropy.coordinates import SkyCoord
from scipy.ndimage import median_filter
warnings.filterwarnings('ignore')
import socket; socket.setdefaulttimeout(180)                  # a stalled archive connection must not hang the run
R = '/home/sarang/tess_search'; TIC, K = int(sys.argv[1]), int(sys.argv[2]); NMAX = int(sys.argv[3]) if len(sys.argv) > 3 else 5
v = json.load(open(f'{R}/deep2/vet/{TIC}_{K}.json')); P, t0, dur, depth = v['P'], v['T0'], v['T14_h'] / 24, v['trap_depth_ppm'] * 1e-6
tg = pd.concat([pd.read_csv(f'{R}/data/{n}.csv') for n in ('targets', 'targets_allsky', 'targets_ffi_faint')]).drop_duplicates('ID').set_index('ID').loc[TIC]; ra0, dec0 = float(tg.ra), float(tg.dec)
D = f'{R}/deep2/pix/{TIC}'; os.makedirs(D, exist_ok=True); t16 = f'{TIC:016d}'; out = dict(tic=TIC, k=K)
from astroquery.mast import Catalogs                      # neighbours from the TESS Input Catalog (Gaia DR2 based, TESS-band magnitudes)
g = None
for _ in range(3):
    try: g = Catalogs.query_region(f'{ra0} {dec0}', radius=120 / 3600, catalog='TIC').to_pandas(); break
    except Exception as ex: print('tic retry', ex)
if g is not None and len(g):
    g['ID'] = g.ID.astype(int); me = g[g.ID == TIC].iloc[0]; g['sep'] = g.dstArcSec.astype(float); g['dT'] = g.Tmag - me.Tmag
    g['need'] = depth * (1 + 10 ** (-.4 * g.dT)) / 10 ** (-.4 * g.dT)                 # eclipse depth the neighbour would need to produce the signal
    cap = g[(g.ID != TIC) & (g.need < 1)].sort_values('sep')
    print(f'TIC {TIC}: Tmag {me.Tmag:.2f} plx={me.plx:.2f} Gaia DR2 {me.GAIA}; neighbours within 120" able to mimic {depth*1e6:.0f} ppm: {len(cap)} (within 42": {(cap.sep < 42).sum()}, within 21": {(cap.sep < 21).sum()})')
    print(cap[cap.sep < 63][['ID', 'Tmag', 'dT', 'sep', 'need', 'plx']].round(3).head(12).to_string(index=False)); g.to_csv(f'{D}/tic_120arcsec.csv', index=False)
    out.update(gaia_dr2=str(me.GAIA), n_mimic_42=int((cap.sep < 42).sum()), n_mimic_21=int((cap.sep < 21).sum()), nearest_mimic=cap[['ID', 'sep', 'dT', 'need']].head(6).values.tolist())
_pl = pd.read_csv(f'{R}/deep2/targets.csv', usecols=['ID', 'plan']).set_index('ID'); PLAN = json.loads(_pl.loc[TIC, 'plan']) if TIC in _pl.index else {}
src = json.load(open(f'{R}/deep2/lc/{TIC}.meta.json'))['src']; secs = sorted([int(s) for s, k in src.items()], reverse=True)


def stat(t, fl, shift, x0, y0, rng=None):
    ph = (((t - t0) / P + .5 + shift) % 1 - .5) * P; nn = np.floor((t - t0) / P + .5 + shift).astype(int)
    intr, near = np.abs(ph) < 0.4 * dur, (np.abs(ph) > 0.7 * dur) & (np.abs(ph) < 3.5 * dur)
    diffs = [np.nanmean(fl[near & (nn == k)], 0) - np.nanmean(fl[intr & (nn == k)], 0) for k in np.unique(nn[intr]) if (intr & (nn == k)).sum() >= 1 and (near & (nn == k)).sum() >= 6]
    if len(diffs) < 3: return None
    diffs = np.array(diffs); oot = np.nanmedian(fl[near], 0); yy, xx = np.indices(oot.shape); m = np.hypot(xx - x0, yy - y0) < 2.5
    def one(dd):
        diff = dd.mean(0); err = dd.std(0) / np.sqrt(len(dd)) + 1e-9; w = 1 / err[m] ** 2
        a = np.sum(w * diff[m] * oot[m]) / np.sum(w * oot[m] ** 2); ae = np.sqrt(1 / np.sum(w * oot[m] ** 2)); c = np.clip(diff[m], 0, None)
        return a, ae, np.sum(xx[m] * c) / c.sum(), np.sum(yy[m] * c) / c.sum()
    a, ae, cx, cy = one(diffs); o = np.clip(oot[m], 0, None); ox, oy = np.sum(xx[m] * o) / o.sum(), np.sum(yy[m] * o) / o.sum()
    bs = np.array([one(diffs[rng.integers(0, len(diffs), len(diffs))])[2:] for _ in range(200)]) if rng is not None else np.zeros((2, 2))
    return a, ae, len(diffs), cx - ox, cy - oy, bs[:, 0].std(), bs[:, 1].std(), diffs.mean(0), oot


res, null, rng, got = [], {}, np.random.default_rng(0), 0
for s in secs:
    if got >= NMAX: break
    n = f'hlsp_tess-spoc_tess_phot_{t16}-s{s:04d}_tess_v1_tp.fits'; url = f'https://archive.stsci.edu/hlsps/tess-spoc/s{s:04d}/target/{t16[0:4]}/{t16[4:8]}/{t16[8:12]}/{t16[12:16]}/{n}'
    if str(PLAN.get(str(s), '')).startswith('tm:'):                     # 2-min target: SPOC target-pixel file on the public bucket
        n = PLAN[str(s)][3:].replace('_lc.fits', '_tp.fits'); x = n.split('-')[2]; url = f'https://stpubdata.s3.amazonaws.com/tess/public/tid/s{s:04d}/{x[0:4]}/{x[4:8]}/{x[8:12]}/{x[12:16]}/{n}'
    p = f'{D}/{n}'
    for _try in range(3):
        if os.path.exists(p) and os.path.getsize(p) > 100000: break
        try: urllib.request.urlretrieve(url, p + '.part'); os.replace(p + '.part', p)
        except Exception as ex: print('  tpf', s, ex)
    if not (os.path.exists(p) and os.path.getsize(p) > 100000): continue
    try: tpf = lk.read(p, quality_bitmask='default')
    except Exception as ex: print('  read', s, ex); continue
    t = tpf.time.value; fl = np.asarray(tpf.flux.value, float); ok = np.all(np.isfinite(fl.reshape(len(t), -1)), 1); t, fl = t[ok], fl[ok]
    cad = np.nanmedian(np.diff(t)); kk = int(round(0.75 / cad)) | 1; fl = fl - median_filter(fl, size=(kk, 1, 1), mode='nearest') + np.nanmedian(fl, 0)
    x0, y0 = tpf.wcs.world_to_pixel(SkyCoord(ra0, dec0, unit='deg')); r = stat(t, fl, 0.0, x0, y0, rng)
    if r is None: continue
    a, ae, ntr, dx, dy, ex_, ey_, diff, oot = r; got += 1
    sn = diff / (np.nanstd(diff[np.hypot(*np.indices(diff.shape)[::-1] - np.array([x0, y0])[:, None, None]) > 3.5]) + 1e-9); py, px = np.unravel_index(np.nanargmax(sn), sn.shape)
    print(f'  S{s}: {ntr} transits; on-target amplitude {a*1e6:.0f} ± {ae*1e6:.0f} ppm; difference-image centroid offset dx={dx*21:+.1f}±{ex_*21:.1f}" dy={dy*21:+.1f}±{ey_*21:.1f}"; brightest diff pixel ({px},{py}) vs target ({x0:.1f},{y0:.1f}), S/N {sn[py, px]:.1f}', flush=True)
    res.append((s, a, ae, ntr, dx * 21, dy * 21, ex_ * 21, ey_ * 21)); np.save(f'{D}/diff_s{s}.npy', np.array([diff, oot]))
    for sh in np.linspace(0.1, 0.9, 41):
        qn = stat(t, fl, sh, x0, y0)
        if qn: null.setdefault(round(sh, 3), []).append(qn[:2])
if res:
    A = np.array(res); w = 1 / A[:, 2] ** 2; wm = np.sum(A[:, 1] * w) / w.sum(); we = 1 / np.sqrt(w.sum())
    nl = np.array([np.sum(np.array(x)[:, 0] / np.array(x)[:, 1] ** 2) / np.sum(1 / np.array(x)[:, 1] ** 2) for x in null.values() if len(x) == len(res)])
    wx, wy = 1 / A[:, 6] ** 2, 1 / A[:, 7] ** 2; mx, my = np.sum(A[:, 4] * wx) / wx.sum(), np.sum(A[:, 5] * wy) / wy.sum(); sx, sy = 1 / np.sqrt(wx.sum()), 1 / np.sqrt(wy.sum())
    emp = (wm - nl.mean()) / nl.std() if len(nl) > 5 else np.nan
    print(f'  combined on-target amplitude {wm*1e6:.0f} ± {we*1e6:.0f} ppm (light-curve depth {depth*1e6:.0f}); null at {len(nl)} other phases: mean {nl.mean()*1e6:.0f}, std {nl.std()*1e6:.0f} -> {emp:.1f} sigma empirical')
    print(f'  mean centroid offset of the dip from the star: dx={mx:+.1f}±{sx:.1f}" dy={my:+.1f}±{sy:.1f}"  ({np.hypot(mx, my):.1f}", {np.hypot(mx / sx, my / sy):.1f} sigma)')
    out.update(amp_ppm=wm * 1e6, amp_err=we * 1e6, null_std=float(nl.std() * 1e6) if len(nl) else None, emp_sigma=float(emp), off_arcsec=float(np.hypot(mx, my)), off_sigma=float(np.hypot(mx / sx, my / sy)), sectors=[int(x) for x in A[:, 0]])
json.dump(out, open(f'{R}/deep2/vet/{TIC}_{K}.pix.json', 'w'), default=float)
