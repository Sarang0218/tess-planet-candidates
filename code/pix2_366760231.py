"""Pixel-level location of the 1.1434-d signal on TIC 366760231: per-transit difference images (out - in) from 2-min and FFI
target-pixel files, fitted with the target's own out-of-transit image. The identical fit at 40 other phases gives the null."""
import os, glob, json, numpy as np, urllib.request, warnings, lightkurve as lk
from astropy.coordinates import SkyCoord
warnings.filterwarnings('ignore')
TIC = 366760231; t16 = f'{TIC:016d}'; D = f'/home/sarang/tess_search/deep/{TIC}'; O = '/home/sarang/tess_search/cand_366760231/two_planet'
eph = json.load(open(f'{O}/ephemeris2.json')); fit = json.load(open(f'{O}/fit2_results.json'))
P, t0, dur = eph['c']['P'], eph['c']['T0_BTJD'] + fit['c']['dt_min'][0] / 1440, fit['c']['T14_h'][0] / 24; Pb, Tb = eph['b']['P'], eph['b']['T0_BTJD']
RA, DEC = 335.11793, 74.46294
tm = {17: 'tess2019279210107-s0017-0000000366760231-0161-s', 18: 'tess2019306063752-s0018-0000000366760231-0162-s', 52: 'tess2022138205153-s0052-0000000366760231-0224-s'}
for s, n in tm.items():
    p = f'{D}/{n}_tp.fits'
    if not os.path.exists(p):
        try: urllib.request.urlretrieve(f'https://stpubdata.s3.amazonaws.com/tess/public/tid/s{s:04d}/{t16[0:4]}/{t16[4:8]}/{t16[8:12]}/{t16[12:16]}/{n}_tp.fits', p)
        except Exception as ex: print('fail', s, ex)
for s in (58, 59, 85):
    n = f'hlsp_tess-spoc_tess_phot_{t16}-s{s:04d}_tess_v1_tp.fits'; p = f'{D}/{n}'
    if not os.path.exists(p):
        try: urllib.request.urlretrieve(f'https://archive.stsci.edu/hlsps/tess-spoc/s{s:04d}/target/{t16[0:4]}/{t16[4:8]}/{t16[8:12]}/{t16[12:16]}/{n}', p)
        except Exception as ex: print('fail', s, ex)
def amp_at(t, fl, x0, y0, shift):
    ph = (((t - t0) / P + .5 + shift) % 1 - .5) * P; nn = np.floor((t - t0) / P + .5 + shift).astype(int)
    intr, near = np.abs(ph) < 0.45 * dur, (np.abs(ph) > 0.7 * dur) & (np.abs(ph) < 3.5 * dur)
    diffs = [np.nanmean(fl[near & (nn == k)], 0) - np.nanmean(fl[intr & (nn == k)], 0) for k in np.unique(nn[intr]) if (intr & (nn == k)).sum() >= 3 and (near & (nn == k)).sum() >= 10]
    if len(diffs) < 5: return None
    diffs = np.array(diffs); diff = diffs.mean(0); err = diffs.std(0) / np.sqrt(len(diffs)); oot = np.nanmedian(fl[near], 0)
    yy, xx = np.indices(diff.shape); m = np.hypot(xx - x0, yy - y0) < 2.5; w = 1 / err[m] ** 2
    a = np.sum(w * diff[m] * oot[m]) / np.sum(w * oot[m] ** 2); ae = np.sqrt(1 / np.sum(w * oot[m] ** 2))
    sn = diff / err; cy = np.sum(yy[m] * np.clip(diff[m], 0, None)) / np.sum(np.clip(diff[m], 0, None)); cx = np.sum(xx[m] * np.clip(diff[m], 0, None)) / np.sum(np.clip(diff[m], 0, None))
    return a, ae, len(diffs), cx, cy
res = []; null = {}
for p in sorted(glob.glob(f'{D}/*_tp.fits')):
    tpf = lk.read(p, quality_bitmask='default'); t = tpf.time.value; fl = tpf.flux.value; ok = np.all(np.isfinite(fl.reshape(len(t), -1)), 1); t, fl = t[ok], fl[ok]
    keep = ~(np.abs(((t - Tb) / Pb + .5) % 1 - .5) * Pb < 0.03); t, fl = t[keep], fl[keep]
    # remove slow per-pixel trends (0.75-d running median) so that differences are not dominated by drifts
    from scipy.ndimage import median_filter
    cad = np.nanmedian(np.diff(t)); k = int(round(0.75 / cad)) | 1; base = median_filter(fl, size=(k, 1, 1), mode='nearest'); fl = fl - base + np.nanmedian(fl, 0)
    x0, y0 = tpf.wcs.world_to_pixel(SkyCoord(RA, DEC, unit='deg')); r = amp_at(t, fl, x0, y0, 0.0)
    if r is None: continue
    a, ae, n, cx, cy = r; print(f'S{tpf.sector}: {n} transits; on-target amplitude {a*1e6:.0f} ± {ae*1e6:.0f} ppm; difference centroid ({cx:.2f},{cy:.2f}) vs target ({x0:.2f},{y0:.2f})')
    res.append((a, ae)); 
    for sh in np.linspace(0.1, 0.9, 41):
        q = amp_at(t, fl, x0, y0, sh)
        if q: null.setdefault(round(sh, 3), []).append(q[:2])
A = np.array(res); wm = np.sum(A[:, 0] / A[:, 1] ** 2) / np.sum(1 / A[:, 1] ** 2); we = 1 / np.sqrt(np.sum(1 / A[:, 1] ** 2))
nl = np.array([np.sum(np.array(v)[:, 0] / np.array(v)[:, 1] ** 2) / np.sum(1 / np.array(v)[:, 1] ** 2) for v in null.values()])
print(f'combined on-target amplitude: {wm*1e6:.0f} ± {we*1e6:.0f} ppm ({wm/we:.1f} sigma formal); light-curve depth {fit["c"]["depth_ppm"][0]:.0f} ppm')
print(f'null (same fit at {len(nl)} other phases): mean {nl.mean()*1e6:.0f}, std {nl.std()*1e6:.0f}, max {nl.max()*1e6:.0f} ppm -> {(wm - nl.mean()) / nl.std():.1f} sigma empirical')
