"""Streaming fetch + preprocess of TESS-SPOC FFI light curves (faint nearby M dwarfs). Resumable; writes out3/{tic}.npz and {tic}.meta.json.
Usage: fetch_proc.py targets.csv  (env NPROC, MAXSEC)"""
import sys, os, json, shutil, urllib.request, warnings, multiprocessing as mp, numpy as np, pandas as pd
warnings.filterwarnings('ignore')
D = '/home/sarang/tess_search'; OUT = f'{D}/out3'; os.makedirs(OUT, exist_ok=True)
S3 = 'https://archive.stsci.edu/hlsps/tess-spoc'
AV = None; MAXSEC = int(os.environ.get('MAXSEC', 8))


def init():
    global AV
    AV = pd.read_pickle(f'{D}/data/avail_spoc_faint.pkl')


def one(tic):
    meta = f'{OUT}/{tic}.meta.json'
    if os.path.exists(meta):
        if not os.environ.get('RETRY_ERR'): return None
        try:
            if json.load(open(meta)).get('status') == 'ok': return None
        except Exception: pass
    import lightkurve as lk
    secs = sorted(AV[tic]); t16 = f'{int(tic):016d}'
    secs = [(x, f'hlsp_tess-spoc_tess_phot_{t16}-s{x:04d}_tess_v1_lc.fits') for x in secs]
    nums = [x for x, _ in secs]
    best = max(((i, [x for x in secs if nums[i] <= x[0] < nums[i] + 14][:MAXSEC]) for i in range(len(secs))), key=lambda t: len(t[1]))[1]
    d = f'/tmp/ff_{tic}'; os.makedirs(d, exist_ok=True); parts = []; res = dict(tic=int(tic), sectors=[s for s, _ in best], all_sectors=nums)
    try:
        for s, fn in best:
            p = f'{d}/{fn}'
            try:
                for _k in range(4):
                    try:
                        urllib.request.urlretrieve(f'{S3}/s{s:04d}/target/{t16[0:4]}/{t16[4:8]}/{t16[8:12]}/{t16[12:16]}/{fn}', p); break
                    except Exception:
                        if _k == 3: raise
                lc = lk.read(p, quality_bitmask='default', flux_column='pdcsap_flux')
            except Exception:
                continue
            finally:
                if os.path.exists(p): os.remove(p)
            lc = lc.remove_nans().normalize(); lc = lc[np.isfinite(lc.flux.value)]
            if len(lc) < 500: continue
            cad = np.nanmedian(np.diff(lc.time.value))
            # flares: >4 sigma above a running median on >=2 consecutive cadences; cut flare + 6 h decay before detrending
            from scipy.ndimage import median_filter
            fv = np.asarray(lc.flux.value, float); tv = np.asarray(lc.time.value, float)
            rm = median_filter(fv, size=int(round(0.25 / cad)) | 1, mode='nearest'); rr = fv - rm
            sg = 1.4826 * np.median(np.abs(rr)); hi = rr > 4 * sg
            fl_idx = np.flatnonzero(hi[:-1] & hi[1:])
            cut = np.zeros(len(tv), bool); ftimes = []
            for i in fl_idx:
                if ftimes and tv[i] - ftimes[-1] < 0.05: continue
                ftimes.append(float(tv[i])); cut |= (tv > tv[i] - 0.02) & (tv < tv[i] + 0.25)
            res.setdefault('flares', []).extend(ftimes)
            lc = lc[~cut]
            if len(lc) < 500: continue
            fl = lc.flatten(window_length=int(round(0.75 / cad)) | 1).remove_outliers(sigma_upper=3, sigma_lower=20)
            parts.append(fl.bin(time_bin_size=10 / 1440).remove_nans())
        if not parts: raise RuntimeError('no usable sectors')
        t = np.concatenate([x.time.value for x in parts]); f = np.concatenate([x.flux.value for x in parts])
        ok = np.isfinite(t) & np.isfinite(f); t, f = t[ok], f[ok]
        np.savez_compressed(f'{OUT}/{tic}.npz', t=t, f=f)
        res.update(status='ok', npts=int(len(t)), used=[s for s, _ in best])
    except Exception as e:
        res['status'] = f'err: {e}'[:200]
    finally:
        shutil.rmtree(d, ignore_errors=True)
    json.dump(res, open(meta, 'w'))
    return tic


if __name__ == '__main__':
    tg = pd.read_csv(sys.argv[1])
    with mp.Pool(int(os.environ.get('NPROC', 8)), initializer=init) as pool:
        for i, _ in enumerate(pool.imap(one, tg.ID.astype(int).tolist(), chunksize=1)):
            if i % 200 == 0: print(i, flush=True)
    print('done')
