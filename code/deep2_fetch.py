"""Deep pass 2 fetch + preprocess: ALL available sectors per star (SPOC 2-min where it is the better product, TESS-SPOC FFI otherwise).
Same cleaning as the first pass (flare cut, 0.75-d detrend, 10-min bins). Adds per-point sector id and per-sector noise.
Resumable; writes deep2/lc/{tic}.npz (t, f, e, s) and {tic}.meta.json.  Usage: deep2_fetch.py targets.csv  (env NPROC)"""
import sys, os, json, shutil, urllib.request, warnings, multiprocessing as mp, numpy as np, pandas as pd
warnings.filterwarnings('ignore')
D = '/home/sarang/tess_search'; OUT = f'{D}/deep2/lc'; TMP = f'{D}/deep2/tmp'; os.makedirs(OUT, exist_ok=True); os.makedirs(TMP, exist_ok=True)
FFI = 'https://archive.stsci.edu/hlsps/tess-spoc'; TM = 'https://stpubdata.s3.amazonaws.com/tess/public/tid'; NTHR = int(os.environ.get('NTHR', 4))


def get(url, p):
    for k in range(3):
        try:
            urllib.request.urlretrieve(url, p); return True
        except Exception as e:
            if 'HTTP Error 404' in str(e) or 'HTTP Error 403' in str(e): return False
    return False


def one(args):
    tic, plan = args
    meta = f'{OUT}/{tic}.meta.json'
    if os.path.exists(meta): return None
    import lightkurve as lk
    from scipy.ndimage import median_filter
    t16 = f'{int(tic):016d}'; d = f'{TMP}/{tic}'; os.makedirs(d, exist_ok=True)
    T, F, E, S = [], [], [], []; res = dict(tic=int(tic), planned=sorted(int(k) for k in plan), used=[], src={}, flares=[])
    try:
        items = []
        for s, src in sorted((int(k), v) for k, v in plan.items()):
            if src == 'ffi':
                fn = f'hlsp_tess-spoc_tess_phot_{t16}-s{s:04d}_tess_v1_lc.fits'; url = f'{FFI}/s{s:04d}/target/{t16[0:4]}/{t16[4:8]}/{t16[8:12]}/{t16[12:16]}/{fn}'
            else:
                fn = src[3:]; x = fn.split('-')[2]; url = f'{TM}/s{s:04d}/{x[0:4]}/{x[4:8]}/{x[8:12]}/{x[12:16]}/{fn}'
            items.append((s, src, url, f'{d}/{fn}'))
        from concurrent.futures import ThreadPoolExecutor
        with ThreadPoolExecutor(NTHR) as ex: got = list(ex.map(lambda it: get(it[2], it[3]), items))
        for (s, src, url, p), g_ in zip(items, got):
            try:
                if not g_: continue
                lc = lk.read(p, quality_bitmask='default', flux_column='pdcsap_flux')
                lc = lc.remove_nans().normalize(); lc = lc[np.isfinite(lc.flux.value)]
                if len(lc) < 300: continue
                cad = np.nanmedian(np.diff(lc.time.value))
                fv = np.asarray(lc.flux.value, float); tv = np.asarray(lc.time.value, float)
                rm = median_filter(fv, size=int(round(0.25 / cad)) | 1, mode='nearest'); rr = fv - rm
                sg = 1.4826 * np.median(np.abs(rr)); hi = rr > 4 * sg
                cut = np.zeros(len(tv), bool); last = -9
                for i in np.flatnonzero(hi[:-1] & hi[1:]):
                    if tv[i] - last < 0.05: continue
                    last = tv[i]; res['flares'].append(float(tv[i])); cut |= (tv > tv[i] - 0.02) & (tv < tv[i] + 0.25)
                lc = lc[~cut]
                if len(lc) < 300: continue
                fl = lc.flatten(window_length=int(round(0.75 / cad)) | 1).remove_outliers(sigma_upper=3, sigma_lower=20)
                if cad < 9.5 / 1440: fl = fl.bin(time_bin_size=10 / 1440).remove_nans()
                t = np.asarray(fl.time.value, float); f = np.ma.filled(np.ma.asarray(fl.flux.value, dtype=float), np.nan)
                ok = np.isfinite(t) & np.isfinite(f); t, f = t[ok], f[ok]
                if len(t) < 200: continue
                e = 1.4826 * np.median(np.abs(f - np.median(f)))
                if not np.isfinite(e) or e <= 0: continue
                T.append(t); F.append(f); E.append(np.full(len(t), e)); S.append(np.full(len(t), s, np.int16))
                res['used'].append(s); res['src'][str(s)] = 'ffi' if src == 'ffi' else '2min'
            except Exception:
                continue
            finally:
                if os.path.exists(p): os.remove(p)
        if not T: raise RuntimeError('no usable sectors')
        t = np.concatenate(T); o = np.argsort(t)
        np.savez_compressed(f'{OUT}/{tic}.npz', t=t[o], f=np.concatenate(F)[o].astype(np.float32), e=np.concatenate(E)[o].astype(np.float32), s=np.concatenate(S)[o])
        res.update(status='ok', npts=int(len(t)))
    except Exception as ex:
        res['status'] = f'err: {ex}'[:200]
    finally:
        shutil.rmtree(d, ignore_errors=True)
    json.dump(res, open(meta, 'w'))
    return tic


if __name__ == '__main__':
    tg = pd.read_csv(sys.argv[1])
    if len(sys.argv) > 2: tg = tg.head(int(sys.argv[2]))
    jobs = [(int(i), json.loads(p)) for i, p in zip(tg.ID, tg.plan)]
    with mp.Pool(int(os.environ.get('NPROC', 8))) as pool:
        for i, _ in enumerate(pool.imap(one, jobs, chunksize=1)):
            if i % 200 == 0: print(i, flush=True)
    print('done')
