"""Streaming GPU searcher for out3/: waits for preprocessed stars, rebuilds the common-mode bad-time mask
from every star seen so far (out/ + out3/), searches new stars in batches. Resumable. Env: K, NK (shard)."""
import sys, os, glob, json, time, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import search_core as sc
D = '/home/sarang/tess_search'; OUT3 = f'{D}/out3'; W = 10 / 1440
K, NK = int(os.environ.get('K', 0)), int(os.environ.get('NK', 1))


def build_bad(files):
    tot, low = {}, {}
    for fn in files:
        try: d = np.load(fn)
        except Exception: continue
        t, f = d['t'], d['f']
        z = (f - np.median(f)) / (1.4826 * np.median(np.abs(f - np.median(f))) + 1e-12)
        b = np.floor(t / W).astype(np.int64)
        for arr, st in ((b, tot), (b[z < -3], low)):
            u, c = np.unique(arr, return_counts=True)
            for k, v in zip(u, c): st[k] = st.get(k, 0) + 1
    keys = np.array(sorted(tot)); n = np.array([tot[k] for k in keys]); l = np.array([low.get(k, 0) for k in keys])
    bad = keys[(n >= 20) & (l >= 4) & (l / np.maximum(n, 1) > 0.03)]
    return np.unique(np.concatenate([bad - 1, bad, bad + 1])) if len(bad) else np.array([], np.int64)


last_n, bad = 0, None
while True:
    ready = sorted(glob.glob(f'{OUT3}/*.npz'))
    todo = [f for f in ready if int(os.path.basename(f)[:-4]) % NK == K and not os.path.exists(f[:-4] + '.res.json')]
    fetch_done = os.path.exists(f'{D}/fetch3.done')
    if not todo:
        if fetch_done: break
        time.sleep(30); continue
    if bad is None or len(ready) > 1.5 * last_n:          # refresh the mask as the sample grows
        bad = build_bad(glob.glob(f'{D}/out/*.npz') + ready); last_n = len(ready)
        print(f'[{time.strftime("%H:%M")}] bad-time mask from {len(ready)} stars: {len(bad)} bins', flush=True)
    for fn in todo[:200]:
        tic = os.path.basename(fn)[:-4]
        try:
            d = np.load(fn); res, t, f = sc.run_star(d['t'], d['f'], bad)
        except Exception as e:
            res = dict(status=f'err: {e}'[:200])
        res['tic'] = tic
        json.dump(res, open(fn[:-4] + '.res.json', 'w'))
print('search3 done')
