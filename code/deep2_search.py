"""Streaming GPU searcher for the deep pass (deep2/lc). Searches stars in target-priority order as they arrive. Resumable."""
import sys, os, glob, json, time, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import deep2_core as dc
dc.limit_gpu(2.0)
D = '/home/sarang/tess_search'; LC = f'{D}/deep2/lc'; W = 10 / 1440
order = {int(t): i for i, t in enumerate(pd.read_csv(f'{D}/deep2/targets.csv', usecols=['ID']).ID)}
old_bad = np.load(f'{D}/data/bad_bins.npy')


def build_bad(files):
    tot, low = {}, {}
    for fn in files:
        try: d = np.load(fn); t, f, e = d['t'], d['f'], d['e']
        except Exception: continue
        b = np.floor(t / W).astype(np.int64); z = (f - 1) / e
        for arr, st in ((b, tot), (b[z < -3], low)):
            u, c = np.unique(arr, return_counts=True)
            for k in u: st[k] = st.get(k, 0) + 1
    keys = np.array(sorted(tot)); n = np.array([tot[k] for k in keys]); l = np.array([low.get(k, 0) for k in keys])
    bad = keys[(n >= 20) & (l >= 4) & (l / np.maximum(n, 1) > 0.03)]
    bad = np.unique(np.concatenate([bad - 1, bad, bad + 1])) if len(bad) else np.array([], np.int64)
    return np.union1d(bad, old_bad)


NB = 10
last_n, bad = 0, old_bad
while True:
    ready = glob.glob(f'{LC}/*.npz')
    todo = sorted((f for f in ready if not os.path.exists(f[:-4] + '.res.json')), key=lambda f: order.get(int(os.path.basename(f)[:-4]), 10 ** 9))
    if not todo:
        if os.path.exists(f'{D}/deep2/fetch.done'): break
        time.sleep(60); continue
    if len(ready) >= 1500 and len(ready) > 1.5 * last_n:
        bad = build_bad(ready); last_n = len(ready); np.save(f'{D}/deep2/bad_bins.npy', bad)
        print(f'[{time.strftime("%H:%M")}] bad-time mask from {len(ready)} stars: {len(bad)} bins', flush=True)
    t0 = time.time()
    for fn in todo[:NB]:
        tic = os.path.basename(fn)[:-4]
        try:
            d = np.load(fn); res, t, f, e, s = dc.run_star(d['t'], d['f'], d['s'], bad)
        except Exception as ex:
            res = dict(status=f'err: {ex}'[:200])
        res['tic'] = tic
        json.dump(res, open(fn[:-4] + '.res.json', 'w'))
    print(f'[{time.strftime("%H:%M")}] searched {min(NB, len(todo))} stars in {time.time() - t0:.0f} s; {len(todo) - min(NB, len(todo))} waiting', flush=True)
print('deep2 search done')
