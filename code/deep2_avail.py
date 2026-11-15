"""Deep pass 2: combined sector availability (SPOC 2-min + TESS-SPOC FFI) for every star already searched,
and the list of stars whose total data is much larger than what the first pass used."""
import pandas as pd, numpy as np, json, glob, os
D = '/home/sarang/tess_search'
a2 = pd.read_pickle(f'{D}/data/avail2min.pkl')
tg = pd.concat([pd.read_csv(f'{D}/data/{n}.csv') for n in ('targets_allsky', 'targets_ffi_faint', 'targets')]).drop_duplicates('ID')
ids = set(tg.ID.astype(int))
ffi = {}
for fn in sorted(glob.glob(f'{D}/data/spoc_targets/s*.csv')):
    s = int(os.path.basename(fn)[1:5])
    x = pd.read_csv(fn, usecols=[0], comment=None); x = x.iloc[:, 0].astype(str).str.lstrip('#').astype(np.int64)
    for t in np.intersect1d(x.values, np.fromiter(ids, np.int64)): ffi.setdefault(int(t), []).append(s)
rows = []
for t in ids:
    try: m = json.load(open(f'{D}/out3/{t}.meta.json'))
    except Exception: m = {}
    used = m.get('used', m.get('sectors', [])) if m.get('status') == 'ok' else []
    tm = {s: fn for s, fn in a2.get(t, [])}; fs = set(ffi.get(t, []))
    plan = {}
    for s in sorted(set(tm) | fs):
        if s in tm and s not in fs: plan[s] = 'tm:' + tm[s]      # FFI product preferred: 10-25x smaller download
        else: plan[s] = 'ffi'
    rows.append(dict(ID=t, nused=len(used), ntot=len(plan), n2=len(tm), nffi=len(fs), span=(max(plan) - min(plan) + 1) if plan else 0, plan=json.dumps(plan)))
r = pd.DataFrame(rows).merge(tg[['ID', 'Tmag', 'Teff', 'rad', 'mass', 'contratio', 'd']], on='ID')
r.to_pickle(f'{D}/deep2/avail_all.pkl')
print(len(r), 'stars;  searched before:', (r.nused > 0).sum())
print(r[['nused', 'ntot', 'n2', 'nffi', 'span']].describe().round(1))
for k in (1.5, 2, 3, 4):
    m = (r.ntot >= k * r.nused.clip(lower=1)) & (r.ntot >= 4)
    print(f'ntot >= {k} x nused and >= 4: {m.sum()} stars, {r.ntot[m].sum()} sector files, median ntot {r.ntot[m].median()}')
