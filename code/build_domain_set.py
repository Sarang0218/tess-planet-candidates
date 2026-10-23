"""Labelled + unlabelled vetting examples from our own search (out3 signals), in exactly the processing used for candidates.
Labels: 1 = matches TOI CP/KP or confirmed planet (PC -> 1 with weight 0.6); 0 = matches TOI FP/FA, EB-catalogue star, or hard physical flag.
Unflagged unknown signals are kept unlabelled (y = -1) for scoring. Writes data/domain_set.npz."""
import os, sys, numpy as np, pandas as pd, multiprocessing as mp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from views import make_views
import search_core
D = '/home/sarang/tess_search'; BAD = np.load(f'{D}/data/bad_bins.npy')
o = pd.read_csv(f'{D}/out3_signals.csv', low_memory=False)
o = o[(o.snr >= 6) & (o.sde >= 6)].copy()
toi = pd.read_csv(f'{D}/data/toi.csv')
disp = {}
for tic, P, dsp in zip(toi['TIC ID'], toi['Period (days)'], toi['TFOPWG Disposition']):
    if np.isfinite(P) and P > 0: disp.setdefault(int(tic), []).append((float(P), str(dsp)))
def label(r):
    fl = set(str(r['flags']).split(';')) if isinstance(r['flags'], str) else set()
    if isinstance(r.known_match, str):
        if 'PS:' in r.known_match: return 1, 1.0
        for P, dsp in disp.get(int(r.tic), []):
            if any(abs(r.P / (P * m) - 1) < .01 for m in (.5, 1, 2)):
                if dsp in ('CP', 'KP'): return 1, 1.0
                if dsp in ('FP', 'FA'): return 0, 1.0
                if dsp in ('PC', 'APC'): return 1, 0.6
        return -1, 0.0
    if r.eb: return 0, 1.0
    hard = {'secondary', 'odd-even', 'deep', 'Rp>15', 'oot-variability', 'rotation-harmonic', 'P-shared', 'syst-period', 'too-long', 'inconsistent-depths'}
    if fl & hard: return 0, 0.8
    if fl: return -1, 0.0            # only soft flags: unknown
    return -1, 0.0
lab = [label(r) for _, r in o.iterrows()]; o['y'] = [a for a, b in lab]; o['w'] = [b for a, b in lab]
print(o.y.value_counts().to_string())
def one(g):
    tic, rows = g
    try:
        d = np.load(f'{D}/out3/{int(tic)}.npz'); t, f, _, _ = search_core.clean(d['t'], d['f'], BAD)
    except Exception: return []
    out = []
    for _, r in rows.iterrows():
        ms = r.Ms if np.isfinite(r.Ms) and r.Ms > 0 else 0.4; rs = r.Rs if np.isfinite(r.Rs) and r.Rs > 0 else 0.4
        G, L, S = make_views(t, f, r.P, r.t0, r.dur, ms, rs)
        sc = np.array([np.log10(max(r.snr, 1)), np.log10(max(r.sde, 1)), np.log10(max(r.ntr, 1)), r.top2, r.medr, min(r.oe, 10) / 10, np.clip(r.sec, -10, 10) / 10,
                       np.log10(max(r.dep, 1)), np.log10(max(r.Rp, .1)) if np.isfinite(r.Rp) else 0., min(r.contam, 5) if np.isfinite(r.contam) else 0.], np.float32)
        out.append((G, L, S, sc, int(r.y), float(r.w), int(r.tic), float(r.P), int(r.k)))
    return out
if __name__ == '__main__':
    groups = list(o.groupby('tic'))
    with mp.Pool(10) as pool:
        res = [x for rr in pool.imap_unordered(one, groups, chunksize=8) for x in rr]
    G, L, S, SC, y, w, tic, P, k = map(np.array, zip(*res))
    np.savez_compressed(f'{D}/data/domain_set.npz', G=G, L=L, S=S, SC=SC, y=y, w=w, tic=tic, P=P, k=k)
    print('saved', len(y), dict(zip(*np.unique(y, return_counts=True))))
