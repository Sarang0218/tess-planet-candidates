"""Build injection-based training set from real flattened light curves (one npz per star)."""
import sys, glob, os, numpy as np, pandas as pd, batman, multiprocessing as mp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from views import make_views, DUR_GRID
TG = pd.read_csv('/home/sarang/tess_search/data/targets.csv').set_index('ID')
OUT = '/home/sarang/tess_search/out'
NP, NEB, NBL, NN = 12, 6, 4, 10

def bm(t, P, t0, k, aR, b, u):
    p = batman.TransitParams(); p.t0, p.per, p.rp, p.a = t0, P, k, aR
    p.inc = np.degrees(np.arccos(min(b / aR, 1))); p.ecc, p.w, p.u, p.limb_dark = 0., 90., list(u), 'quadratic'
    return batman.TransitModel(p, t).light_curve(p)

def snap(d): return DUR_GRID[np.argmin(np.abs(np.log(DUR_GRID / max(d, 1e-3))))]

def t14(model, t, P, t0):
    ph = ((t - t0) / P + 0.5) % 1 - 0.5
    i = model < 1 - 1e-7
    return (np.ptp(ph[i]) * P) if i.sum() > 2 else 0.02

def depth_for_snr(snr, mad, n_in): return snr * mad / np.sqrt(max(n_in, 1))

def one_star(arg):
    fn, rep = arg
    tic = int(os.path.basename(fn)[:-4]); rng = np.random.default_rng([tic, rep])
    d = np.load(fn); t, f = d['t'], d['f']; B = np.ptp(t)
    if len(t) < 2000 or B < 20: return []
    st = TG.loc[tic]; Ms, Rs = float(st['mass']), float(st['rad'])
    if not (np.isfinite(Ms) and np.isfinite(Rs)): Ms, Rs = 0.4, 0.4
    mad = 1.4826 * np.median(np.abs(f - np.median(f))); cad = np.median(np.diff(t))
    pmax = min(40., B / 3.); out = []
    def aR(P, M=Ms): return 4.2084 * P ** (2 / 3) * M ** (1 / 3) / Rs
    def jitter(P, t0, d):
        return P * (1 + rng.normal(0, d / (6 * B))), t0 + rng.normal(0, d / 8), snap(d * rng.uniform(.8, 1.25))
    def add(tt, ff, P, t0, d, y, kind):
        G, L, S = make_views(tt, ff, P, t0, d, Ms, Rs); out.append((G, L, S, y, kind, tic))
    u = (rng.uniform(.15, .45), rng.uniform(.2, .45))
    # --- planets
    for _ in range(NP):
        P = np.exp(rng.uniform(np.log(.5), np.log(pmax))); t0 = t.min() + rng.uniform(0, P)
        a = aR(P); b = rng.uniform(0, .9)
        m0 = bm(t, P, t0, .05, a, b, u); T = t14(m0, t, P, t0)
        n_in = T / cad * B / P
        dep = depth_for_snr(np.exp(rng.uniform(np.log(7), np.log(60))), mad, n_in)
        k = np.sqrt(dep / .9)
        if k > .3: continue
        m = bm(t, P, t0, k, a, b, u); add(t, f * m, *jitter(P, t0, T), 1, 'planet')
    # --- eclipsing binaries (diluted), shown at true P or at P/2
    for _ in range(NEB):
        P = np.exp(rng.uniform(np.log(1.), np.log(min(80., B / 2)))); t0 = t.min() + rng.uniform(0, P)
        k = rng.uniform(.3, 1.); a = aR(P, Ms * (1 + k ** 2)); b = rng.uniform(0, 1 + k * .8)
        prim = 1 - bm(t, P, t0, k, a, min(b, a - 1), u)
        if prim.max() <= 0: continue
        ecc_off = .5 if rng.random() < .7 else rng.uniform(.38, .62)
        sec = 1 - bm(t, P, t0 + ecc_off * P, k, a, min(b, a - 1), u)
        r = rng.uniform(.05, 1.); shape = prim + r * sec * prim.max() / max(sec.max(), 1e-9)
        T = t14(1 - prim, t, P, t0); n_in = T / cad * B / P
        dep = depth_for_snr(np.exp(rng.uniform(np.log(7), np.log(80))), mad, n_in)
        ff = f * (1 - dep * shape / prim.max())
        if r > .6 and ecc_off == .5 and rng.random() < .6: add(t, ff, *jitter(P / 2, t0, T), 0, 'eb_half')
        else: add(t, ff, *jitter(P, t0, T), 0, 'eb')
    # --- background blends: big companion, grazing V-shape, no visible secondary
    for _ in range(NBL):
        P = np.exp(rng.uniform(np.log(.5), np.log(pmax))); t0 = t.min() + rng.uniform(0, P)
        k = rng.uniform(.3, .9); a = aR(P) * rng.uniform(.5, 2); b = rng.uniform(1 - k * .3, 1 + k * .7)
        m0 = bm(t, P, t0, k, a, min(b, a - 1), u); s = 1 - m0
        if s.max() <= 0: continue
        T = t14(m0, t, P, t0); n_in = T / cad * B / P
        dep = depth_for_snr(np.exp(rng.uniform(np.log(7), np.log(60))), mad, n_in)
        add(t, f * (1 - dep * s / s.max()), *jitter(P, t0, T), 0, 'blend')
    # --- false alarms: deepest dip at a random period (mini-BLS), incl. inverted flux
    for j in range(NN):
        ff = f if j % 2 == 0 else 2 - f
        P = np.exp(rng.uniform(np.log(.5), np.log(pmax))); dd = rng.choice(DUR_GRID)
        ph = (t % P); nb = max(int(P / dd), 4); idx = (ph / P * nb).astype(int)
        s = np.bincount(idx, ff - 1, nb) / np.maximum(np.bincount(idx, None, nb), 1)
        cnt = np.bincount(idx, None, nb); s[cnt < 3] = 0
        kb = np.argmin(s); add(t, ff, P, t.min() - t.min() % P + (kb + .5) * P / nb, dd, 0, 'noise' if j % 2 == 0 else 'inverted')
    return out

if __name__ == '__main__':
    files = sorted(glob.glob(f'{OUT}/*.npz'))
    if len(sys.argv) > 2: files = files[:int(sys.argv[2])]
    with mp.Pool(int(os.environ.get('NPROC', 4))) as pool:
        args = [(f, r) for f in files for r in range(int(os.environ.get('REPS', 1)))]
        res = [r for rr in pool.imap_unordered(one_star, args, chunksize=2) for r in rr]
    G, L, S, y, kind, tic = map(np.array, zip(*res))
    np.savez_compressed(sys.argv[1], G=G, L=L, S=S, y=y, kind=kind, tic=tic)
    print(len(y), dict(zip(*np.unique(kind, return_counts=True))))
