"""AstroNet-style phase-folded views shared by training-set generation and inference."""
import numpy as np
NG, NL = 201, 81
DUR_GRID = np.array([0.02, 0.03, 0.045, 0.065, 0.09, 0.13])

def _binned(x, y, lo, hi, n):
    e = np.linspace(lo, hi, n + 1); i = np.digitize(x, e) - 1
    ok = (i >= 0) & (i < n); i, y = i[ok], y[ok]
    v = np.full(n, np.nan)
    if len(i):
        o = np.argsort(i); i, y = i[o], y[o]
        u, st = np.unique(i, return_index=True)
        for a, b, k in zip(st, list(st[1:]) + [len(i)], u): v[k] = np.median(y[a:b])
    g = np.isfinite(v)
    if g.sum() == 0: return np.zeros(n), 0.0
    if g.sum() < n: v[~g] = np.interp(np.flatnonzero(~g), np.flatnonzero(g), v[g])
    return v, g.mean()

def expected_dur(P, Ms, Rs):
    a = 4.2084 * P ** (2 / 3) * Ms ** (1 / 3) / Rs   # a/R*, P in days
    return P / np.pi * np.arcsin(min(1.0, 1.0 / a))

def make_views(t, f, P, t0, dur, Ms, Rs):
    ph = ((t - t0) / P + 0.5) % 1 - 0.5           # orbital phase, transit at 0
    n = np.floor((t - t0) / P + 0.5).astype(int)
    w = min(3 * dur / P, 0.5)
    g, _ = _binned(ph, f, -0.5, 0.5, NG)
    loc, cov = _binned(ph, f, -w, w, NL)
    odd, _ = _binned(ph[n % 2 == 1], f[n % 2 == 1], -w, w, NL)
    even, _ = _binned(ph[n % 2 == 0], f[n % 2 == 0], -w, w, NL)
    ph2 = ((t - t0) / P) % 1 - 0.5
    sec, _ = _binned(ph2, f, -w, w, NL)
    out_of = np.abs(np.linspace(-w, w, NL)) > 1.5 * dur / P
    base = np.median(loc[out_of]) if out_of.any() else np.median(g)
    oot = loc[out_of] if out_of.sum() > 5 else g
    noise = 1.4826 * np.median(np.abs(oot - np.median(oot)))
    depth = max(base - loc[NL // 2 - 2: NL // 2 + 3].mean(), 2 * noise, 1e-6)
    L = np.stack([loc, odd, even, sec]) - base
    L /= depth
    G = (g - np.median(g)) / depth
    scal = np.array([np.log10(dur / expected_dur(P, Ms, Rs)), np.log10(dur / P), cov], dtype=np.float32)
    return np.clip(G, -3, 3).astype(np.float32), np.clip(L, -3, 3).astype(np.float32), scal
