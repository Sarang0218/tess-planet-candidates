"""Noise-floor test for a deep-pass signal: the identical search on (a) the sign-inverted cleaned light curve and (b) N copies in which
each sector's cleaned flux is shifted cyclically by a random offset (keeps the noise and any one-off features, destroys multi-sector
phase coherence).  Reports the strongest PERIODIC (multi-transit) peak of each run.   Usage: deep2_null.py TIC [N]"""
import sys, os, numpy as np
sys.path.insert(0, '/home/sarang/tess_search/scripts'); import deep2_core as dc
dc.limit_gpu(2.0); bad = np.load('/home/sarang/tess_search/deep2/bad_bins.npy')
tic = sys.argv[1]; N = int(sys.argv[2]) if len(sys.argv) > 2 else 20; d = np.load(f'/home/sarang/tess_search/deep2/lc/{tic}.npz')
t, f, e, s, _ = dc.prepare(d['t'], d['f'], d['s'], bad)
def top(res):
    m = [g for g in res['signals'] if dc.is_multi(g)]; return max(m, key=lambda g: g['snr']) if m else dict(P=np.nan, snr=0.0, sde=0.0, dur=0, depth_ppm=0)
def show(lab, res): print(lab + ': ' + ' | '.join(f"P={g['P']:.4f} dep={g['depth_ppm']:.0f} snr={g['snr']:.1f} sde={g['sde']:.1f} ntr={g['ntr']}{'' if dc.is_multi(g) else ' (one-off)'}" for g in res['signals']), flush=True)
real = dc.search_prepared(t, f, e, s, nsig=2); show('real    ', real); inv = dc.search_prepared(t, 2 - f, e, s, nsig=1); show('inverted', inv)
# remove the real periodic signals first, so the scrambles measure the noise and not a partly re-aligned signal
keep = np.ones(len(t), bool)
for g in real['signals']:
    if dc.is_multi(g) and g['snr'] >= 7.5: keep &= np.abs(((t - g['t0']) / g['P'] + .5) % 1 - .5) * g['P'] > 0.75 * g['dur']
t, f, e, s = t[keep], f[keep], e[keep], s[keep]; print(f'masked {int((~keep).sum())} in-transit points of the real signals before scrambling')
rng = np.random.default_rng(1); tops = []
for it in range(N):
    f2 = f.copy()
    for u in np.unique(s):
        m = np.flatnonzero(s == u); f2[m] = np.roll(f[m], rng.integers(50, len(m) - 50))
    g = top(dc.search_prepared(t, f2, e, s, nsig=1)); tops.append((g['snr'], g['sde'], g['P']))
A = np.array(tops); print('scrambled, strongest periodic peak (snr, sde):', ' '.join(f'({a:.1f},{b:.1f},P={c:.3f})' for a, b, c in tops))
print(f'scrambled max SNR {A[:,0].max():.1f}, max SDE {A[:,1].max():.1f}; inverted {top(inv)["snr"]:.1f}; real {top(real)["snr"]:.1f} (SDE {top(real)["sde"]:.1f})')
