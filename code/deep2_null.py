"""Noise-floor test for a deep-pass signal: run the identical search on (a) the sign-inverted light curve and
(b) the light curve with known signals masked and sector-wise time-reversed blocks, and report the top SNR/SDE found.
Usage: deep2_null.py TIC"""
import sys, numpy as np
sys.path.insert(0, '/home/sarang/tess_search/scripts'); import deep2_core as dc
import os
bad = np.load('/home/sarang/tess_search/deep2/bad_bins.npy') if os.path.exists('/home/sarang/tess_search/deep2/bad_bins.npy') else np.load('/home/sarang/tess_search/data/bad_bins.npy')
tic = sys.argv[1]; d = np.load(f'/home/sarang/tess_search/deep2/lc/{tic}.npz')
def show(lab, res):
    print(lab + ': ' + ' | '.join(f"P={g['P']:.4f} dur={g['dur']*24:.2f}h dep={g['depth_ppm']:.0f} snr={g['snr']:.1f} sde={g['sde']:.1f} sec+={g.get('nsec_pos')}/{g.get('nsec_in')}" for g in res['signals']))
res, t, f, e, s = dc.run_star(d['t'], d['f'], d['s'], bad); show('real    ', res)
res_i, *_ = dc.run_star(d['t'], 2 - d['f'], d['s'], bad); show('inverted', res_i)
# scramble: keep each sector's flux series but shift it cyclically by a random offset (destroys multi-sector phase coherence, keeps the noise)
rng = np.random.default_rng(1); tops = []
for it in range(int(sys.argv[2]) if len(sys.argv) > 2 else 6):
    f2 = d['f'].copy()
    for u in np.unique(d['s']):
        m = np.flatnonzero(d['s'] == u); f2[m] = np.roll(d['f'][m], rng.integers(50, len(m) - 50))
    r2, *_ = dc.run_star(d['t'], f2, d['s'], bad, nsig=1); tops.append((r2['signals'][0]['snr'], r2['signals'][0]['sde']))
print('scrambled top (snr, sde):', ' '.join(f'({a:.1f},{b:.1f})' for a, b in tops))
