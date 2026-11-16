"""TRICERATOPS false-positive probability for a deep-pass signal, using the apertures of the TESS-SPOC pixel files already
downloaded by deep2_pix.py and the phase-folded all-sector light curve.   Usage: deep2_tric.py TIC K [nruns]"""
import sys, os, json, glob, numpy as np, pandas as pd, warnings, lightkurve as lk
warnings.filterwarnings('ignore')
R = '/home/sarang/tess_search'; TIC, K = int(sys.argv[1]), int(sys.argv[2]); NR = int(sys.argv[3]) if len(sys.argv) > 3 else 3
v = json.load(open(f'{R}/deep2/vet/{TIC}_{K}.json')); P, t0, dur, depth = v['P'], v['T0'], v['T14_h'] / 24, v['trap_depth_ppm'] * 1e-6
z = np.load(f'{R}/deep2/vet/{TIC}.prep.npz'); t, f, e, s = z['t'], z['f'], z['e'], z['s']        # cleaned light curve written by deep2_core.prepare (GPU environment)
sig = pd.read_csv(f'{R}/deep2_signals.csv', low_memory=False)
for _, o in sig[(sig.tic == TIC) & (sig.k != K) & (sig.snr >= 7.5) & (sig.ntr >= 3)].iterrows():
    if min(abs(o.P / P - q) for q in (1, .5, 2)) < 0.01: continue
    keep = np.abs(((t - o.t0) / o.P + .5) % 1 - .5) * o.P > 0.75 * o.dur; t, f, e = t[keep], f[keep], e[keep]
ph = (((t - t0) / P + .5) % 1 - .5) * P; m = np.abs(ph) < max(3 * dur, 0.12); ph, fl, er = ph[m], f[m], e[m]; o = np.argsort(ph); ph, fl, er = ph[o], fl[o], er[o]
nb = int(np.clip(len(ph) / 12, 30, 90)); ed = np.linspace(ph.min(), ph.max(), nb + 1); idx = np.digitize(ph, ed[1:-1])
tb = np.array([ph[idx == j].mean() for j in range(nb) if (idx == j).sum() > 2]); fb = np.array([np.average(fl[idx == j], weights=1 / er[idx == j] ** 2) for j in range(nb) if (idx == j).sum() > 2])
eb = float(np.median([np.sqrt(1 / np.sum(1 / er[idx == j] ** 2)) for j in range(nb) if (idx == j).sum() > 2])); fb = fb / np.median(fb[np.abs(tb) > dur])
aps, secs, tpp = [], [], {}
for p in sorted(glob.glob(f'{R}/deep2/pix/{TIC}/*_tp.fits')):
    tpf = lk.read(p); r_, c_ = np.where(tpf.pipeline_mask)
    if len(r_): aps.append(np.array([[tpf.column + cc, tpf.row + rr] for rr, cc in zip(r_, c_)])); secs.append(int(tpf.sector)); tpp[int(tpf.sector)] = p
import triceratops.triceratops as tr
import lightkurve as _lk                         # MAST's cutout service can be down: hand TRICERATOPS the pixel files already on disk (same WCS and pixel grid)
class _Local:
    def __init__(self, path): self.path = path
    def download_all(self, **kw): return [_lk.read(self.path)]
_lk.search_tesscut = lambda target=None, sector=None, **kw: _Local(tpp[int(sector)])
tg = tr.target(ID=TIC, sectors=np.array(secs)); tg.calc_depths(tdepth=depth, all_ap_pixels=aps)
print(f'TIC {TIC} k={K}: P={P:.5f} depth {depth*1e6:.0f} ppm; sectors {secs}; {len(tb)} phase bins, error {eb*1e6:.0f} ppm')
print(tg.stars[['ID', 'Tmag', 'sep (arcsec)', 'fluxratio', 'tdepth']].head(6).to_string()); res = []
for run in range(NR):
    tg.calc_probs(time=tb, flux_0=fb, flux_err_0=eb, P_orb=P, N=300000, parallel=True); res.append((float(tg.FPP), float(tg.NFPP))); print(f'run {run}: FPP = {tg.FPP:.4f}  NFPP = {tg.NFPP:.4f}', flush=True)
print(tg.probs[['ID', 'scenario', 'prob']].sort_values('prob', ascending=False).head(8).to_string())
A = np.array(res); print(f'TRICERATOPS FPP = {A[:,0].mean():.4f} ± {A[:,0].std():.4f}; NFPP = {A[:,1].mean():.4f} (max {A[:,1].max():.4f})')
json.dump(dict(FPP=A[:, 0].tolist(), NFPP=A[:, 1].tolist(), sectors=secs), open(f'{R}/deep2/vet/{TIC}_{K}.tric.json', 'w')); tg.probs.to_csv(f'{R}/deep2/vet/{TIC}_{K}.tric_probs.csv', index=False)
