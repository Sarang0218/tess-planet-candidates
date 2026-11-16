"""Stage-1 vetting of a deep-pass signal from the prepared all-sector light curve (no downloads):
per-sector depths, odd/even, secondary and other-phase dips, depth at 2P, trapezoid shape, duration vs stellar density, figure.
Usage: deep2_vet.py TIC K [TIC K ...]   (K = signal index in deep2_signals.csv)"""
import sys, os, json, numpy as np, pandas as pd, warnings, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.optimize import least_squares
warnings.filterwarnings('ignore'); sys.path.insert(0, '/home/sarang/tess_search/scripts'); import deep2_core as dc
R = '/home/sarang/tess_search'; dc.limit_gpu(1.0)
bad = np.load(f'{R}/deep2/bad_bins.npy'); sig = pd.read_csv(f'{R}/deep2_signals.csv', low_memory=False)


def depth(t, f, e, P, T0, d, shift=0.0):
    ph = (((t - T0) / P + .5 + shift) % 1 - .5) * P; a = np.abs(ph) < d / 2; b = (np.abs(ph) > 0.75 * d) & (np.abs(ph) < 3 * d)
    if a.sum() < 3 or b.sum() < 5: return np.nan, np.nan
    return (f[b].mean() - f[a].mean()) * 1e6, np.sqrt(np.mean(e[a] ** 2) / a.sum() + np.mean(e[b] ** 2) / b.sum()) * 1e6


def trap(ph, dep, T14, q):                    # q = ingress fraction T12/T14 (0 box, 0.5 V)
    x = np.abs(ph); tin = T14 * q; flat = T14 / 2 - tin
    return 1 - dep * np.clip((T14 / 2 - x) / max(tin, 1e-6), 0, 1) if tin > 0 else 1 - dep * (x < T14 / 2)


def vet(tic, k, sig=None):
    sig = globals()['sig'] if sig is None else sig
    g = sig[(sig.tic == tic) & (sig.k == k)].iloc[0]; P, T0, d = float(g.P), float(g.t0), float(g.dur)
    z = np.load(f'{R}/deep2/lc/{tic}.npz'); t, f, e, s, _ = dc.prepare(z['t'], z['f'], z['s'], bad)
    # mask other strong signals of this star
    for _, o in sig[(sig.tic == tic) & (sig.k != k) & (sig.snr >= 7.5) & (sig.ntr >= 3)].iterrows():
        if abs(o.P / P - 1) < 0.01 or abs(o.P / P - 0.5) < 0.01 or abs(o.P / P - 2) < 0.01: continue
        keep = np.abs(((t - o.t0) / o.P + .5) % 1 - .5) * o.P > 0.75 * o.dur; t, f, e, s = t[keep], f[keep], e[keep], s[keep]
    # refine T0 and duration on a fine grid
    best = (-1, T0, d)
    for dd in d * np.array([0.5, 0.65, 0.8, 1.0, 1.25, 1.5]):
        for dt in np.linspace(-0.6 * d, 0.6 * d, 25):
            a, b = depth(t, f, e, P, T0 + dt, dd)
            if np.isfinite(a) and a / b > best[0]: best = (a / b, T0 + dt, dd)
    _, T0, d = best; d0, e0 = depth(t, f, e, P, T0, d)
    secs = np.unique(s); ps = {int(u): depth(t[s == u], f[s == u], e[s == u], P, T0, d) for u in secs}; v = np.array([x for x in ps.values() if np.isfinite(x[0])])
    w = 1 / v[:, 1] ** 2; dm = np.sum(v[:, 0] * w) / w.sum(); chi2 = np.sum(w * (v[:, 0] - dm) ** 2); npos = int((v[:, 0] > 0).sum())
    n = np.floor((t - T0) / P + .5).astype(int); od = depth(t[n % 2 == 1], f[n % 2 == 1], e[n % 2 == 1], P, T0, d); ev = depth(t[n % 2 == 0], f[n % 2 == 0], e[n % 2 == 0], P, T0, d)
    oe = abs(od[0] - ev[0]) / np.hypot(od[1], ev[1]); sec = depth(t, f, e, P, T0, d, 0.5)
    oo = np.abs(((t - T0) / P + .5) % 1 - .5) * P > 0.75 * d            # scan the other phases with the transit itself removed, so its wings do not register as structure
    sh = np.linspace(0.06, 0.94, 221); sc = np.array([depth(t[oo], f[oo], e[oo], P, T0, d, x) for x in sh]); zsc = sc[:, 0] / sc[:, 1]
    p2a, p2b = depth(t, f, e, 2 * P, T0, d), depth(t, f, e, 2 * P, T0 + P, d)
    # trapezoid fit on the fold
    ph = (((t - T0) / P + .5) % 1 - .5) * P; m = np.abs(ph) < 3 * d
    fit = least_squares(lambda p: (f[m] - trap(ph[m] - p[3], p[0], p[1], np.clip(p[2], 0.02, 0.5))) / e[m], [d0 * 1e-6, d * 1.1, 0.2, 0.0], bounds=([0, d * 0.3, 0.02, -d], [1, d * 3, 0.5, d]))
    dep_t, T14, q, dt = fit.x
    Rs, Ms = float(g.Rs), float(g.Ms); aR = 4.2084 * P ** (2 / 3) * Ms ** (1 / 3) / Rs; rp = np.sqrt(max(dep_t, 1e-9)); Tmax = P / np.pi * np.arcsin(min(1, (1 + rp) / aR))
    ntr = len(np.unique(n[np.abs(ph) < d / 2]))
    # per-transit depths
    pt = []
    for kk in np.unique(n[np.abs(ph) < d / 2]):
        mm = n == kk; a, b = depth(t[mm], f[mm], e[mm], P, T0, d)
        if np.isfinite(a): pt.append((a, b))
    pt = np.array(pt); zt = pt[:, 0] / pt[:, 1]
    out = dict(tic=int(tic), k=int(k), P=P, T0=float(T0), dur_h=float(d * 24), depth_ppm=float(d0), depth_err=float(e0), snr_box=float(d0 / e0), nsec=int(len(v)), nsec_pos=npos, chi2_sec=float(chi2),
               odd=od, even=ev, oe_sig=float(oe), sec_ppm=sec, sec_sig=float(sec[0] / sec[1]), max_other_phase_sig=float(np.nanmax(zsc)), min_other_phase_sig=float(np.nanmin(zsc)), scan_std=float(np.nanstd(zsc)),
               dep2P=(p2a, p2b), trap_depth_ppm=float(dep_t * 1e6), T14_h=float(T14 * 24), T14_max_h=float(Tmax * 24), ingress_frac=float(np.clip(q, 0.02, 0.5)), Rp_Re=float(rp * Rs * 109.2), aR=float(aR),
               ntr=int(ntr), frac_tr_pos=float(np.mean(zt > 0)), med_tr_z=float(np.median(zt)), max_tr_z=float(zt.max()), Tmag=float(g.Tmag), contam=float(g.contam), Rs=Rs)
    json.dump(out, open(f'{R}/deep2/vet/{tic}_{k}.json', 'w'), default=float)
    print(f"\nTIC {tic} k={k}: P={P:.6f} d  T0={T0:.4f}  box {d*24:.2f} h  depth {d0:.0f} ± {e0:.0f} ppm ({d0/e0:.1f} sigma)  Tmag {g.Tmag:.2f} R*={Rs:.2f} contam={g.contam:.3f}")
    print(f"   sectors +{npos}/{len(v)}, chi2 {chi2:.1f}/{len(v)-1};  odd {od[0]:.0f}±{od[1]:.0f} even {ev[0]:.0f}±{ev[1]:.0f} ({oe:.1f} sig);  secondary {sec[0]:.0f}±{sec[1]:.0f} ({sec[0]/sec[1]:.1f} sig)")
    print(f"   other phases: z std {np.nanstd(zsc):.2f}, max {np.nanmax(zsc):.1f}, min {np.nanmin(zsc):.1f};  at 2P: {p2a[0]:.0f}±{p2a[1]:.0f} / {p2b[0]:.0f}±{p2b[1]:.0f}")
    print(f"   trapezoid: depth {dep_t*1e6:.0f} ppm, T14 {T14*24:.2f} h (central-transit max {Tmax*24:.2f} h), ingress frac {np.clip(q,0.02,0.5):.2f}, Rp {rp*Rs*109.2:.2f} Re, a/R* {aR:.1f}")
    print(f"   transits: {ntr}, {np.mean(zt>0)*100:.0f}% positive, median z {np.median(zt):.2f}, max z {zt.max():.1f}")
    fig, ax = plt.subplots(1, 4, figsize=(19, 4)); phh = ph * 24
    for a, lim, nb in ((ax[0], P * 12, 80), (ax[1], 4 * d * 24, 40)):
        mm = np.abs(phh) < lim; ed = np.linspace(-lim, lim, nb + 1); idx = np.digitize(phh[mm], ed)
        a.plot(phh[mm], (f[mm] - 1) * 1e6, '.', color='.85', ms=1, rasterized=True)
        y = np.array([(f[mm][idx == j].mean() - 1) * 1e6 if (idx == j).sum() > 2 else np.nan for j in range(1, nb + 1)]); ye = np.array([f[mm][idx == j].std() / np.sqrt((idx == j).sum()) * 1e6 if (idx == j).sum() > 2 else np.nan for j in range(1, nb + 1)])
        a.errorbar((ed[:-1] + ed[1:]) / 2, y, ye, fmt='ok', ms=3); a.set_ylim(np.nanmin(y) - 3 * np.nanmedian(ye) - 0.3 * d0, np.nanmax(y) + 3 * np.nanmedian(ye) + 0.3 * d0); a.set_xlabel('hours from mid-transit'); a.set_ylabel('ppm')
    xx = np.linspace(-4 * d, 4 * d, 400); ax[1].plot(xx * 24, (trap(xx - dt, dep_t, T14, np.clip(q, 0.02, 0.5)) - 1) * 1e6, 'r-')
    ax[0].set_title(f'TIC {tic}  P = {P:.5f} d  full phase'); ax[1].set_title(f'depth {d0:.0f} ppm, T14 {T14*24:.2f} h, {ntr} transits')
    ss = sorted(ps); vv = np.array([ps[u] for u in ss]); ax[2].errorbar(range(len(ss)), vv[:, 0], vv[:, 1], fmt='ok'); ax[2].axhline(0, color='k', lw=.5); ax[2].axhline(dm, color='r', ls='--')
    ax[2].set_xticks(range(len(ss))); ax[2].set_xticklabels([f'S{u}' for u in ss], rotation=60, fontsize=7); ax[2].set_title('depth per sector (ppm)')
    for lab, mk, c in (('odd', n % 2 == 1, 'C0'), ('even', n % 2 == 0, 'C1')):
        mm = mk & (np.abs(phh) < 3 * d * 24); ed = np.linspace(-3 * d * 24, 3 * d * 24, 22); idx = np.digitize(phh[mm], ed)
        ax[3].errorbar((ed[:-1] + ed[1:]) / 2, [(f[mm][idx == j].mean() - 1) * 1e6 if (idx == j).sum() > 2 else np.nan for j in range(1, 22)], [f[mm][idx == j].std() / np.sqrt(max((idx == j).sum(), 1)) * 1e6 for j in range(1, 22)], fmt='o-', color=c, ms=3, label=lab)
    ax[3].legend(); ax[3].set_title(f'odd vs even ({oe:.1f} sigma); secondary {sec[0]/sec[1]:.1f} sigma'); ax[3].set_xlabel('hours')
    plt.tight_layout(); plt.savefig(f'{R}/deep2/vet/{tic}_{k}.png', dpi=75); plt.close()
    return out


if __name__ == '__main__':
    for i in range(1, len(sys.argv), 2):
        try: vet(int(sys.argv[i]), int(sys.argv[i + 1]))
        except Exception as ex: print('FAIL', sys.argv[i], ex)
