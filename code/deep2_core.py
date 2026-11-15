"""Deep pass 2 search core: multi-year, all-sector, inverse-variance-weighted GPU BLS.
Differences from search_core: per-sector noise weights; rotation prewhitening per contiguous season; period grid to 80 d;
coarse grid (phase drift <= 2/3 of the minimum duration) followed by a full-resolution local refinement."""
import numpy as np, torch
import search_core as sc
BINW = 0.02 / 3
WINS = np.array([0.02, 0.03, 0.045, 0.065, 0.09, 0.13])
W10 = 10 / 1440


def period_grid(B, pmin=0.5, pmax=80.0, coarse=2.0):
    pmax = min(pmax, B / 3.0); P, ps = pmin, []
    while P < pmax:
        ps.append(P); P += coarse * P * 0.03 * P ** (1 / 3) / (3.0 * B)
    return np.array(ps)


def bls_power(t, y, w, periods, dev='cuda', budget=3.0e7):
    """Weighted BLS. y: flux - baseline; w: weights in units where a typical point has w=1 (noise sig0 per unit weight).
    Returns per-period best SNR (in units of sig0), depth, t0 offset, duration."""
    t = torch.tensor(t - t.min(), dtype=torch.float32, device=dev)
    w_ = torch.tensor(w, dtype=torch.float32, device=dev); wy = torch.tensor(w * y, dtype=torch.float32, device=dev)
    N = len(y); chunk = int(max(32, min(1024, budget // N))); Wtot = float(w_.sum()); Stot = float(wy.sum())
    wb = np.maximum(np.round(WINS / BINW).astype(int), 1)
    out = {k: np.zeros(len(periods)) for k in ('snr', 'depth', 't0', 'dur')}
    for c in range(0, len(periods), chunk):
        Pc = periods[c:c + chunk]; nb = np.maximum(np.round(Pc / BINW).astype(int), 8); NB = int(nb.max()); C = len(Pc)
        Pt = torch.tensor(Pc, dtype=torch.float32, device=dev); nbt = torch.tensor(nb, device=dev)
        x = t[None, :] / Pt[:, None]; x = x - torch.floor(x)
        b = torch.clamp(torch.minimum(torch.floor(x * nbt[:, None].float()), nbt[:, None].float() - 1), min=0).long()
        b += (torch.arange(C, device=dev) * NB)[:, None]; b = b.ravel(); del x
        S = torch.zeros(C * NB, device=dev).scatter_add_(0, b, wy.repeat(C)).view(C, NB)
        Nn = torch.zeros(C * NB, device=dev).scatter_add_(0, b, w_.repeat(C)).view(C, NB); del b
        Wm = int(wb.max()); ar = torch.arange(NB + Wm, device=dev)[None, :]
        src = torch.remainder(ar, nbt[:, None])
        cS = torch.nn.functional.pad(torch.gather(S, 1, src).cumsum(1), (1, 0)); cN = torch.nn.functional.pad(torch.gather(Nn, 1, src).cumsum(1), (1, 0))
        best = torch.full((C,), -1e9, device=dev); bd = torch.zeros_like(best); bt = torch.zeros_like(best); bw = torch.zeros_like(best)
        valid = ar[:, :NB] < nbt[:, None]
        for wi in np.unique(wb):
            sin = cS[:, wi:wi + NB] - cS[:, :NB]; nin = cN[:, wi:wi + NB] - cN[:, :NB]; nout = Wtot - nin
            dep = -(sin / nin.clamp(min=1e-3) - (Stot - sin) / nout.clamp(min=1e-3))
            snr = dep / torch.sqrt(1 / nin.clamp(min=1e-3) + 1 / nout.clamp(min=1e-3))
            snr = torch.where(valid & (nin >= 3), snr, torch.full_like(snr, -1e9))
            m, j = snr.max(1); upd = m > best
            best = torch.where(upd, m, best); bd = torch.where(upd, dep.gather(1, j[:, None])[:, 0], bd)
            bt = torch.where(upd, (j.double() + wi / 2).float(), bt); bw = torch.where(upd, torch.full_like(bw, float(wi)), bw)
        out['snr'][c:c + chunk] = best.cpu().numpy(); out['depth'][c:c + chunk] = bd.cpu().numpy()
        out['t0'][c:c + chunk] = (bt.cpu().numpy() / nb) * Pc; out['dur'][c:c + chunk] = bw.cpu().numpy() * BINW
    return out


def refine(t, y, w, P, t0, dur, B, coarse=2.0):
    ph = ((t - t0) / P + .5) % 1 - .5; sel = np.abs(ph) * P < 0.3
    ts, ys, ws = t[sel], y[sel], w[sel]; Wtot, Stot = w.sum(), (w * y).sum(); best = (-1e9, P, t0, dur)
    tm = np.median(t)                                    # pivot the period scan about the middle of the data
    n0 = np.round((tm - t0) / P); tc = t0 + n0 * P
    for p in P * (1 + np.linspace(-1, 1, 61) * 1.5 * coarse * 0.03 * P ** (1 / 3) / (3 * B)):
        ph0 = (((ts - tc) / p + .5) % 1 - .5) * p
        for d in WINS:
            for dt in np.linspace(-0.75 * d, 0.75 * d, 9):
                a = np.abs(ph0 - dt) < d / 2; wa = ws[a].sum()
                if wa < 3: continue
                sa = (ws[a] * ys[a]).sum(); dd = -(sa / wa - (Stot - sa) / (Wtot - wa)); s = dd / np.sqrt(1 / wa + 1 / (Wtot - wa))
                if s > best[0]: best = (s, p, tc + dt, d)
    return best


def stats(t, y, w, P, t0, dur):
    """Weighted versions of the first-pass statistics; y in flux units, w = 1/sigma^2 (flux^-2)."""
    ph = ((t - t0) / P + .5) % 1 - .5; n = np.floor((t - t0) / P + .5).astype(int)
    intr = np.abs(ph) < dur / (2 * P); near = (np.abs(ph) > dur / P) & (np.abs(ph) < 3 * dur / P)
    def dep(a, b):
        if a.sum() < 2 or b.sum() < 3: return np.nan, np.nan
        return np.average(y[b], weights=w[b]) - np.average(y[a], weights=w[a]), np.sqrt(1 / w[a].sum() + 1 / w[b].sum())
    d_all = dep(intr, near); d_odd = dep(intr & (n % 2 == 1), near & (n % 2 == 1)); d_even = dep(intr & (n % 2 == 0), near & (n % 2 == 0))
    ph2 = ((t - t0) / P) % 1 - .5; d_half = dep(np.abs(ph2) < dur / (2 * P), (np.abs(ph2) > dur / P) & (np.abs(ph2) < 3 * dur / P))
    per = []
    for k in np.unique(n[intr]):
        a = intr & (n == k); b = near & (n == k)
        if a.sum() >= 2 and b.sum() >= 3:
            dk, ek = dep(a, b); per.append((dk, 1 / ek ** 2))
    per = np.array(per) if per else np.zeros((0, 2))
    if len(per) >= 2:
        chi = per[:, 0] ** 2 * per[:, 1]; top2 = np.sort(chi)[-2:].sum() / chi.sum()
        medr = np.median(per[:, 0]) / d_all[0] if d_all[0] > 0 else 0
    else: top2, medr = 1.0, 0.0
    return dict(depth_ppm=d_all[0] * 1e6, depth_e=d_all[1] * 1e6, odd=d_odd[0] * 1e6, odd_e=d_odd[1] * 1e6, even=d_even[0] * 1e6, even_e=d_even[1] * 1e6,
                half=d_half[0] * 1e6, half_e=d_half[1] * 1e6, ntr=int(len(per)), top2_frac=float(top2), med_ratio=float(medr))


def prepare(t, f, s, bad):
    """Bad-time mask, isolated-outlier cut, per-season rotation prewhitening, per-sector noise. Returns t, f, e, s, removed."""
    keep = ~np.isin(np.floor(t / W10).astype(np.int64), bad); t, f, s = t[keep], f[keep], s[keep]
    def sec_sig(f, s):
        e = np.zeros(len(f))
        for u in np.unique(s):
            m = s == u; e[m] = 1.4826 * np.median(np.abs(f[m] - np.median(f[m]))) + 1e-9
        return e
    def iso_cut(t, f, s):
        z = (f - 1) / sec_sig(f, s); zl, zr = np.r_[0, z[:-1]], np.r_[z[1:], 0]
        k = ~((z < -4) & (zl > -1.5) & (zr > -1.5)); return t[k], f[k], s[k]
    t, f, s = iso_cut(t, f, s); f = f.astype(float).copy(); removed = []
    edges = np.r_[0, np.flatnonzero(np.diff(t) > 30) + 1, len(t)]        # seasons: split at gaps > 30 d
    for a, b in zip(edges[:-1], edges[1:]):
        if b - a < 300: continue
        f[a:b], rem = sc.prewhiten(t[a:b], f[a:b]); removed += [(float(p), float(am)) for p, am in rem]
    t, f, s = iso_cut(t, f, s)
    for u in np.unique(s): f[s == u] /= np.median(f[s == u])
    return t, f, sec_sig(f, s), s, removed


def run_star(t, f, s, bad, nsig=3):
    t, f, e, s, removed = prepare(t, f, s, bad)
    res = dict(npts=int(len(t)), nsec=int(len(np.unique(s))), prewhitened=removed)
    if len(t) < 1000 or np.ptp(t) < 10:
        res['status'] = 'too_few'; return res, t, f, e, s
    res.update(search_prepared(t, f, e, s, nsig)); res['status'] = 'ok'
    return res, t, f, e, s


def search_prepared(t, f, e, s, nsig=3, stop_snr=6.0):
    """BLS iterations on an already cleaned light curve (used by run_star and by the injection / scrambled-null runs)."""
    res = {}
    B = np.ptp(t); periods = period_grid(B); sig0 = float(np.median(e)); w = (sig0 / e) ** 2; y = (f - 1) / sig0
    res['mad_ppm'] = sig0 * 1e6; res['baseline'] = float(B); res['nper'] = len(periods)
    m = np.ones(len(t), bool); sigs = []
    for k in range(nsig):
        r = bls_power(t[m], y[m], w[m], periods)
        p = r['snr']; med = np.median(p); mad = 1.4826 * np.median(np.abs(p - med)); j = int(np.argmax(p))
        snr, P, t0, dur = refine(t[m], y[m], w[m], periods[j], r['t0'][j] + t[m].min(), r['dur'][j], B)
        sg = dict(P=float(P), t0=float(t0), dur=float(dur), snr=float(snr), sde=float((p[j] - med) / mad), snr_coarse=float(p[j]))
        sg.update(stats(t[m], f[m], 1 / e[m] ** 2, P, t0, dur))
        ph = ((t - t0) / P + .5) % 1 - .5; intr = np.abs(ph) < dur / (2 * P); near = (np.abs(ph) > dur / P) & (np.abs(ph) < 4 * dur / P)
        ds = []
        for u in np.unique(s):
            a, b = intr & (s == u) & m, near & (s == u) & m
            if a.sum() >= 2 and b.sum() >= 5: ds.append((int(u), float(f[b].mean() - f[a].mean()), float(e[a][0] * np.sqrt(1 / a.sum() + 1 / b.sum()))))
        if ds:
            d_ = np.array(ds); wv = 1 / d_[:, 2] ** 2; h = len(d_) // 2
            def comb(x): return (float(np.sum(x[:, 1] / x[:, 2] ** 2) / np.sum(1 / x[:, 2] ** 2) * np.sqrt(np.sum(1 / x[:, 2] ** 2))) if len(x) else 0.0)
            sg.update(nsec_in=len(d_), nsec_pos=int((d_[:, 1] > 0).sum()), snr_h1=comb(d_[:h]), snr_h2=comb(d_[h:]),
                      sec_chi2=float(np.sum(wv * (d_[:, 1] - np.sum(wv * d_[:, 1]) / wv.sum()) ** 2)))
        sigs.append(sg)
        if snr < stop_snr: break
        m &= ~(np.abs(ph) < 0.75 * dur / P)
    res['signals'] = sigs
    return res
