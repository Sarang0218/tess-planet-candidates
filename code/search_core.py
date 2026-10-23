"""Shared search core: cleaning, GPU BLS iteration, refinement, per-signal vetting statistics."""
import numpy as np
import gpu_bls32 as g
import torch


def prewhiten(t, f, nmax=4, pmin=0.02, pmax=5.0, amp_snr=8.0):
    """Remove strong quasi-sinusoidal stellar variability (rotation) with GPU Lomb-Scargle:
    fit and subtract the top peak with 3 harmonics, repeat while the peak amplitude is significant.
    Returns cleaned flux and list of (period, amplitude) removed."""
    dev = 'cuda'; tt = torch.tensor(t - t.min(), dtype=torch.float32, device=dev)
    y = torch.tensor(f - np.median(f), dtype=torch.float32, device=dev)
    B = float(tt.max()); fr = torch.arange(1 / pmax, 1 / pmin, 1 / (5 * B), device=dev)
    sig = float(1.4826 * np.median(np.abs(f - np.median(f)))); removed = []
    for _ in range(nmax):
        pw = torch.zeros_like(fr)
        for i in range(0, len(fr), 2048):
            w = 2 * np.pi * fr[i:i + 2048, None] * tt[None, :]
            c, s_ = (y[None] * torch.cos(w)).mean(1), (y[None] * torch.sin(w)).mean(1)
            pw[i:i + 2048] = 2 * torch.sqrt(c ** 2 + s_ ** 2)            # approx sinusoid amplitude
        j = int(torch.argmax(pw)); amp = float(pw[j]); f0 = float(fr[j])
        if amp < amp_snr * sig / np.sqrt(len(f)) * 2 or amp < 3e-4: break
        cols = [torch.ones_like(tt)]
        for h in (1, 2, 3):
            cols += [torch.cos(2 * np.pi * h * f0 * tt), torch.sin(2 * np.pi * h * f0 * tt)]
        A = torch.stack(cols, 1)
        # spots evolve: fit amplitudes/phases separately per continuous segment (split at gaps > 0.5 d)
        edges = np.r_[0, np.flatnonzero(np.diff(t) > 0.5) + 1, len(t)]
        for a, b in zip(edges[:-1], edges[1:]):
            if b - a < 50: continue
            coef = torch.linalg.lstsq(A[a:b], y[a:b, None]).solution
            y[a:b] = y[a:b] - (A[a:b] @ coef)[:, 0]
        removed.append((1 / f0, amp))
    return y.cpu().numpy() + np.median(f), removed
W = 10 / 1440

def clean(t, f, bad):
    keep = ~np.isin(np.floor(t / W).astype(np.int64), bad)
    t, f = t[keep], f[keep]
    sig = 1.4826 * np.median(np.abs(f - np.median(f))); z = (f - np.median(f)) / sig
    # isolated low outliers: < -4 sigma with both neighbours > -1.5 sigma (a real transit spans >=3 bins)
    zl, zr = np.r_[0, z[:-1]], np.r_[z[1:], 0]
    iso = (z < -4) & (zl > -1.5) & (zr > -1.5)
    return t[~iso], f[~iso], int((~keep).sum()), int(iso.sum())

def stats(t, f, P, t0, dur, sig):
    ph = ((t - t0) / P + .5) % 1 - .5; n = np.floor((t - t0) / P + .5).astype(int)
    intr = np.abs(ph) < dur / (2 * P); near = (np.abs(ph) > dur / P) & (np.abs(ph) < 3 * dur / P)
    def dep(m):
        a, b = intr & m, near & m
        if a.sum() < 2 or b.sum() < 3: return np.nan, np.nan
        return np.median(f[b]) - f[a].mean(), sig * np.sqrt(1 / a.sum() + 1 / b.sum())
    d_all = dep(np.ones_like(intr)); d_odd = dep(n % 2 == 1); d_even = dep(n % 2 == 0)
    ph2 = ((t - t0) / P) % 1 - .5; i2 = np.abs(ph2) < dur / (2 * P); o2 = (np.abs(ph2) > dur / P) & (np.abs(ph2) < 3 * dur / P)
    d_half = (np.median(f[o2]) - f[i2].mean(), sig * np.sqrt(1 / max(i2.sum(), 1) + 1 / max(o2.sum(), 1))) if i2.sum() > 1 else (np.nan, np.nan)
    per = []
    for k in np.unique(n[intr]):
        a = intr & (n == k); b = near & (n == k)
        if a.sum() >= 2 and b.sum() >= 3: per.append((np.median(f[b]) - f[a].mean(), a.sum()))
    per = np.array(per) if per else np.zeros((0, 2))
    # fraction of total chi2 contributed by the top-2 epochs, and median per-transit depth / global depth
    if len(per) >= 2:
        chi = per[:, 0] ** 2 * per[:, 1]; top2 = np.sort(chi)[-2:].sum() / chi.sum()
        medr = np.median(per[:, 0]) / d_all[0] if d_all[0] > 0 else 0
    else: top2, medr = 1.0, 0.0
    return dict(depth_ppm=d_all[0] * 1e6, odd=d_odd[0] * 1e6, odd_e=d_odd[1] * 1e6, even=d_even[0] * 1e6, even_e=d_even[1] * 1e6,
                half=d_half[0] * 1e6, half_e=d_half[1] * 1e6, ntr=int(len(per)), top2_frac=float(top2), med_ratio=float(medr))

def refine(t, f, P, t0, dur):
    # local fine grid around the peak (period, epoch, duration), using only points near the transit
    B = np.ptp(t); sig = 1.4826 * np.median(np.abs(f - np.median(f))); base = np.median(f)
    ph = ((t - t0) / P + .5) % 1 - .5
    sel = np.abs(ph) * P < 0.3
    ts, fs = t[sel], f[sel]; N = len(f); best = (-1e9, P, t0, dur)
    for p in P * (1 + np.linspace(-1, 1, 41) * dur / (3 * B)):
        ph0 = (((ts - t0) / p + .5) % 1 - .5) * p
        for d in g.WINS:
            for dt in np.linspace(-d / 2, d / 2, 7):
                a = np.abs(ph0 - dt) < d / 2; na = a.sum()
                if na < 3: continue
                dd = base - fs[a].mean(); s = dd / (sig * np.sqrt(1 / na + 1 / (N - na)))
                if s > best[0]: best = (s, p, t0 + dt, d)
    return best


def run_star(t, f, bad, nsig=3):
    t, f, nbad, niso = clean(t, f, bad)
    res = dict(n_badtime=nbad, n_iso=niso, npts=int(len(t)))
    if len(t) > 100:
        f, rem = prewhiten(t, f); res['prewhitened'] = [[float(p), float(a)] for p, a in rem]
        t, f, _, n2 = clean(t, f, np.array([], np.int64)); res['n_iso'] += n2
    if len(t) < 1000 or np.ptp(t) < 10:
        res['status'] = 'too_few'; return res, t, f
    B = np.ptp(t); periods = g.period_grid(B); sig = 1.4826 * np.median(np.abs(f - np.median(f)))
    res['mad_ppm'] = sig * 1e6; res['baseline'] = B
    m = np.ones(len(t), bool); sigs = []
    for k in range(nsig):
        r = g.bls_power(t[m], f[m], periods)
        p = r['snr']; med = np.median(p); mad = 1.4826 * np.median(np.abs(p - med)); j = int(np.argmax(p))
        snr, P, t0, dur = refine(t[m], f[m], periods[j], r['t0'][j] + t[m].min(), r['dur'][j])
        s = dict(P=float(P), t0=float(t0), dur=float(dur), snr=float(snr), sde=float((p[j] - med) / mad))
        s.update(stats(t[m], f[m], P, t0, dur, sig)); sigs.append(s)
        ph = ((t - t0) / P + .5) % 1 - .5; m &= ~(np.abs(ph) < 0.75 * dur / P)
    res['signals'] = sigs; res['status'] = 'ok'
    return res, t, f
