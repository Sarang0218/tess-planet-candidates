"""Batched box-least-squares on GPU (torch). Fixed phase-bin width ~ dur_min/3; cyclic sliding windows."""
import numpy as np, torch
BINW = 0.02 / 3
WINS = np.array([0.02, 0.03, 0.045, 0.065, 0.09, 0.13])

def period_grid(B, pmin=0.5, pmax=None):
    pmax = pmax or min(40.0, B / 2.0)
    P, ps = pmin, []
    while P < pmax:
        ps.append(P); P += P * 0.03 * P ** (1 / 3) / (3.0 * B)
    return np.array(ps)

@torch.compile(dynamic=True)
def _fold(t, inv, nbf, off):
    x = t[None, :] * inv[:, None]
    fr = x - torch.floor(x)
    b = torch.clamp(torch.minimum(torch.floor(fr * nbf[:, None]), nbf[:, None] - 1), min=0)
    return (b + off[:, None]).long().ravel()

def bls_power(t, f, periods, dev='cuda', chunk=1024):
    """Returns per-period best SNR, depth, t0, duration (constant per-point noise)."""
    t = torch.tensor(t - t.min(), dtype=torch.float32, device=dev)
    y = torch.tensor(f - np.median(f), dtype=torch.float32, device=dev)
    sig = float(1.4826 * np.median(np.abs(f - np.median(f))))
    N = len(f); ytot = y.sum(); ones = torch.ones(N * chunk, device=dev)
    wb = np.maximum(np.round(WINS / BINW).astype(int), 1)
    out = {k: np.zeros(len(periods)) for k in ('snr', 'depth', 't0', 'dur')}
    for c in range(0, len(periods), chunk):
        Pc = periods[c:c + chunk]; nb = np.maximum(np.round(Pc / BINW).astype(int), 8); NB = int(nb.max())
        Pt = torch.tensor(Pc, dtype=torch.float32, device=dev); nbt = torch.tensor(nb, device=dev)
        C = len(Pc)
        flat = _fold(t, (1 / Pt).float(), nbt.float(), (torch.arange(C, device=dev) * NB).float())
        S = torch.zeros(C * NB, device=dev).scatter_add_(0, flat, y.repeat(C))
        Nn = torch.zeros(C * NB, device=dev).scatter_add_(0, flat, ones[:N * C])
        S, Nn = S.view(C, NB), Nn.view(C, NB)
        # cyclic wrap: append the first max(w) bins of each period's own cycle after its last valid bin
        W = int(wb.max()); ar = torch.arange(NB + W, device=dev)[None, :]
        src = torch.remainder(ar, nbt[:, None])
        S2, N2 = torch.gather(S, 1, src), torch.gather(Nn, 1, src)
        cS = torch.nn.functional.pad(S2.cumsum(1), (1, 0)); cN = torch.nn.functional.pad(N2.cumsum(1), (1, 0))
        best = torch.full((len(Pc),), -1e9, device=dev); bd = torch.zeros_like(best); bt = torch.zeros_like(best); bw = torch.zeros_like(best)
        valid = ar[:, :NB] < nbt[:, None]
        for w in np.unique(wb):
            sin = cS[:, w:w + NB] - cS[:, :NB]; nin = cN[:, w:w + NB] - cN[:, :NB]
            nout = N - nin
            dep = -(sin / nin.clamp(min=1) - (ytot - sin) / nout.clamp(min=1))
            snr = dep / (sig * torch.sqrt(1 / nin.clamp(min=1) + 1 / nout.clamp(min=1)))
            snr = torch.where(valid & (nin >= 3), snr, torch.full_like(snr, -1e9))
            m, j = snr.max(1)
            upd = m > best
            best = torch.where(upd, m, best); bd = torch.where(upd, dep.gather(1, j[:, None])[:, 0], bd)
            bt = torch.where(upd, (j.double() + w / 2).float(), bt); bw = torch.where(upd, torch.full_like(bw, float(w)), bw)
        out['snr'][c:c + chunk] = best.cpu().numpy(); out['depth'][c:c + chunk] = bd.cpu().numpy()
        out['t0'][c:c + chunk] = (bt.cpu().numpy() / nb) * Pc + float(t.min()) * 0
        out['dur'][c:c + chunk] = bw.cpu().numpy() * BINW
    return out
