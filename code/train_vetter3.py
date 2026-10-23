"""Vetter v3: (A) ensemble of view-based CNNs trained on injections + real TOI/EB labels + our own labelled search signals;
(B) gradient boosting on our search diagnostics + CNN score (star-grouped CV). Scores all domain signals -> data/vet3_scores.csv"""
import sys, os, numpy as np, torch, torch.nn as nn, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cnn_model import Net, auc
dev = 'cuda'; D = '/home/sarang/tess_search/data'
def is_val(t): return (int(t) * 2654435761 % 2 ** 32) % 5 == 0
inj, real, dom = np.load(f'{D}/train_inj.npz'), np.load(f'{D}/real_labels.npz'), np.load(f'{D}/domain_set.npz')
def part(d, idx, wt): return [d['G'][idx], d['L'][idx], d['S'][idx], d['y'][idx].astype(np.float32), np.full(len(idx), wt, np.float32)]
rv = np.array([is_val(t) for t in real['tic']]); dv = np.array([is_val(t) for t in dom['tic']]); dl = dom['y'] >= 0
tr_parts = [part(inj, np.arange(len(inj['y'])), 1.0), part(real, np.flatnonzero(~rv), 6.0)]
di = np.flatnonzero(~dv & dl); dp = [dom['G'][di], dom['L'][di], dom['S'][di], dom['y'][di].astype(np.float32), dom['w'][di].astype(np.float32) * np.where(dom['y'][di] == 1, 40., 4.)]
tr_parts.append(dp)
G, L, S, y, w = [np.concatenate([p[i] for p in tr_parts]) for i in range(5)]
T = lambda a: torch.tensor(np.nan_to_num(a), dtype=torch.float32, device=dev)
Gt, Lt, St, yt, wt = T(G), T(L), T(S), T(y), T(w)
print('train examples', len(y), 'positives', int(y.sum()))
def predict(net, G_, L_, S_):
    net.eval(); out = []
    with torch.no_grad():
        for i in range(0, len(G_), 8192): out.append(torch.sigmoid(net(T(G_[i:i + 8192]), T(L_[i:i + 8192]), T(S_[i:i + 8192]))).cpu().numpy())
    return np.concatenate(out)
nets = []
for seed in range(5):
    torch.manual_seed(seed); np.random.seed(seed); net = Net().to(dev)
    opt = torch.optim.AdamW(net.parameters(), 1e-3, weight_decay=1e-2); E = 25; nb = len(y) // 512 + 1
    sch = torch.optim.lr_scheduler.OneCycleLR(opt, 2e-3, total_steps=E * nb)
    for ep in range(E):
        net.train(); perm = torch.randperm(len(y), device=dev)
        for i in range(0, len(y), 512):
            j = perm[i:i + 512]; g, l = Gt[j], Lt[j]
            if np.random.rand() < .5: g, l = g.flip(-1), l.flip(-1)
            loss = (nn.functional.binary_cross_entropy_with_logits(net(g, l, St[j]), yt[j], reduction='none') * wt[j]).sum() / wt[j].sum()
            opt.zero_grad(); loss.backward(); opt.step(); sch.step()
    nets.append(net); print('model', seed, 'done', flush=True)
def ens(G_, L_, S_): return np.mean([predict(n, G_, L_, S_) for n in nets], 0)
# --- evaluation on held-out stars
pr = ens(real['G'][rv], real['L'][rv], real['S'][rv]); yr = real['y'][rv]; kr = real['kind'][rv]; toi = np.char.startswith(kr.astype(str), 'toi')
print(f'\nREAL held-out: AUC(all) {auc(pr, yr):.3f}  AUC(TOI planet vs FP) {auc(pr[toi], yr[toi]):.3f}')
for k in np.unique(kr): print(f'   {k:8s} n={np.sum(kr == k):4d} pass(>0.5) {np.mean(pr[kr == k] > .5):.3f}')
pa = ens(dom['G'], dom['L'], dom['S']); m = dv & dl
print(f'DOMAIN held-out (our M-dwarf search signals): AUC {auc(pa[m], dom["y"][m]):.3f}; known planets pass {np.mean(pa[m & (dom["y"] == 1)] > .5):.3f}, '
      f'labelled FPs pass {np.mean(pa[m & (dom["y"] == 0)] > .5):.3f}  (n_pos={np.sum(m & (dom["y"] == 1))}, n_neg={np.sum(m & (dom["y"] == 0))})')
np.savez(f'{D}/vet3_cnn_scores.npz', pa=pa)
for i, n in enumerate(nets): torch.save(n.state_dict(), f'{D}/vet3_cnn_{i}.pt')
