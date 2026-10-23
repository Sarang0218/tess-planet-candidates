import numpy as np, torch, torch.nn as nn

class Net(nn.Module):
    def __init__(s):
        super().__init__()
        def blk(i, o, n): return [nn.Conv1d(i, o, 5, padding=2), nn.GELU(), nn.Conv1d(o, o, 5, padding=2), nn.GELU(), nn.MaxPool1d(n, n)]
        s.g = nn.Sequential(*blk(1, 16, 2), *blk(16, 32, 2), *blk(32, 64, 2), *blk(64, 128, 2), nn.AdaptiveAvgPool1d(4), nn.Flatten())
        s.l = nn.Sequential(*blk(4, 32, 2), *blk(32, 64, 2), *blk(64, 128, 2), nn.AdaptiveAvgPool1d(4), nn.Flatten())
        s.h = nn.Sequential(nn.Linear(512 + 512 + 3, 256), nn.GELU(), nn.Dropout(.3), nn.Linear(256, 128), nn.GELU(), nn.Dropout(.3), nn.Linear(128, 1))
    def forward(s, g, l, sc):
        return s.h(torch.cat([s.g(g[:, None]), s.l(l), sc], 1))[:, 0]

def auc(p, y):
    o = np.argsort(p); r = np.empty(len(p)); r[o] = np.arange(len(p)); n1 = y.sum(); n0 = len(y) - n1
    return (r[y == 1].sum() - n1 * (n1 - 1) / 2) / (n1 * n0)
