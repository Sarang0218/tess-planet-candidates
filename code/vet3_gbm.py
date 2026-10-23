import numpy as np, pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.model_selection import GroupKFold
D = '/home/sarang/tess_search/data'
def auc(p, y):
    o = np.argsort(p); r = np.empty(len(p)); r[o] = np.arange(len(p)); n1 = y.sum(); n0 = len(y) - n1
    return (r[y == 1].sum() - n1 * (n1 - 1) / 2) / (n1 * n0)
dom = np.load(f'{D}/domain_set.npz'); pa = np.load(f'{D}/vet3_cnn_scores.npz')['pa']
X = np.column_stack([dom['SC'], dom['S'], pa]); yy = dom['y']; grp = dom['tic']; lab = yy >= 0
oof = np.full(len(yy), np.nan); final = np.zeros(len(yy))
for trn, tst in GroupKFold(5).split(X[lab], yy[lab], grp[lab]):
    ii, jj = np.flatnonzero(lab)[trn], np.flatnonzero(lab)[tst]
    clf = HistGradientBoostingClassifier(max_iter=300, learning_rate=.05, max_leaf_nodes=15, class_weight='balanced').fit(X[ii], yy[ii])
    oof[jj] = clf.predict_proba(X[jj])[:, 1]; final += clf.predict_proba(X)[:, 1] / 5
print(f'(B) diagnostics+CNN, star-grouped CV: AUC {auc(oof[lab], yy[lab]):.3f}; known planets >0.5: {np.mean(oof[lab & (yy == 1)] > .5):.3f}; FPs >0.5: {np.mean(oof[lab & (yy == 0)] > .5):.4f}')
sc = np.where(lab, oof, final)
pd.DataFrame(dict(tic=dom['tic'], P=dom['P'], k=dom['k'], y=yy, cnn=pa, vet3=sc)).to_csv(f'{D}/vet3_scores.csv', index=False)
