"""Noise-calibrated vetter for the deep pass.
Train: gradient-boosted trees on search statistics of injected transits (positives) vs peaks found in block-shuffled real light curves (negatives).
Stars are split into train/test by TIC (no star in both). Reports recall of injections at fixed false-alarm rates per SNR bin, then scores the
real deep-pass signals (deep2_signals.csv -> deep2_scored.csv).  Usage: deep2_vetml.py [train|score|both]"""
import sys, glob, json, os, numpy as np, pandas as pd
D = '/home/sarang/tess_search'
tg = pd.concat([pd.read_csv(f'{D}/data/{n}.csv') for n in ('targets', 'targets_allsky', 'targets_ffi_faint')]).drop_duplicates('ID').set_index('ID')
FEATS = ['snr', 'sde', 'snr_sde', 'logP', 'dur_h', 'dur_ratio', 'depth_sig', 'ntr', 'top2', 'medr', 'oe', 'sec', 'nsec_in', 'pos_frac', 'h_min', 'h_min_rel', 'h_asym',
         'chi2_red', 'snr_stat', 'coarse_gap', 'lognpts', 'nsec', 'logB', 'k']


def feats(g, k, npts, nsec, baseline, mad_ppm, rad, ms):
    P, dur, snr = g['P'], g['dur'], g['snr']; aR = 4.2084 * P ** (2 / 3) * ms ** (1 / 3) / rad; tmax = P / np.pi * np.arcsin(min(1.0, 1.0 / aR))
    oe = abs(g['odd'] - g['even']) / np.hypot(g['odd_e'], g['even_e']) if np.isfinite(g.get('odd_e', np.nan)) and np.isfinite(g.get('even_e', np.nan)) else 0.0
    sec = g['half'] / g['half_e'] if np.isfinite(g.get('half_e', np.nan)) and g.get('half_e', 0) > 0 else 0.0
    nin = g.get('nsec_in', 0) or 0; h1, h2 = g.get('snr_h1', 0.0), g.get('snr_h2', 0.0)
    return dict(snr=snr, sde=g['sde'], snr_sde=snr / max(g['sde'], 0.1), logP=np.log10(P), dur_h=dur * 24, dur_ratio=dur / tmax, depth_sig=g['depth_ppm'] / mad_ppm, ntr=g['ntr'],
                top2=g['top2_frac'], medr=g['med_ratio'], oe=oe, sec=sec, nsec_in=nin, pos_frac=(g.get('nsec_pos', 0) / nin) if nin else np.nan,
                h_min=min(h1, h2) if nin >= 2 else np.nan, h_min_rel=(min(h1, h2) / snr) if nin >= 2 else np.nan, h_asym=(abs(h1 - h2) / snr) if nin >= 2 else np.nan,
                chi2_red=(g.get('sec_chi2', np.nan) / (nin - 1)) if nin >= 2 else np.nan, snr_stat=g['depth_ppm'] / g['depth_e'] if g.get('depth_e', 0) > 0 else np.nan,
                coarse_gap=snr - g.get('snr_coarse', snr), lognpts=np.log10(npts), nsec=nsec, logB=np.log10(baseline), k=k)


def load_mlset():
    rows = []
    for fn in glob.glob(f'{D}/deep2/mlset/*.json'):
        r = json.load(open(fn)); inj = r['inj']
        for k, g in enumerate(r['signals']):
            if not all(np.isfinite(g.get(x, np.nan)) for x in ('depth_ppm', 'snr', 'P')): continue
            lab = 0
            if inj:
                rt = g['P'] / inj['P']
                if abs(rt - 1) < 0.003: lab = 1
                elif any(abs(rt / q - 1) < 0.01 for q in (0.5, 2, 1 / 3, 3, 1.5, 2 / 3, 0.25, 4)): lab = -1           # harmonic of the injection: exclude
            rows.append(dict(tic=r['tic'], label=lab, inj_snr=inj['snr_target'] if inj else np.nan, **feats(g, k, r['npts'], r['nsec'], r['baseline'], r['mad_ppm'], r['rad'], r['mass'])))
    return pd.DataFrame(rows)


def train():
    import pickle
    from sklearn.ensemble import HistGradientBoostingClassifier as HGB
    df = load_mlset(); df = df[df.label >= 0].reset_index(drop=True); test = (df.tic % 5 == 0)
    print(f'{df.tic.nunique()} stars, {len(df)} signals: {int((df.label == 1).sum())} injected-recovered, {int((df.label == 0).sum())} noise peaks; test stars {df[test].tic.nunique()}')
    P = dict(max_iter=300, learning_rate=0.05, max_leaf_nodes=15, min_samples_leaf=20, l2_regularization=2.0, random_state=0)
    m = HGB(**P).fit(df.loc[~test, FEATS], df.loc[~test, 'label']); s = m.predict_proba(df.loc[test, FEATS])[:, 1]; te = df[test].assign(score=s)
    from sklearn.metrics import roc_auc_score
    print(f'held-out stars: AUC {roc_auc_score(te.label, te.score):.4f};  SNR-only AUC {roc_auc_score(te.label, te.snr):.4f}')
    neg = te[te.label == 0]; pos = te[te.label == 1]; nst = te.tic.nunique()
    for far in (0.01, 0.003):                        # false alarms per star searched
        nfa = max(int(round(far * nst)), 1); thr_s = np.sort(neg.score.values)[-nfa]; thr_snr = np.sort(neg.snr.values)[-nfa]
        print(f'\nat ~{far*100:.1f} false alarms per 100 stars (test: {nfa} of {len(neg)} noise peaks pass): score > {thr_s:.3f}  |  equivalent pure-SNR cut: SNR > {thr_snr:.2f}')
        for lo, hi in ((5, 6), (6, 6.5), (6.5, 7), (7, 7.5), (7.5, 8.5), (8.5, 12)):
            b = pos[(pos.snr >= lo) & (pos.snr < hi)]
            if len(b): print(f'   injections recovered at SNR {lo}-{hi}: n={len(b):4d}  kept by vetter {np.mean(b.score > thr_s)*100:5.1f}%   kept by SNR cut {np.mean(b.snr > thr_snr)*100:5.1f}%')
        print(f'   all recovered injections: vetter {np.mean(pos.score > thr_s)*100:.1f}%  vs SNR cut {np.mean(pos.snr > thr_snr)*100:.1f}%')
    from sklearn.inspection import permutation_importance
    pi = permutation_importance(m, df.loc[test, FEATS], df.loc[test, 'label'], n_repeats=5, random_state=0, scoring='roc_auc'); imp = pd.Series(pi.importances_mean, FEATS).sort_values(ascending=False)
    print('\ntop features (permutation, AUC drop):', ', '.join(f'{k} {v:.3f}' for k, v in imp.head(8).items()))
    # final models: 5-fold by star, used as an ensemble; out-of-fold scores of the noise peaks calibrate the false-alarm rate
    models, oof = [], np.zeros(len(df))
    for kf in range(5):
        tr = df.tic % 5 != kf; mk = HGB(**P).fit(df.loc[tr, FEATS], df.loc[tr, 'label']); models.append(mk); oof[~tr.values] = mk.predict_proba(df.loc[~tr, FEATS])[:, 1]
    pickle.dump(dict(models=models, feats=FEATS, null_scores=np.sort(oof[df.label.values == 0]), n_null_stars=int(df.tic.nunique())), open(f'{D}/deep2/vetml.pkl', 'wb'))


def score():
    import pickle
    M = pickle.load(open(f'{D}/deep2/vetml.pkl', 'rb')); rows = []
    for fn in glob.glob(f'{D}/deep2/lc/*.res.json'):
        r = json.load(open(fn))
        if r.get('status') != 'ok' or int(r['tic']) not in tg.index: continue
        st = tg.loc[int(r['tic'])]; rad = float(st.rad) if np.isfinite(st.rad) else 0.4; ms = float(st['mass']) if np.isfinite(st['mass']) else 0.4
        for k, g in enumerate(r['signals']):
            if all(np.isfinite(g.get(x, np.nan)) for x in ('depth_ppm', 'snr', 'P')): rows.append(dict(tic=int(r['tic']), **feats(g, k, r['npts'], r['nsec'], r['baseline'], r['mad_ppm'], rad, ms)))
    f = pd.DataFrame(rows); f['vetml'] = np.mean([m.predict_proba(f[M['feats']])[:, 1] for m in M['models']], 0)
    # false alarms per 100 stars expected above this score, from the out-of-fold scrambled-null scores
    f['fa_per100'] = [(np.sum(M['null_scores'] >= v) + 0.5) / M['n_null_stars'] * 100 for v in f.vetml]
    o = pd.read_csv(f'{D}/deep2_signals.csv', low_memory=False).merge(f[['tic', 'k', 'vetml', 'fa_per100']], on=['tic', 'k'], how='left'); o.to_csv(f'{D}/deep2_scored.csv', index=False)
    hard = {'odd-even', 'secondary', 'deep', 'Rp>15', 'P-shared', 'syst-period', 'too-long', 'rotation-harmonic'}
    c = o[(o.vetml >= 0.5) & o.known_match.isna() & ~o.eb & ~o['flags'].fillna('').apply(lambda x: bool(set(x.split(';')) & hard))].sort_values('vetml', ascending=False)
    cols = ['tic', 'k', 'P', 'dur_h', 'dep', 'snr', 'sde', 'ntr', 'nsec_pos', 'nsec_in', 'h1', 'h2', 'Rp', 'Tmag', 'pass1_snr', 'vetml', 'fa_per100', 'flags', 'tce_match', 'public', 'known_other']
    print(f'{o.tic.nunique()} stars scored; unknown signals with score >= 0.5 and no hard flag: {len(c)}  (of which SNR < 7.5: {(c.snr < 7.5).sum()})')
    print(c[cols].head(80).to_string(index=False))


if __name__ == '__main__':
    mode = sys.argv[1] if len(sys.argv) > 1 else 'both'
    if mode in ('train', 'both'): train()
    if mode in ('score', 'both'): score()
