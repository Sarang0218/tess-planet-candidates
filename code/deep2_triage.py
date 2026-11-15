"""Triage for the deep (all-sector) pass: deep2/lc/*.res.json -> deep2_signals.csv. Same flags as triage3 plus multi-sector
consistency flags; cross-match with TOI/CTOI (fresh), confirmed planets, SPOC TCEs, EB catalog, public community lists, and the first pass."""
import glob, json, os, numpy as np, pandas as pd
D = '/home/sarang/tess_search/data'; R = '/home/sarang/tess_search'
tg = pd.concat([pd.read_csv(f'{D}/targets.csv'), pd.read_csv(f'{D}/targets_allsky.csv'), pd.read_csv(f'{D}/targets_ffi_faint.csv')]).drop_duplicates('ID').set_index('ID')
toi = pd.read_csv(f'{D}/toi_20261116.csv'); ctoi = pd.read_csv(f'{D}/ctoi_20261116.csv')
ps = pd.read_csv(f'{D}/ps.csv'); ps['tic'] = ps.tic_id.str.replace('TIC ', '').astype(float)
eb = pd.read_csv(f'{D}/tess_eb.tsv', sep='\t', comment='#', dtype=str); eb = eb[eb.TIC.str.strip().str.isdigit()]; eb_tic = set(eb.TIC.astype(int))
known = {}
for df, col, pcol, lab in [(toi, 'TIC ID', 'Period (days)', 'TOI'), (ctoi, 'TIC ID', 'Period (days)', 'CTOI'), (ps, 'tic', 'pl_orbper', 'PS')]:
    for tic, P in zip(df[col], df[pcol]):
        P = pd.to_numeric(P, errors='coerce')
        if np.isfinite(tic): known.setdefault(int(tic), []).append((lab, float(P) if np.isfinite(P) else np.nan))
tce = pd.read_csv(f'{D}/tce_all.csv'); tce_k = {}
for tic, P in zip(tce.tic, tce.P):
    if np.isfinite(P) and P > 0: tce_k.setdefault(int(tic), []).append(float(P))
pub = pd.read_csv(f'{D}/public/public_candidates.csv'); pub_k = {}
for tic, P, src in zip(pub.tic, pub.P, pub.source): pub_k.setdefault(int(tic), []).append((src, float(P)))
if os.path.exists(f'{D}/public/lastovka2026.csv'):
    for tic, P in pd.read_csv(f'{D}/public/lastovka2026.csv')[['tic', 'P']].values: pub_k.setdefault(int(tic), []).append(('Lastovka+2026', float(P)))
o3 = pd.read_csv(f'{R}/out3_signals.csv', low_memory=False); p1 = {t: list(zip(g.P, g.snr)) for t, g in o3[o3.snr >= 5].groupby('tic')}
RATIOS = sorted({p / q for p in range(1, 5) for q in range(1, 5)})
rows = []
for fn in glob.glob(f'{R}/deep2/lc/*.res.json'):
    r = json.load(open(fn))
    if r.get('status') != 'ok': continue
    prot = sorted({round(p, 4) for p, a in r.get('prewhitened', [])})
    for k, s in enumerate(r['signals']): rows.append(dict(tic=int(r['tic']), k=k, mad=r['mad_ppm'], prot=prot, nsec=r['nsec'], baseline=r['baseline'], **s))
df = pd.DataFrame(rows); df = df[np.isfinite(df.depth_ppm) & np.isfinite(df.snr) & np.isfinite(df.P)].reset_index(drop=True)
sP = df[df.sde > 7]
def n_shared(P, t0, dur, tc):
    o = sP[(np.abs(sP.P / P - 1) < 0.002) & (sP.tic != tc)]
    if not len(o): return 0
    dphase = np.abs(((o.t0.values - t0) / P + .5) % 1 - .5) * P
    return o[dphase < 2 * max(dur, 0.02)].tic.nunique()
df['n_sameP'] = [n_shared(P, t0, du, tc) for P, t0, du, tc in zip(df.P, df.t0, df.dur, df.tic)]
def ratio_match(P, Pk): return np.isfinite(Pk) and Pk > 0 and any(abs(P / (Pk * m) - 1) < 0.01 for m in (0.5, 1, 2, 1 / 3, 3))
out = []
for _, s in df.iterrows():
    st = tg.loc[s.tic] if s.tic in tg.index else None; rad = float(st.rad) if st is not None else np.nan; fl = []
    oe = abs(s.odd - s.even) / np.hypot(s.odd_e, s.even_e) if np.isfinite(s.odd_e) and np.isfinite(s.even_e) else 0
    sec = s.half / s.half_e if np.isfinite(s.half_e) and s.half_e > 0 else 0
    rp = np.sqrt(max(s.depth_ppm, 0) * 1e-6) * rad * 109.1
    if oe > 3: fl.append('odd-even')
    if sec > 3 and s.half > 0.1 * s.depth_ppm: fl.append('secondary')
    if s.depth_ppm > 30000: fl.append('deep')
    if rp > 15: fl.append('Rp>15')
    if s.n_sameP >= 3: fl.append('P-shared')
    if any(abs(s.P / x - 1) < 0.01 for x in (13.7, 27.4, 1.0, 0.5, 2.0)): fl.append('syst-period')
    if s.ntr < 3: fl.append('ntr<3')
    ms = float(st['mass']) if st is not None and np.isfinite(st['mass']) else 0.4
    if np.isfinite(rad) and rad > 0:
        aR = 4.2084 * s.P ** (2 / 3) * ms ** (1 / 3) / rad; tmax = s.P / np.pi * np.arcsin(min(1.0, 1.2 / aR))
        if s.dur > 2.0 * tmax + 0.01: fl.append('too-long')
    if s.top2_frac > 0.6: fl.append('few-epochs')
    if s.med_ratio < 0.4: fl.append('inconsistent-depths')
    if any(abs(s.P / (pr * q) - 1) < 0.01 for pr in s.prot for q in RATIOS): fl.append('rotation-harmonic')
    nin = s.get('nsec_in', np.nan)
    if np.isfinite(nin) and nin >= 4:
        if min(s.snr_h1, s.snr_h2) < 2.0: fl.append('one-half-only')
        if s.nsec_pos / nin < 0.7: fl.append('few-sectors-positive')
        if s.sec_chi2 / (nin - 1) > 3.0: fl.append('inconsistent-sectors')
    kn = known.get(s.tic, []); f1 = [(P, sn) for P, sn in p1.get(s.tic, []) if ratio_match(s.P, P)]
    out.append(dict(tic=s.tic, k=s.k, P=s.P, t0=s.t0, dur=s.dur, dur_h=round(s.dur * 24, 2), dep=round(s.depth_ppm), snr=round(s.snr, 1), sde=round(s.sde, 1), ntr=s.ntr,
                    oe=round(oe, 1), sec=round(sec, 1), top2=round(s.top2_frac, 2), medr=round(s.med_ratio, 2), nsec=s.nsec, nsec_in=nin, nsec_pos=s.get('nsec_pos', np.nan),
                    h1=round(s.get('snr_h1', np.nan), 1), h2=round(s.get('snr_h2', np.nan), 1), Rp=round(rp, 2), Tmag=st.Tmag if st is not None else np.nan, Rs=rad,
                    Ms=st['mass'] if st is not None else np.nan, contam=st.contratio if st is not None else np.nan, d=st.d if st is not None else np.nan,
                    prot=';'.join(f'{p:.3f}' for p in s.prot[:4]), flags=fl, pass1_snr=max([sn for _, sn in f1], default=np.nan),
                    known_match=';'.join(f'{l}:{P:.3f}' for l, P in kn if ratio_match(s.P, P)),
                    public=';'.join(f'{src}:{P:.3f}' for src, P in pub_k.get(s.tic, []) if (not np.isfinite(P)) or ratio_match(s.P, P)),
                    tce_match=';'.join(sorted({f'{P:.3f}' for P in tce_k.get(s.tic, []) if ratio_match(s.P, P)})),
                    known_other=';'.join(f'{l}:{P:.3f}' for l, P in kn if not ratio_match(s.P, P)), eb=s.tic in eb_tic))
o = pd.DataFrame(out)
rep = o.groupby('tic').P.transform(lambda p: [int(np.sum(np.abs(p / x - 1) < .01)) - 1 for x in p])
o['flags'] = [';'.join(f + (['repeat-P-in-star'] if r >= 1 else [])) for f, r in zip(o['flags'], rep)]
o.to_csv(f'{R}/deep2_signals.csv', index=False)
T1, T2 = float(os.environ.get('SNR', 7.5)), float(os.environ.get('SDE', 7.0)); strong = o[(o.snr >= T1) & (o.sde >= T2)]
kh = {t for t in o.tic.unique() if t in known}
print(f'{o.tic.nunique()} stars searched, {len(o)} signals, {len(strong)} strong (SNR>={T1}, SDE>={T2}); known hosts: {len(kh)}, with a strong matching signal: {strong[strong.known_match != ""].tic.nunique()}')
cols = ['tic', 'k', 'P', 'dur_h', 'dep', 'snr', 'sde', 'ntr', 'oe', 'sec', 'top2', 'medr', 'nsec_pos', 'nsec_in', 'h1', 'h2', 'Rp', 'Tmag', 'contam', 'pass1_snr', 'prot', 'known_other', 'tce_match', 'public']
clean = strong[(strong.known_match == '') & (strong['flags'] == '') & ~strong.eb]
print('\nCLEAN, not TOI/CTOI/confirmed:', len(clean)); print(clean.sort_values('snr', ascending=False)[cols].to_string(index=False))
soft = strong[(strong.known_match == '') & (strong['flags'] != '') & ~strong.eb & strong['flags'].apply(lambda f: set(f.split(';')) <= {'few-epochs', 'ntr<3', 'repeat-P-in-star'})]
print('\nSOFT-flag only (few-epochs / ntr<3 / repeat-P):', len(soft)); print(soft.sort_values('snr', ascending=False)[cols + ['flags']].head(40).to_string(index=False))
rej = strong[(strong.known_match == '') & ((strong['flags'] != '') | strong.eb)]
print('\nflag counts among rejected:'); print(rej['flags'].str.split(';').explode().value_counts().to_string())
print('\nknown planets/TOIs recovered:'); print(strong[strong.known_match != ''].sort_values('snr', ascending=False)[['tic', 'P', 'dep', 'snr', 'sde', 'pass1_snr', 'known_match', 'flags']].head(60).to_string(index=False))
