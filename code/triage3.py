"""Triage for the unified all-sky pipeline (out3/*.res.json). Catalog cross-match + vetting flags incl. rotation harmonics.
Writes out3_signals.csv; prints recovered knowns and clean candidates."""
import glob, json, numpy as np, pandas as pd
D = '/home/sarang/tess_search/data'
tg = pd.concat([pd.read_csv(f'{D}/targets.csv'), pd.read_csv(f'{D}/targets_allsky.csv'), pd.read_csv(f'{D}/targets_ffi_faint.csv')]).drop_duplicates('ID').set_index('ID')
toi = pd.read_csv(f'{D}/toi.csv'); ctoi = pd.read_csv(f'{D}/ctoi.csv')
ps = pd.read_csv(f'{D}/ps.csv'); ps['tic'] = ps.tic_id.str.replace('TIC ', '').astype(float)
eb = pd.read_csv(f'{D}/tess_eb.tsv', sep='\t', comment='#', dtype=str)
eb = eb[eb.TIC.str.strip().str.isdigit()]; eb_tic = set(eb.TIC.astype(int))
known = {}
for df, col, pcol, lab in [(toi, 'TIC ID', 'Period (days)', 'TOI'), (ctoi, 'TIC ID', 'Period (days)', 'CTOI'), (ps, 'tic', 'pl_orbper', 'PS')]:
    for tic, P in zip(df[col], df[pcol]):
        P = pd.to_numeric(P, errors='coerce')
        if np.isfinite(tic) and np.isfinite(P): known.setdefault(int(tic), []).append((lab, float(P)))
tce = pd.read_csv(f'{D}/tce_all.csv'); tce_k = {}
for tic, P in zip(tce.tic, tce.P):
    if np.isfinite(P) and P > 0: tce_k.setdefault(int(tic), []).append(float(P))
pub = pd.read_csv(f'{D}/public/public_candidates.csv'); pub_k = {}
for tic, P, src in zip(pub.tic, pub.P, pub.source): pub_k.setdefault(int(tic), []).append((src, float(P)))
RATIOS = sorted({p / q for p in range(1, 5) for q in range(1, 5)})

rows = []
for fn in glob.glob('/home/sarang/tess_search/out3/*.res.json'):
    r = json.load(open(fn))
    if r.get('status') != 'ok': continue
    prot = [p for p, a in r.get('prewhitened', [])]
    for k, s in enumerate(r['signals']):
        rows.append(dict(tic=int(r['tic']), k=k, mad=r['mad_ppm'], prot=prot, nsig=len(r['signals']), **s))
df = pd.DataFrame(rows)
df = df[np.isfinite(df.depth_ppm) & np.isfinite(df.snr) & np.isfinite(df.P)].reset_index(drop=True)
sP = df[df.sde > 7]
def n_shared(P, t0, dur, tc):
    # systematics share the period AND the transit times; chance period matches among ~1e5 signals do not
    o = sP[(np.abs(sP.P / P - 1) < 0.002) & (sP.tic != tc)]
    if not len(o): return 0
    dphase = np.abs(((o.t0.values - t0) / P + .5) % 1 - .5) * P
    return o[dphase < 2 * max(dur, 0.02)].tic.nunique()
df['n_sameP'] = [n_shared(P, t0, du, tc) for P, t0, du, tc in zip(df.P, df.t0, df.dur, df.tic)]

def ratio_match(P, Pk):
    return np.isfinite(Pk) and Pk > 0 and any(abs(P / (Pk * m) - 1) < 0.01 for m in (0.5, 1, 2, 1 / 3, 3))

import os
LC = {}
out = []
for _, s in df.iterrows():
    st = tg.loc[s.tic] if s.tic in tg.index else None
    rad = float(st.rad) if st is not None else np.nan; fl = []
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
        aR = 4.2084 * s.P ** (2 / 3) * ms ** (1 / 3) / rad
        tmax = s.P / np.pi * np.arcsin(min(1.0, 1.2 / aR))          # central transit, generous k
        if s.dur > 2.0 * tmax + 0.01: fl.append('too-long')
    if s.top2_frac > 0.6: fl.append('few-epochs')
    if s.med_ratio < 0.4: fl.append('inconsistent-depths')
    if any(abs(s.P / (pr * q) - 1) < 0.01 for pr in s.prot for q in RATIOS): fl.append('rotation-harmonic')
    if s.snr >= float(os.environ.get('THR', 7)) and s.sde >= float(os.environ.get('THR', 7)):
        if s.tic not in LC:
            dd = np.load(f'/home/sarang/tess_search/out3/{int(s.tic)}.npz'); LC.clear(); LC[s.tic] = (dd['t'], dd['f'])
        tt, ff = LC[s.tic]; ph = (((tt - s.t0) / s.P + .5) % 1 - .5); oot = np.abs(ph * s.P) > 1.5 * s.dur
        k = np.minimum((((ph[oot] + .5) * 40).astype(int)), 39); prof = np.array([np.median(ff[oot][k == j]) for j in range(40) if np.sum(k == j) > 20])
        if len(prof) > 10 and (np.percentile(prof, 95) - np.percentile(prof, 5)) / 2 > 0.3 * s.depth_ppm * 1e-6: fl.append('oot-variability')
    kn = known.get(s.tic, [])
    out.append(dict(tic=s.tic, k=s.k, P=s.P, t0=s.t0, dur=s.dur, dur_h=round(s.dur * 24, 2), dep=round(s.depth_ppm), snr=round(s.snr, 1),
                    sde=round(s.sde, 1), ntr=s.ntr, oe=round(oe, 1), sec=round(sec, 1), top2=round(s.top2_frac, 2), medr=round(s.med_ratio, 2),
                    Rp=round(rp, 2), Tmag=st.Tmag if st is not None else np.nan, Rs=rad, Ms=st['mass'] if st is not None else np.nan,
                    contam=st.contratio if st is not None else np.nan, prot=';'.join(f'{p:.3f}' for p in s.prot), flags=fl,
                    known_match=';'.join(f'{l}:{P:.3f}' for l, P in kn if ratio_match(s.P, P)),
                    public=';'.join(f'{src}:{P:.3f}' for src, P in pub_k.get(s.tic, []) if (not np.isfinite(P)) or ratio_match(s.P, P)),
                    tce_match=';'.join(sorted({f'{P:.3f}' for P in tce_k.get(s.tic, []) if ratio_match(s.P, P)})),
                    known_other=';'.join(f'{l}:{P:.3f}' for l, P in kn if not ratio_match(s.P, P)), eb=s.tic in eb_tic))
o = pd.DataFrame(out)
rep = o.groupby('tic').P.transform(lambda p: [int(np.sum(np.abs(p / x - 1) < .01)) - 1 for x in p])
o['flags'] = [';'.join(f + (['repeat-P-in-star'] if r >= 1 else [])) for f, r in zip(o['flags'], rep)]
o.to_csv('/home/sarang/tess_search/out3_signals.csv', index=False)
import os
THR = float(os.environ.get('THR', 7))
strong = o[(o.sde >= THR) & (o.snr >= THR)]
kn_hosts = {t for t in o.tic.unique() if t in known}
print(f'{o.tic.nunique()} stars searched, {len(o)} signals, {len(strong)} strong; known-planet hosts searched: {len(kn_hosts)}, '
      f'with a strong recovered signal: {strong[strong.known_match != ""].tic.nunique()}')
clean = strong[(strong.known_match == '') & (strong['flags'] == '') & ~strong.eb]
print('clean signals also in public hobbyist lists (ahrazzle/comdex4):', (clean.public != '').sum())
clean = clean[clean.public == '']
cols = ['tic', 'public', 'P', 'dur_h', 'dep', 'snr', 'sde', 'ntr', 'oe', 'top2', 'medr', 'Rp', 'Tmag', 'contam', 'prot', 'known_other', 'tce_match']
print('\nclean, not TOI/CTOI/confirmed:', len(clean))
print('  of which SPOC TCE at same period (detected, not promoted):', (clean.tce_match != '').sum())
print(clean[clean.tce_match != ''].sort_values('snr', ascending=False)[cols].to_string(index=False))
print('\n  TRULY NEW (no TOI/CTOI/confirmed/TCE):', (clean.tce_match == '').sum())
print(clean[clean.tce_match == ''].sort_values('snr', ascending=False)[cols].to_string(index=False))
rej = strong[(strong.known_match == '') & ((strong['flags'] != '') | strong.eb)]
print('\nflag counts among rejected:'); print(rej['flags'].str.split(';').explode().value_counts().to_string())
