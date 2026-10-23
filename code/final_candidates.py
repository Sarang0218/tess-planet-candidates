import pandas as pd, numpy as np, sys
sys.path.insert(0, '/home/sarang/tess_search/scripts'); import ffi_tce
o = pd.read_csv('/home/sarang/tess_search/out3_signals.csv', low_memory=False); v = pd.read_csv('/home/sarang/tess_search/data/vet3_scores.csv')
m = o.merge(v[['tic', 'k', 'cnn', 'vet3']], on=['tic', 'k'], how='left')
soft = {'few-epochs', 'ntr<3', 'repeat-P-in-star'}
def ok_flags(f): return (not isinstance(f, str)) or set(f.split(';')) <= soft
c = m[(m.snr >= 7) & (m.sde >= 6.5) & m.known_match.isna() & m.tce_match.isna() & m.public.isna() & ~m.eb & m['flags'].apply(ok_flags) & (m.vet3 >= 0.5)].copy()
c['ffi_tce'] = ''
for i, r in c.iterrows():
    ps = ffi_tce.periods(int(r.tic)); hit = [p for _, _, p in ps if any(abs(r.P / (p * mm) - 1) < .01 for mm in (.5, 1, 2))]
    c.at[i, 'ffi_tce'] = 'MATCH ' + ';'.join(f'{p:.4f}' for p in hit) if hit else ('host(other P)' if ps else '')
cols = ['tic', 'P', 'dur_h', 'dep', 'snr', 'sde', 'ntr', 'oe', 'top2', 'medr', 'Rp', 'Rs', 'Tmag', 'contam', 'cnn', 'vet3', 'flags', 'ffi_tce']
c = c.sort_values('vet3', ascending=False); c.to_csv('/home/sarang/tess_search/final_candidates.csv', index=False)
print('candidates:', len(c), '| FFI-TCE matches:', c.ffi_tce.str.startswith('MATCH').sum())
print(c[cols].round(3).to_string(index=False))
