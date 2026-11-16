"""Automatic stage-1 vetting of every unknown strong periodic signal of the deep pass -> deep2_shortlist.csv.
Cuts: duration physically possible for the star, most sectors positive and mutually consistent, no odd/even difference, no secondary,
no other dips or strong variability at other phases, most individual transits positive."""
import sys, io, contextlib, numpy as np, pandas as pd
sys.path.insert(0, '/home/sarang/tess_search/scripts'); import deep2_vet as dv
import os
R = '/home/sarang/tess_search'; s = pd.read_csv(f'{R}/deep2_scored.csv' if os.path.exists(f'{R}/deep2_scored.csv') else f'{R}/deep2_signals.csv', low_memory=False)
if 'vetml' not in s: s['vetml'] = np.nan
hard = {'repeat-P-in-star', 'odd-even', 'secondary', 'deep', 'Rp>15', 'P-shared', 'syst-period', 'too-long', 'rotation-harmonic', 'ntr<3'}
fl = s['flags'].fillna('').apply(lambda x: set(x.split(';')))
# strong signals, plus weaker ones (SNR >= 6.5) that the noise-calibrated classifier ranks above its 1-false-alarm-per-100-stars level
c = s[(((s.snr >= 7.5) & (s.sde >= 7)) | ((s.vetml >= 0.85) & (s.snr >= 6.5))) & (s.ntr >= 3) & (s.top2 < 0.85) & s.known_match.isna() & ~s.eb.astype(bool) & ~fl.apply(lambda x: bool(x & hard))]
print(f'{s.tic.nunique()} stars; unknown strong periodic signals without a hard flag: {len(c)}'); rows = []
for _, g in c.iterrows():
    try:
        with contextlib.redirect_stdout(io.StringIO()): o = dv.vet(int(g.tic), int(g.k), s)
    except Exception as ex: print('fail', g.tic, ex); continue
    why = []
    if o['T14_h'] > 1.35 * o['T14_max_h']: why.append('too-long')
    if o['nsec_pos'] < 0.75 * o['nsec']: why.append('sectors-negative')
    if o['nsec'] > 1 and o['chi2_sec'] / (o['nsec'] - 1) > 4: why.append('sectors-inconsistent')
    if o['oe_sig'] > 3: why.append('odd-even')
    if o['sec_sig'] > 3: why.append('secondary')
    if o['min_other_phase_sig'] < -4.5 or o['scan_std'] > 1.7: why.append('other-phase-structure')
    if o['frac_tr_pos'] < 0.6: why.append('transits-not-consistent')
    if o['snr_box'] < (7 if g.snr >= 7.5 else 5.5): why.append('weak-on-refit')
    if o['ingress_frac'] > 0.42 and o['Rp_Re'] > 4: why.append('V-shaped-large')
    rows.append(dict(tic=int(g.tic), k=int(g.k), P=g.P, T0=o['T0'], T14_h=o['T14_h'], T14_max_h=o['T14_max_h'], depth_ppm=o['trap_depth_ppm'], snr=g.snr, sde=g.sde, vetml=g.vetml, snr_box=o['snr_box'], Rp=o['Rp_Re'], ntr=o['ntr'],
                     nsec=o['nsec'], nsec_pos=o['nsec_pos'], Tmag=o['Tmag'], contam=o['contam'], pass1_snr=g.pass1_snr, tce=g.tce_match, public=g.public, known_other=g.known_other, flags=g['flags'], reject=';'.join(why)))
o = pd.DataFrame(rows).sort_values('snr', ascending=False); o.to_csv(f'{R}/deep2_shortlist.csv', index=False)
ok = o[o.reject == '']; pd.set_option('display.width', 260); pd.set_option('display.max_colwidth', 40)
print(f'pass all stage-1 cuts: {len(ok)} of {len(o)};  rejected for: ' + ', '.join(f'{k} {v}' for k, v in o.reject[o.reject != ''].str.split(';').explode().value_counts().items()))
print(ok.round(3).to_string(index=False))
