"""Deep pass 2 target list, in priority order.
P0: stars with an unexplained first-pass signal (SNR>=6, SDE>=6) and sectors the first pass did not use.
P1: total sectors >= 3x used (and >= 4).  P2: >= 2x used (and >= 4).  Within a tier: larger gain and smaller star first."""
import pandas as pd, numpy as np
D = '/home/sarang/tess_search'
r = pd.read_pickle(f'{D}/deep2/avail_all.pkl')
o = pd.read_csv(f'{D}/out3_signals.csv', low_memory=False)
lead = o[(o.snr >= 6) & (o.snr < 15) & (o.sde >= 6) & o.known_match.isna() & ~o.eb & (o.dep < 30000)]
lead_best = lead.groupby('tic').snr.max()
r['lead_snr'] = r.ID.map(lead_best).fillna(0)
cvz = set(pd.read_csv(f'{D}/data/targets.csv').ID.astype(int))            # first-pass CVZ stars were searched with one year of 2-min data
r.loc[r.ID.isin(cvz) & (r.nused == 0), 'nused'] = r.n2.clip(upper=13)
r['gain'] = np.sqrt(r.ntot / r.nused.clip(lower=1))
r['tier'] = 9
r.loc[(r.ntot >= 2 * r.nused.clip(lower=1)) & (r.ntot >= 4), 'tier'] = 2
r.loc[(r.ntot >= 3 * r.nused.clip(lower=1)) & (r.ntot >= 4), 'tier'] = 1
r.loc[(r.lead_snr > 0) & (r.ntot >= r.nused + 2) & (r.gain >= 1.25), 'tier'] = 0
r = r[r.tier < 9].copy()
r['score'] = np.where(r.tier == 0, r.lead_snr * r.gain, r.gain / r.rad.fillna(0.5).clip(lower=0.1))
r['heavy'] = (r.ntot > 24).astype(int)                                  # very long light curves go last within a tier
r = r.sort_values(['tier', 'heavy', 'score'], ascending=[True, True, False])
r.to_csv(f'{D}/deep2/targets.csv', index=False)
print(r.groupby('tier').agg(n=('ID', 'size'), files=('ntot', 'sum'), med_ntot=('ntot', 'median'), med_T=('Tmag', 'median')))
print(r.head(5)[['ID', 'tier', 'nused', 'ntot', 'Tmag', 'rad', 'lead_snr', 'score']])
