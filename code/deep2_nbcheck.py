"""Is a known TOI / CTOI / published candidate sitting next to a lead?  Positional match (3 arcmin) against the TOI and CTOI lists and
ID match of every TIC neighbour within 120" (deep2/pix/{TIC}/tic_120arcsec.csv) against all stored public lists.  Usage: deep2_nbcheck.py TIC ..."""
import sys, pandas as pd, numpy as np, os
from astropy.coordinates import SkyCoord
import astropy.units as u
R = '/home/sarang/tess_search'
toi = pd.read_csv(f'{R}/data/toi_20261116.csv', low_memory=False); ctoi = pd.read_csv(f'{R}/data/ctoi_20261116.csv', low_memory=False)
tc = SkyCoord(toi['RA'].astype(str), toi['Dec'].astype(str), unit=(u.hourangle, u.deg))
cra = [c for c in ctoi.columns if c.startswith('RA')][0]; cde = [c for c in ctoi.columns if c.startswith('Dec')][0]
try: cc = SkyCoord(pd.to_numeric(ctoi[cra]), pd.to_numeric(ctoi[cde]), unit='deg')
except Exception: cc = SkyCoord(ctoi[cra].astype(str), ctoi[cde].astype(str), unit=(u.hourangle, u.deg))
ids = {}
for f, c in (('lastovka2026.csv', 'tic'), ('t16_roth2026_candidates.csv', 'TICID'), ('t16_roth2026_single.csv', 'TICID'), ('diamante_v2.csv', 'ticID'), ('public_candidates.csv', 'tic')):
    ids[f] = set(pd.to_numeric(pd.read_csv(f'{R}/data/public/{f}', low_memory=False)[c], errors='coerce').dropna().astype(int))
ids['TOI'] = set(toi['TIC ID']); ids['CTOI'] = set(ctoi['TIC ID'])
tg = pd.concat([pd.read_csv(f'{R}/data/{n}.csv') for n in ('targets', 'targets_allsky', 'targets_ffi_faint')]).drop_duplicates('ID').set_index('ID')
print('TOI list RA/Dec parsed:', tc[0].to_string('hmsdms'), '| CTOI:', cc[0].to_string('hmsdms'), '| columns', cra, cde)
for t in [int(x) for x in sys.argv[1:]]:
    c0 = SkyCoord(float(tg.loc[t].ra), float(tg.loc[t].dec), unit='deg'); st = c0.separation(tc).arcsec; sc = c0.separation(cc).arcsec
    msg = [f"TOI {toi['TOI'].iloc[i]} (TIC {toi['TIC ID'].iloc[i]}) at {st[i]:.0f}\"" for i in np.where(st < 180)[0]] + [f"CTOI TIC {ctoi['TIC ID'].iloc[i]} at {sc[i]:.0f}\"" for i in np.where(sc < 180)[0]]
    p = f'{R}/deep2/pix/{t}/tic_120arcsec.csv'; nb = 'no neighbour file'
    if os.path.exists(p):
        nbid = set(pd.read_csv(p).ID.astype(int)) - {t}; h = [f'{k}:{sorted(nbid & v)}' for k, v in ids.items() if nbid & v]; nb = f'{len(nbid)} neighbours; in lists: {h if h else "none"}'
    print(t, '| within 3 arcmin:', msg if msg else 'none', '|', nb)
