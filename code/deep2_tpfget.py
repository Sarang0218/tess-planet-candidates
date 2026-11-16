"""Download target-pixel files (TESS-SPOC FFI stamps, or 2-min SPOC files where the star was a 2-min target) for every sector of a star
into deep2/pix/{TIC}/, newest first, and the TIC neighbour list if missing.   Usage: deep2_tpfget.py TIC [max_sectors]"""
import sys, os, json, socket, urllib.request, numpy as np, pandas as pd
socket.setdefaulttimeout(180)
R = '/home/sarang/tess_search'; TIC = int(sys.argv[1]); NMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 99
D = f'{R}/deep2/pix/{TIC}'; os.makedirs(D, exist_ok=True); t16 = f'{TIC:016d}'
if not os.path.exists(f'{D}/tic_120arcsec.csv'):
    from astroquery.mast import Catalogs
    tg = pd.concat([pd.read_csv(f'{R}/data/{n}.csv') for n in ('targets', 'targets_allsky', 'targets_ffi_faint')]).drop_duplicates('ID').set_index('ID').loc[TIC]
    Catalogs.query_region(f'{float(tg.ra)} {float(tg.dec)}', radius=120 / 3600, catalog='TIC').to_pandas().to_csv(f'{D}/tic_120arcsec.csv', index=False)
_pl = pd.read_csv(f'{R}/deep2/targets.csv', usecols=['ID', 'plan']).set_index('ID'); PLAN = json.loads(_pl.loc[TIC, 'plan']) if TIC in _pl.index else {}
src = json.load(open(f'{R}/deep2/lc/{TIC}.meta.json'))['src']; secs = sorted([int(s) for s in src], reverse=True); got = 0
for s in secs:
    if got >= NMAX: break
    n = f'hlsp_tess-spoc_tess_phot_{t16}-s{s:04d}_tess_v1_tp.fits'; url = f'https://archive.stsci.edu/hlsps/tess-spoc/s{s:04d}/target/{t16[0:4]}/{t16[4:8]}/{t16[8:12]}/{t16[12:16]}/{n}'
    if str(PLAN.get(str(s), '')).startswith('tm:'):
        n = PLAN[str(s)][3:].replace('_lc.fits', '_tp.fits'); x = n.split('-')[2]; url = f'https://stpubdata.s3.amazonaws.com/tess/public/tid/s{s:04d}/{x[0:4]}/{x[4:8]}/{x[8:12]}/{x[12:16]}/{n}'
    p = f'{D}/{n}'
    for _try in range(2):
        if os.path.exists(p) and os.path.getsize(p) > 100000: break
        try: urllib.request.urlretrieve(url, p + '.part'); os.replace(p + '.part', p)
        except Exception as ex: print('  tpf', s, ex, flush=True)
    got += os.path.exists(p) and os.path.getsize(p) > 100000
print(f'TIC {TIC}: {got} of {len(secs)} sectors have pixel files', flush=True)
