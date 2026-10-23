"""Real-label vetting set: TOI CP/KP (1) vs FP/FA (0), plus TESS EBs (0). Views at catalog ephemerides."""
import sys, os, json, numpy as np, pandas as pd, urllib.request, warnings, multiprocessing as mp, shutil
warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from views import make_views
D = '/home/sarang/tess_search/data'
AV = pd.read_pickle(f'{D}/avail2min.pkl')
S3 = 'https://stpubdata.s3.amazonaws.com/tess/public/tid'

def targets():
    toi = pd.read_csv(f'{D}/toi.csv')
    toi = toi[toi['TFOPWG Disposition'].isin(['CP', 'KP', 'FP', 'FA']) & toi['TIC ID'].isin(AV)]
    toi = toi.dropna(subset=['Period (days)', 'Epoch (BJD)', 'Duration (hours)'])
    toi = toi[(toi['Period (days)'] > 0.2) & (toi['Period (days)'] < 100)]
    a = pd.DataFrame(dict(tic=toi['TIC ID'].astype(int), P=toi['Period (days)'], t0=toi['Epoch (BJD)'] - 2457000,
                          dur=toi['Duration (hours)'] / 24, y=toi['TFOPWG Disposition'].isin(['CP', 'KP']).astype(int),
                          kind='toi_' + toi['TFOPWG Disposition'], Ms=toi['Stellar Mass (M_Sun)'], Rs=toi['Stellar Radius (R_Sun)'],
                          name='TOI-' + toi['TOI'].astype(str)))
    eb = pd.read_csv(f'{D}/tess_eb_full.tsv', sep='\t', comment='#', skiprows=[], low_memory=False)
    eb = eb[pd.to_numeric(eb['TIC'], errors='coerce').notna()]
    for c in ['TIC', 'BJD0', 'Per', 'Wp-pf', 'Morph']: eb[c] = pd.to_numeric(eb[c], errors='coerce')
    eb = eb.dropna(subset=['TIC', 'BJD0', 'Per', 'Wp-pf'])
    eb = eb[eb.TIC.astype(int).isin(AV) & (eb.Per > 0.2) & (eb.Per < 100) & (eb['Wp-pf'] > 0) & (eb['Wp-pf'] < .3)]
    eb = eb.sample(frac=1, random_state=0)
    b = pd.DataFrame(dict(tic=eb.TIC.astype(int), P=eb.Per, t0=eb.BJD0, dur=eb['Wp-pf'] * eb.Per, y=0,
                          kind='tess_eb', Ms=np.nan, Rs=np.nan, name='EB'))
    return pd.concat([a, b], ignore_index=True)

RD = f'{D}/real'

def one(row):
    key = f'{RD}/{int(row.tic)}_{row.P:.6f}.npz'
    if os.path.exists(key) or os.path.exists(key[:-4] + '.none'): return 'skip'
    r = _one(row)
    if r is not None:
        G, L, Sc, y, kind, tic, name = r
        np.savez(key, G=G, L=L, S=Sc, y=y, kind=kind, tic=tic, name=name)
    else: open(key[:-4] + '.none', 'w').close()
    return r

def _one(row):
    tic, P, t0, dur, Ms, Rs = int(row.tic), row.P, row.t0, row.dur, row.Ms, row.Rs
    if not np.isfinite(Ms) or Ms <= 0: Ms = 1.0
    if not np.isfinite(Rs) or Rs <= 0: Rs = 1.0
    import lightkurve as lk
    files = sorted(AV[tic])[-4:]                      # latest <=4 sectors
    d = f'/tmp/rl_{tic}_{os.getpid()}'; os.makedirs(d, exist_ok=True); parts = []
    try:
        for s, fn in files:
            t16 = fn.split('-')[2]
            u = f'{S3}/s{s:04d}/{t16[0:4]}/{t16[4:8]}/{t16[8:12]}/{t16[12:16]}/{fn}'
            p = f'{d}/{fn}'
            try: urllib.request.urlretrieve(u, p); lc = lk.read(p, quality_bitmask='default')
            except Exception: continue
            lc = lc.remove_nans().normalize(); lc = lc[np.isfinite(lc.flux.value)]
            if len(lc) < 500: continue
            cad = np.nanmedian(np.diff(lc.time.value)); win = max(0.75, 3 * dur)
            # mask the known transits so flattening does not eat them
            ph = ((lc.time.value - t0) / P + .5) % 1 - .5
            msk = np.abs(ph * P) < dur
            fl = lc.flatten(window_length=int(round(win / cad)) | 1, mask=msk).remove_outliers(sigma_upper=3, sigma_lower=50)
            parts.append(fl.bin(time_bin_size=10 / 1440).remove_nans())
    finally:
        shutil.rmtree(d, ignore_errors=True)
    if not parts: return None
    t = np.concatenate([x.time.value for x in parts]); f = np.concatenate([x.flux.value for x in parts])
    ok = np.isfinite(t) & np.isfinite(f); t, f = t[ok], f[ok]
    ph = ((t - t0) / P + .5) % 1 - .5
    if np.sum(np.abs(ph * P) < dur / 2) < 5: return None          # no in-transit coverage
    dur_s = min(max(dur, 0.02), 0.5)
    G, L, Sc = make_views(t, f, P, t0, dur_s, Ms, Rs)
    return G, L, Sc, int(row.y), row.kind, tic, row.name

if __name__ == '__main__':
    tg = targets(); print(tg.kind.value_counts().to_string(), flush=True)
    rows = [r for _, r in tg.iterrows()]
    os.makedirs(RD, exist_ok=True)
    with mp.Pool(int(os.environ.get('NPROC', 8))) as pool:
        for i, r in enumerate(pool.imap_unordered(one, rows, chunksize=4)):
            if i % 250 == 0: print(i, flush=True)
    import glob
    res = [np.load(f) for f in glob.glob(f'{RD}/*.npz')]
    G, L, S = (np.stack([r[k] for r in res]) for k in ('G', 'L', 'S'))
    y, kind, tic, name = (np.array([r[k].item() for r in res]) for k in ('y', 'kind', 'tic', 'name'))
    np.savez_compressed(f'{D}/real_labels.npz', G=G, L=L, S=S, y=y, kind=kind, tic=tic, name=name)
    print('saved', len(y), dict(zip(*np.unique(kind, return_counts=True))))
