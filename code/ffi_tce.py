"""Look up TESS-SPOC FFI TCE periods for a star (downloads its DV XML reports on demand; cached)."""
import json, os, re, glob, warnings
warnings.filterwarnings('ignore')
IDX = json.load(open('/home/sarang/tess_search/data/spoc_ffi_tce_index.json')); CACHE = '/home/sarang/tess_search/data/spoc_dv/xml'
def periods(tic):
    if str(int(tic)) not in IDX: return []
    from astroquery.mast import Observations
    out = []
    for rng, _ in IDX[str(int(tic))]:
        fn = f'hlsp_tess-spoc_tess_phot_{int(tic):016d}-{rng}_tess_v1_dvr.xml'; loc = glob.glob(f'{CACHE}/**/{fn}', recursive=True)
        if not loc:
            try:
                Observations.download_file(f'mast:HLSP/tess-spoc/{rng}/target/{int(tic):016d}'[:0] + f'mast:HLSP/tess-spoc/{rng}/target/' + '/'.join([f'{int(tic):016d}'[i:i + 4] for i in (0, 4, 8, 12)]) + '/' + fn, local_path=f'{CACHE}/{fn}')
            except Exception: pass
            loc = glob.glob(f'{CACHE}/**/{fn}', recursive=True)
        if loc:
            s = open(loc[0]).read()
            for m in re.finditer(r'<dv:planetResults[^>]*planetNumber="(\d+)".*?orbitalPeriodDays" uncertainty="[^"]*" value="([^"]+)"', s, re.S):
                out.append((rng, int(m.group(1)), float(m.group(2))))
    return out
if __name__ == '__main__':
    import sys; os.makedirs(CACHE, exist_ok=True); print(periods(int(sys.argv[1])))
