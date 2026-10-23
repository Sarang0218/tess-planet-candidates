import json, urllib.request, concurrent.futures as cf, os
SMAP = json.load(open('/home/sarang/tess_search/data/sector_map.json'))
BASE = 'https://stpubdata.s3.amazonaws.com/tess/public/tid'
def url(tic, s):
    ts, sc = SMAP[str(s)]; t16 = f'{int(tic):016d}'
    return f'{BASE}/s{int(s):04d}/{t16[0:4]}/{t16[4:8]}/{t16[8:12]}/{t16[12:16]}/tess{ts}-s{int(s):04d}-{t16}-{sc}-s_lc.fits'
def exists(u):
    for _ in range(3):
        try:
            urllib.request.urlopen(urllib.request.Request(u, method='HEAD'), timeout=20); return True
        except urllib.error.HTTPError as e:
            if e.code in (403, 404): return False
        except Exception: pass
    return False
def available(tic):
    secs = sorted(int(s) for s in SMAP)
    with cf.ThreadPoolExecutor(16) as ex:
        ok = list(ex.map(lambda s: exists(url(tic, s)), secs))
    return [s for s, o in zip(secs, ok) if o]
def fetch(tic, secs, d):
    os.makedirs(d, exist_ok=True)
    def one(s):
        p = f'{d}/s{s}.fits'
        for _ in range(4):
            try: urllib.request.urlretrieve(url(tic, s), p); return p
            except Exception: pass
        return None
    with cf.ThreadPoolExecutor(8) as ex:
        return [p for p in ex.map(one, secs) if p]

QLP = 'https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:HLSP/qlp'
QLP_SECTORS = list(range(1, 14)) + list(range(27, 40)) + list(range(61, 70)) + list(range(87, 99))
def qlp_url(tic, s):
    t16 = f'{int(tic):016d}'
    return f'{QLP}/s{s:04d}/{t16[0:4]}/{t16[4:8]}/{t16[8:12]}/{t16[12:16]}/hlsp_qlp_tess_ffi_s{s:04d}-{t16}_tess_v01_llc.fits'
def qlp_available(tic):
    def ok(s):
        for _ in range(3):
            try:
                r = urllib.request.urlopen(urllib.request.Request(qlp_url(tic, s), method='HEAD'), timeout=30); return r.status == 200
            except urllib.error.HTTPError as e:
                if e.code in (403, 404): return False
            except Exception: pass
        return False
    with cf.ThreadPoolExecutor(16) as ex:
        res = list(ex.map(ok, QLP_SECTORS))
    return [s for s, o in zip(QLP_SECTORS, res) if o]
def qlp_fetch(tic, secs, d):
    os.makedirs(d, exist_ok=True)
    def one(s):
        p = f'{d}/q{s}.fits'
        for _ in range(4):
            try: urllib.request.urlretrieve(qlp_url(tic, s), p); return p
            except Exception: pass
        return None
    with cf.ThreadPoolExecutor(8) as ex:
        return [p for p in ex.map(one, secs) if p]
