"""Training set for the noise-calibrated vetter of the deep pass.
For each real light curve: (1) block-shuffle the flux inside every sector (1-day blocks) -> no periodic signal survives, the noise does;
(2) with probability 0.65 inject a limb-darkened-like transit with a physical duration and an expected SNR drawn from U(4.5, 11);
(3) run the identical search. Signals matching the injection are positives; signals from non-injected (or unrelated) peaks are negatives.
Usage: deep2_mlset.py NSTARS [seed]"""
import sys, os, glob, json, time, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import deep2_core as dc
dc.limit_gpu(2.0)
D = '/home/sarang/tess_search'; LC = f'{D}/deep2/lc'; OUT = f'{D}/deep2/mlset'; os.makedirs(OUT, exist_ok=True)
N = int(sys.argv[1]); SEED = int(sys.argv[2]) if len(sys.argv) > 2 else 0
tg = pd.concat([pd.read_csv(f'{D}/data/{n}.csv') for n in ('targets', 'targets_allsky', 'targets_ffi_faint')]).drop_duplicates('ID').set_index('ID')
bad = np.load(f'{D}/data/bad_bins.npy')


def transit(t, P, t0, T14, depth, b, k):
    """Smooth trapezoid with ingress/egress set by k and b, mild limb-darkening curvature at the bottom."""
    x = np.abs(((t - t0) / P + .5) % 1 - .5) * P; tau = np.clip(T14 * k / max(1 - b ** 2, 0.05), 0.04 * T14, 0.5 * T14)
    flat = T14 / 2 - tau; z = np.clip((T14 / 2 - x) / tau, 0, 1); z = np.where(x < flat, 1.0, z)
    ld = 1 + 0.15 * (1 - np.clip(x / (T14 / 2), 0, 1) ** 2) - 0.10             # centre slightly deeper than limb
    return 1 - depth * z * ld


done = {os.path.basename(f)[:-5] for f in glob.glob(f'{OUT}/*.json')}; n_done = len(done); t_start = time.time()
while n_done < N:
    files = sorted(glob.glob(f'{LC}/*.npz')); rng0 = np.random.default_rng(SEED); rng0.shuffle(files)
    todo = [f for f in files if os.path.basename(f)[:-4] not in done]
    if not todo: time.sleep(120); continue
    n_before = n_done
    for fn in todo:
        tic = os.path.basename(fn)[:-4]; done.add(tic)
        try:
            rj = fn[:-4] + '.res.json'
            if not os.path.exists(rj): done.discard(tic); continue           # wait for the real search of this star
            rr = json.load(open(rj))
            # noise calibration uses quiet stars only: a real strong signal would leave non-periodic deep dips after shuffling
            if rr.get('status') != 'ok' or max([g['snr'] for g in rr['signals'] if dc.is_multi(g)], default=0) >= 7.2: continue
            d = np.load(fn)
            if len(d['t']) > 70000 or int(tic) not in tg.index: continue
            rng = np.random.default_rng(int(tic) % (2 ** 31) + SEED)
            t, f, e, s, _ = dc.prepare(d['t'], d['f'], d['s'], bad)
            if len(t) < 1000 or np.ptp(t) < 10: continue
            f2 = f.copy()
            for u in np.unique(s):                                    # block shuffle inside each sector
                idx = np.flatnonzero(s == u); cad = np.median(np.diff(t[idx])); nb = max(int(round(1.0 / cad)), 20)
                blocks = [idx[i:i + nb] for i in range(0, len(idx), nb)]; order = rng.permutation(len(blocks))
                f2[idx] = np.concatenate([f[blocks[j]] for j in order])
            st = tg.loc[int(tic)]; rad = float(st.rad) if np.isfinite(st.rad) else 0.4; ms = float(st['mass']) if np.isfinite(st['mass']) else 0.4
            inj = None
            if rng.random() < 0.65:
                P = float(np.exp(rng.uniform(np.log(0.6), np.log(30)))); t0 = t.min() + rng.uniform(0, P); b = rng.uniform(0, 0.9)
                aR = 4.2084 * P ** (2 / 3) * ms ** (1 / 3) / rad; T14 = P / np.pi * np.arcsin(min(1.0, np.sqrt(1 - b ** 2) / aR)) * rng.uniform(0.9, 1.1)
                sig0 = np.median(e); w = (sig0 / e) ** 2; intr = np.abs(((t - t0) / P + .5) % 1 - .5) * P < T14 / 2
                ntr = len(np.unique(np.floor((t[intr] - t0) / P + .5)))
                if intr.sum() >= 6 and ntr >= 3:
                    snr_t = rng.uniform(4.5, 11); depth = snr_t * sig0 / np.sqrt(w[intr].sum()); k = np.sqrt(depth)
                    if depth < 0.03:
                        f2 = f2 * transit(t, P, t0, T14, depth, b, k); inj = dict(P=P, t0=float(t0), T14=float(T14), depth=float(depth), b=float(b), snr_target=float(snr_t), ntr=int(ntr))
            res = dc.search_prepared(t, f2, e, s, nsig=1)           # top peak only: that is what the real search sees on a quiet star
            res.update(tic=int(tic), inj=inj, npts=int(len(t)), nsec=int(len(np.unique(s))), rad=rad, mass=ms, Tmag=float(st.Tmag))
            json.dump(res, open(f'{OUT}/{tic}.json', 'w')); n_done += 1
            if n_done % 50 == 0: print(f'[{time.strftime("%H:%M")}] {n_done} stars, {(time.time() - t_start) / max(n_done - len(done) + len(todo), 1):.1f} s/star', flush=True)
            if n_done >= N: break
        except Exception as ex:
            print('err', tic, str(ex)[:150], flush=True)
    if n_done == n_before: time.sleep(120)
print('mlset done', n_done)
