"""Which of two close stars dims?  For a signal (deep2/vet/{TIC}_{K}.json) and a neighbour NB within the same TESS pixel stamp:
fit the out-of-transit image with two pixel-integrated Gaussians at the catalogue positions (common shift, width, fluxes free), then
(a) fit the per-transit difference image with amplitudes on the two stars, (b) fit one free-position source to it and express its position
along the line target->neighbour (0 = target, 1 = neighbour).  Errors from bootstrap over transits.  Uses only pixel files already in
deep2/pix/{TIC}/.   Usage: deep2_pix2.py TIC K NB_TIC [other_signal_K ...]"""
import sys, os, json, glob, numpy as np, pandas as pd, warnings, lightkurve as lk
from astropy.coordinates import SkyCoord
from scipy.ndimage import median_filter
from scipy.optimize import least_squares
from scipy.special import erf
warnings.filterwarnings('ignore')
R = '/home/sarang/tess_search'; TIC, K, NB = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]); OTH = [int(x) for x in sys.argv[4:]]
v = json.load(open(f'{R}/deep2/vet/{TIC}_{K}.json')); P, t0, dur = v['P'], v['T0'], v['T14_h'] / 24
oth = [json.load(open(f'{R}/deep2/vet/{TIC}_{k}.json')) for k in OTH]
D = f'{R}/deep2/pix/{TIC}'; g = pd.read_csv(f'{D}/tic_120arcsec.csv').drop_duplicates('ID').set_index('ID'); A, B = g.loc[TIC], g.loc[NB]
sep = SkyCoord(A.ra, A.dec, unit='deg').separation(SkyCoord(B.ra, B.dec, unit='deg')).arcsec
print(f'TIC {TIC} (T={A.Tmag:.2f}) vs TIC {NB} (T={B.Tmag:.2f}), separation {sep:.1f}";  signal P={P:.6f} d, T14 {dur*24:.2f} h')
rng = np.random.default_rng(1)

def psf(xx, yy, x, y, s):                                   # pixel-integrated circular Gaussian, unit flux
    q = np.sqrt(2) * s
    return .25 * (erf((xx + .5 - x) / q) - erf((xx - .5 - x) / q)) * (erf((yy + .5 - y) / q) - erf((yy - .5 - y) / q))

rows = []; sk = []
for p in sorted(glob.glob(f'{D}/*_tp.fits')):
    try: tpf = lk.read(p, quality_bitmask='default')
    except Exception as ex: print('  read', os.path.basename(p), ex); continue
    t = tpf.time.value; fl = np.asarray(tpf.flux.value, float); ok = np.all(np.isfinite(fl.reshape(len(t), -1)), 1); t, fl = t[ok], fl[ok]
    cad = np.nanmedian(np.diff(t)); kk = int(round(0.75 / cad)) | 1; fl = fl - median_filter(fl, size=(kk, 1, 1), mode='nearest') + np.nanmedian(fl, 0)
    keep = np.ones(len(t), bool)
    for o in oth: keep &= np.abs(((t - o['T0']) / o['P'] + .5) % 1 - .5) * o['P'] > 0.75 * o['T14_h'] / 24
    t, fl = t[keep], fl[keep]
    xa, ya = [float(z) for z in tpf.wcs.world_to_pixel(SkyCoord(A.ra, A.dec, unit='deg'))]; xb, yb = [float(z) for z in tpf.wcs.world_to_pixel(SkyCoord(B.ra, B.dec, unit='deg'))]
    ph = (((t - t0) / P + .5) % 1 - .5) * P; nn = np.floor((t - t0) / P + .5).astype(int)
    intr, near = np.abs(ph) < 0.4 * dur, (np.abs(ph) > 0.7 * dur) & (np.abs(ph) < 3.5 * dur)
    diffs = np.array([np.nanmean(fl[near & (nn == k)], 0) - np.nanmean(fl[intr & (nn == k)], 0) for k in np.unique(nn[intr]) if (intr & (nn == k)).sum() >= 1 and (near & (nn == k)).sum() >= 6])
    if len(diffs) < 3: continue
    oot = np.nanmedian(fl[near], 0); yy, xx = np.indices(oot.shape); m = np.hypot(xx - xa, yy - ya) < 3.0
    bg = np.nanpercentile(oot, 10)
    rA = 1 / (1 + 10 ** (-.4 * (B.Tmag - A.Tmag))); rB = 1 - rA                    # flux ratio fixed to the catalogue (free ratio is degenerate with the shift)
    def r_oot(q): return (q[3] * (rA * psf(xx, yy, xa + q[0], ya + q[1], q[2]) + rB * psf(xx, yy, xb + q[0], yb + q[1], q[2])) + q[4] - oot)[m]
    tot = (oot - bg)[m].sum(); fo = least_squares(r_oot, [0, 0, 0.8, tot, bg], bounds=([-1.5, -1.5, 0.4, 0, -np.inf], [1.5, 1.5, 2.0, np.inf, np.inf])).x
    dx, dy, sg = fo[:3]; fA, fB = fo[3] * rA, fo[3] * rB; pa, pb = psf(xx, yy, xa + dx, ya + dy, sg)[m], psf(xx, yy, xb + dx, yb + dy, sg)[m]
    M = np.vstack([pa, pb, np.ones(m.sum())]).T; u = np.array([xb - xa, yb - ya]); L2 = u @ u
    def one(dd):
        d = dd.mean(0)[m]; amp = np.linalg.lstsq(M, d, rcond=None)[0]
        c2 = [np.sum((d - np.vstack([q, np.ones(m.sum())]).T @ np.linalg.lstsq(np.vstack([q, np.ones(m.sum())]).T, d, rcond=None)[0]) ** 2) for q in (pa, pb)]
        fr = least_squares(lambda q: (q[2] * psf(xx, yy, q[0], q[1], sg)[m] + q[3] - d), [xa + dx, ya + dy, d.sum(), 0]).x
        w = np.array([fr[0] - xa - dx, fr[1] - ya - dy]); return amp[0], amp[1], (w @ u) / L2, (w[0] * u[1] - w[1] * u[0]) / np.sqrt(L2), c2[0], c2[1], fr[0] - dx, fr[1] - dy
    r = np.array(one(diffs)); bs = np.array([one(diffs[rng.integers(0, len(diffs), len(diffs))]) for _ in range(150)]); e = bs.std(0)
    def sky(px, py):                                        # fitted pixel position -> offset from the target, arcsec East / North
        w_ = tpf.wcs.pixel_to_world(np.atleast_1d(px), np.atleast_1d(py)); return (w_.ra.deg - A.ra) * np.cos(np.radians(A.dec)) * 3600, (w_.dec.deg - A.dec) * 3600
    e0, n0 = sky(r[6], r[7]); eb, nb_ = sky(bs[:, 6], bs[:, 7]); sk.append([e0[0], n0[0], np.std(eb), np.std(nb_)])
    sec = int(tpf.sector); rows.append([sec, len(diffs), r[0] / fA, e[0] / fA, r[1] / fB, e[1] / fB, r[2], e[2], r[3] * 21, e[3] * 21, r[4], r[5], fA / fB])
    print(f'  S{sec}: {len(diffs)} transits; fitted flux ratio target/neighbour {fA/fB:.2f}, PSF sigma {sg:.2f} px, shift ({dx*21:+.1f}",{dy*21:+.1f}"); depth if on target {r[0]/fA*1e6:.0f}±{e[0]/fA*1e6:.0f} ppm, if on neighbour {r[1]/fB*1e6:.0f}±{e[1]/fB*1e6:.0f} ppm; '
          f'dip position along target->neighbour {r[2]:+.2f}±{e[2]:.2f} (0=target, 1=neighbour), across {r[3]*21:+.1f}±{e[3]*21:.1f}"; chi2(on target)/chi2(on neighbour) {r[4]/r[5]:.2f}', flush=True)
if rows:
    Z = np.array(rows); wm = lambda a, e: (np.sum(a / e ** 2) / np.sum(1 / e ** 2), 1 / np.sqrt(np.sum(1 / e ** 2)))
    dA, dB, s, c = wm(Z[:, 2], Z[:, 3]), wm(Z[:, 4], Z[:, 5]), wm(Z[:, 6], Z[:, 7]), wm(Z[:, 8], Z[:, 9])
    print(f'  combined ({len(Z)} sectors): two-star fit depth on target {dA[0]*1e6:.0f}±{dA[1]*1e6:.0f} ppm ({dA[0]/dA[1]:.1f} sig), on neighbour {dB[0]*1e6:.0f}±{dB[1]*1e6:.0f} ppm ({dB[0]/dB[1]:.1f} sig)')
    print(f'  dip position along target->neighbour: {s[0]:+.2f}±{s[1]:.2f}  = {s[0]*sep:+.1f}±{s[1]*sep:.1f}" from the target ({abs(s[0])/s[1]:.1f} sig from target, {abs(1-s[0])/s[1]:.1f} sig from neighbour); across the line {c[0]:+.1f}±{c[1]:.1f}"')
    S = np.array(sk); E, N = wm(S[:, 0], S[:, 2]), wm(S[:, 1], S[:, 3]); gg = g.reset_index(); gg['dE'] = (gg.ra - A.ra) * np.cos(np.radians(A.dec)) * 3600; gg['dN'] = (gg.dec - A.dec) * 3600
    gg['dist_to_dip'] = np.hypot(gg.dE - E[0], gg.dN - N[0]); gg['nsig'] = np.hypot((gg.dE - E[0]) / E[1], (gg.dN - N[0]) / N[1])
    print(f'  dip position on the sky relative to the target: {E[0]:+.1f}±{E[1]:.1f}" East, {N[0]:+.1f}±{N[1]:.1f}" North  ({np.hypot(E[0], N[0]):.1f}" away, {np.hypot(E[0]/E[1], N[0]/N[1]):.1f} sigma from the target)')
    print('  catalogue stars nearest to the dip position:'); print(gg.sort_values('nsig')[['ID', 'Tmag', 'dE', 'dN', 'dist_to_dip', 'nsig']].head(5).round(1).to_string(index=False))
    json.dump(dict(tic=TIC, k=K, nb=NB, dip_E=E[0], dip_E_err=E[1], dip_N=N[0], dip_N_err=N[1], sep=sep, depth_on_target_ppm=dA[0] * 1e6, depth_on_target_err=dA[1] * 1e6, depth_on_nb_ppm=dB[0] * 1e6, depth_on_nb_err=dB[1] * 1e6, pos_along=s[0], pos_along_err=s[1], pos_across_arcsec=c[0], pos_across_err=c[1],
                   sectors=[int(x) for x in Z[:, 0]], flux_ratio=float(np.median(Z[:, 12]))), open(f'{R}/deep2/vet/{TIC}_{K}.pix2.json', 'w'), default=float)
