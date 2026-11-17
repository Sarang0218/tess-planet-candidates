"""Candidate page (README), transit table and hand-off copy for a deep-pass candidate.  Usage: deep2_writeup.py TIC K "star name"   (caveats in cand_{TIC}/caveats.md)"""
import sys, json, numpy as np, pandas as pd, os, shutil, re, warnings
warnings.filterwarnings('ignore')
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_sun
import astropy.units as u
from astropy.utils import iers; iers.conf.auto_download = False; iers.conf.auto_max_age = None
R = '/home/sarang/tess_search'; V = f'{R}/deep2/vet'; tic, k, name = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
loc = EarthLocation(lat=35.19 * u.deg, lon=129.10 * u.deg, height=50 * u.m); now = Time.now().jd
f = json.load(open(f'{V}/{tic}_{k}.fit.json')); fr = json.load(open(f'{V}/{tic}_{k}.fit_freerho.json')); v = json.load(open(f'{V}/{tic}_{k}.json')); tr = json.load(open(f'{V}/{tic}_{k}.tric.json'))
ch = np.load(f'{V}/{tic}_{k}.chain.npy'); nl = open(f'{V}/null_{tic}.log').read(); m = re.search(r'scrambled max SNR ([\d.]+), max SDE ([\d.]+); inverted ([\d.]+); real ([\d.]+) \(SDE ([\d.]+)\)', nl); nscr = len(re.findall(r',P=', nl))
g = pd.read_csv(f'{R}/deep2/pix/{tic}/tic_120arcsec.csv'); me = g[g.ID == tic].iloc[0]; c = SkyCoord(me.ra, me.dec, unit='deg'); P, T0 = f['period_d'][0], f['T0_BJD'][0]
lc = json.load(open(f'{V}/{tic}_{k}.loc.json'))
if 'sep' not in g: g['sep'] = np.hypot((g.ra - me.ra) * np.cos(np.radians(me.dec)), g.dec - me.dec) * 3600
if 'need' not in g: g['need'] = v['trap_depth_ppm'] * 1e-6 * (1 + 10 ** (-.4 * (g.Tmag - me.Tmag))) / 10 ** (-.4 * (g.Tmag - me.Tmag))
LOCALT = (f'The least-excluded other catalogued star that could mimic the signal is TIC {lc["alt_tic"]} (T = {lc["alt_Tmag"]:.1f}, {lc["alt_sep"]:.1f}″ away), disfavoured at {lc["alt_nsig"]:.1f} sigma.' if lc.get('alt_tic') else 'No other catalogued star within 60″ is bright enough to mimic the signal.')
nb = g[(g.ID != tic) & (g.need < 1) & (g.sep < 42)].sort_values('sep'); O = f'{R}/cand_{tic}'; os.makedirs(O, exist_ok=True); rows = []; n0 = int(np.ceil((now - T0) / P))
for n in range(n0, n0 + int(45 / P) + 1):
    x = ch[:, 1] + (n + round((T0 - 2457000 - np.median(ch[:, 1])) / P)) * ch[:, 0]; lo, hi = np.percentile(x, [16, 84]); sig = (hi - lo) / 2 * 1440
    t = Time(T0 + n * P, format='jd', scale='tdb'); tt = t + np.array([-1, 0, 1]) * (f['T14_h'][0] / 2 + 0.5) * u.hour
    alt = c.transform_to(AltAz(obstime=tt, location=loc)).alt.deg; sun = get_sun(tt).transform_to(AltAz(obstime=tt, location=loc)).alt.deg; rows.append((t, sig, alt, sun))
out = [f'# TIC {tic}.01: upcoming transits\n', f'Ephemeris: P = {P:.6f} d, T0 = BJD_TDB {T0:.5f}. Times are mid-transit in KST; altitude is from Busan at the start / middle / end of a window covering the transit ± 0.5 h.\n', '| Mid-transit (KST) | ± min | Altitude (deg) | Observable from Busan |', '|---|---|---|---|']
for t, sig, alt, sun in rows: out.append(f'| {(t.utc + 9 * u.hour).strftime("%Y-%m-%d %H:%M")} | {sig:.1f} | {alt[0]:.0f} / {alt[1]:.0f} / {alt[2]:.0f} | {"yes" if (alt.min() > 25 and sun.max() < -12) else "no"} |')
open(f'{O}/upcoming_transits.md', 'w').write('\n'.join(out) + '\n'); vis = [r for r in rows if r[2].min() > 25 and r[3].max() < -12]; e = lambda key: max(f[key][1:])
vm = f'{me.Vmag:.2f}' if np.isfinite(me.Vmag) else 'not listed'
md = f'''# TIC {tic}.01: planet candidate (not confirmed)

Found in the all-sector ("deep") re-search of TESS M dwarfs, 2026-11-16. Star: {name}.

| | |
|---|---|
| RA, Dec (J2000) | {c.ra.to_string(unit=u.hour, sep="hms", precision=2)}, {c.dec.to_string(sep="dms", precision=1, alwayssign=True)} |
| TESS mag / Gaia G / V | {me.Tmag:.2f} / {me.GAIAmag:.2f} / {vm} |
| Distance | {me.d:.1f} pc |
| Star (TIC) | Teff {me.Teff:.0f} K, R = {f["star"]["Rs"]:.3f} R☉, M = {f["star"]["Ms"]:.3f} M☉ |
| Period | {P:.6f} ± {e("period_d"):.6f} d |
| Mid-transit T0 | BJD_TDB {T0:.5f} ± {e("T0_BJD"):.5f} |
| Depth | {f["depth_ppm"][0]:.0f} (+{f["depth_ppm"][1]:.0f} / −{f["depth_ppm"][2]:.0f}) ppm |
| Duration T14 | {f["T14_h"][0]:.2f} ± {e("T14_h"):.2f} h |
| Rp/R★, impact parameter | {f["Rp_Rs"][0]:.4f} ± {e("Rp_Rs"):.4f}, {f["b"][0]:.2f} (+{f["b"][1]:.2f} / −{f["b"][2]:.2f}) |
| Planet radius | {f["Rp_Rearth"][0]:.2f} ± {e("Rp_Rearth"):.2f} R⊕ |
| a/R★, equilibrium temperature | {f["a_Rs"][0]:.1f}, {f["Teq_K"][0]:.0f} K (zero albedo) |
| Data | TESS sectors {", ".join(str(s) for s in f["sectors"])}; {v["ntr"]} transits |

## Evidence

- Search: SNR {float(m.group(4)):.1f}, SDE {float(m.group(5)):.1f}. Same search on {nscr} sector-scrambled copies with this signal removed: strongest periodic peak SNR {float(m.group(1)):.1f} (SDE {float(m.group(2)):.1f}); on the inverted light curve {float(m.group(3)):.1f}.
- Positive depth in {v["nsec_pos"]} of {v["nsec"]} sectors (chi2 {v["chi2_sec"]:.1f} for {v["nsec"]-1} dof).
- Odd vs even transits: {v["odd"][0]:.0f} ± {v["odd"][1]:.0f} vs {v["even"][0]:.0f} ± {v["even"][1]:.0f} ppm ({v["oe_sig"]:.1f} sigma). Secondary eclipse at phase 0.5: {v["sec_ppm"][0]:.0f} ± {v["sec_ppm"][1]:.0f} ppm.
- Where the dip is (joint fit of the difference images of {lc["nsec"]} sectors, {lc["ntr"]} transits): a point-source dip is detected at {lc["det_sigma"]:.1f} sigma, {np.hypot(lc["best_E"], lc["best_N"]):.1f}″ from the star's position at the time of observation, which is within the 4″ accuracy of the method (the star is {lc["target_nsig"]:.1f} sigma from the best position). Measured in the pixels as a dip of the star itself, the depth is {lc["depth_on_target_ppm"]:.0f} ± {lc["depth_on_target_err"]:.0f} ppm (the pixel measurement averages over the central 80 % of the transit and uses a simple PSF, so agreement within about 25 % is expected). {LOCALT}
- Stellar density from the transit shape with no prior: {fr["rho_star_gcc"][0]:.1f} (+{fr["rho_star_gcc"][1]:.1f} / −{fr["rho_star_gcc"][2]:.1f}) g/cc; catalogue value {f["star"]["rho_prior"][0]:.1f} g/cc.
- TRICERATOPS, {len(tr["FPP"])} runs: FPP = {np.mean(tr["FPP"]):.3f} ± {np.std(tr["FPP"]):.3f}, NFPP = {np.mean(tr["NFPP"]):.4f}. Background-star scenarios were not evaluated (no field-star population available), so the true FPP is somewhat higher.
- Neighbours within 42″ bright enough to mimic the signal (TIC): {"; ".join(f"TIC {int(r.ID)} at {r.sep:.1f}″, T = {r.Tmag:.2f}, would need a {r.need*100:.0f} % eclipse" for _, r in nb.iterrows()) if len(nb) else "none"}.
- Not in: TOI and CTOI lists (2026-11-16), Lastovka et al. 2026, T16 (Roth et al. 2026), DIAmante v1/v2, NEMESIS (Feliz et al. 2021), ahrazzle/comdex4 hobbyist lists. NOVELTY_EXTRA

## Caveats

{open(f"{O}/caveats.md").read().strip()}

## Files

- `EVIDENCE_TIC{tic}.png`: one-page figure
- `fit.json`, `fit_freerho.json`: transit-fit posteriors (density prior / free density)
- `upcoming_transits.md`: next 45 days, KST, with altitude from Busan
- `triceratops_probs.csv`: scenario probabilities of the last run
'''
extra = open(f'{O}/novelty_extra.txt').read().strip() if os.path.exists(f'{O}/novelty_extra.txt') else ''
open(f'{O}/README_TIC{tic}.md', 'w').write(md.replace('NOVELTY_EXTRA', extra))
for a, b in ((f'{V}/evid/EVIDENCE_TIC{tic}.png', f'{O}/EVIDENCE_TIC{tic}.png'), (f'{V}/{tic}_{k}.fit.json', f'{O}/fit.json'), (f'{V}/{tic}_{k}.fit_freerho.json', f'{O}/fit_freerho.json'), (f'{V}/{tic}_{k}.tric_probs.csv', f'{O}/triceratops_probs.csv'), (f'{V}/{tic}_{k}.fit.png', f'{O}/transit_fit.png')): shutil.copy(a, b)
H = f'/mnt/c/Users/saran/tess_candidates/deep_pass/TIC{tic}'; os.makedirs(H, exist_ok=True)
for fn in os.listdir(O):
    if fn not in ('caveats.md', 'novelty_extra.txt'): shutil.copy(f'{O}/{fn}', f'{H}/{fn}')
print(tic, 'written;', len(vis), 'observable transits from Busan in 45 d; next:', (vis[0][0].utc + 9 * u.hour).strftime('%m/%d %H:%M KST') if vis else '-', f'± {vis[0][1]:.1f} min' if vis else '', '| max altitude', round(90 - abs(35.19 - me.dec)), 'deg | T', me.Tmag, 'V', me.Vmag)
