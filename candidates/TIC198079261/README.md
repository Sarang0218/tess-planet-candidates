# TIC 198079261: a candidate Earth-sized planet on a 2.82-day orbit around a nearby ultracool (M6–M7) dwarf

**Status:** new planet candidate from an independent search of public TESS full-frame-image (FFI) light curves. Not confirmed.
Checked on 2026-10-21. The signal is not a TOI, a CTOI or a confirmed planet. It is not in any SPOC 2-min TCE catalog (132 runs through sector 96), and it is not a TESS-SPOC FFI TCE: the star has no data-validation products in any of the 51 TESS-SPOC DV releases. It is not in the public candidate lists of *ahrazzle/exoplanets* or *tess-transit-hunter*. A web search on its TIC, 2MASS and Gaia identifiers returns nothing. ExoFOP lists no follow-up observations.

## Host star (TIC v8.2 / Gaia DR3)

| | |
|---|---|
| IDs | TIC 198079261, 2MASS J04160640−5935386, Gaia DR3 4678628654309752576, WISE J041606.56−593537.5 |
| RA, Dec (J2000) | 64.02671, −59.59404; pm 116.7, 93.5 mas/yr |
| TESS / Gaia G / K mag | 15.03 / 16.66 / 12.02 |
| Teff, R★, M★ | 2825 ± 157 K, 0.123 ± 0.004 R☉, 0.099 ± 0.020 M☉ (ultracool dwarf, ~M6–M7) |
| Distance | 31.55 ± 0.09 pc (Gaia parallax 31.64 mas) |
| Gaia RUWE / non_single_star | 1.02 / 0, so no astrometric sign of binarity |
| TIC contamination ratio | 0.007 |

## Candidate parameters
MCMC (batman + emcee) on unbinned TESS-SPOC 200-s data from S62–S69, with a stellar-density prior from the TIC. The ephemeris comes from a linear fit to per-sector mid-times over 20 sectors (2018–2023), with uncertainties scaled by √(χ²/dof).

| Parameter | Value |
|---|---|
| **Period** | **2.817636 ± 0.000018 d** |
| **Mid-transit epoch** | **2460071.6821 ± 0.0024 BJD_TDB** |
| Depth (model, at mid-transit) | 6520 ± 1440 ppm (box depth over all sectors: 7314 ppm) |
| Duration T14 | 0.66 ± 0.08 h |
| Rp/R★ | 0.076 ± 0.009 |
| Impact parameter | 0.44 (+0.23/−0.28) |
| a/R★ | 31.6 ± 2.3 |
| **Planet radius** | **1.02 ± 0.13 R⊕** |
| **Teq (A = 0), insolation** | **355 ± 13 K, 2.7 ± 0.4 S⊕** |
| Detection SNR (all sectors) | 13.3 |
| Expected RV semi-amplitude | ~2–3 m/s for ~1 M⊕ (0.1 M☉ host) |

**Stellar density from the transit shape alone:** a fit without the density prior gives ρ★ = 90 (+64/−66) g cm⁻³, consistent with the TIC value of 75 ± 17. In other words, the transit shape and duration fit this specific ultracool host.

## Evidence

| Test | Result | Figure |
|---|---|---|
| Discovery | GPU BLS with rotation prewhitening, flare masking and bad-time masking, on TESS-SPOC FFI sectors 28–39. The vetter v3 ensemble ranked it at the top of the unknown signals (score 1.00). | — |
| **Independent confirmation** | 13 sectors not used for detection (S2–12 at 30-min cadence, S62–69 at 200-s), at the detection-only ephemeris with the ±3.3 h timing window scanned. Depth 5314 ± 824 ppm. 0 of 150 identical scans at random phases reach this depth (random mean 2902 ± 495 ppm, ~4.9σ). | fig2 |
| Depth consistency | 19 of 20 sectors show a dip of 4900–12700 ppm, each consistent with ~7–8 ppt. S68 shows 500 ± 3300 ppm (2.5σ low). | — |
| Odd vs even | 6936 ± 775 vs 7696 ± 779 ppm (0.7σ) | — |
| Secondary eclipse | None detected | — |
| **Pixel-level location** | The stacked difference images over 5 TESS-SPOC FFI TPF sectors, fitted with the target's own PSF, give **7195 ± 903 ppm (8.0σ)** at the target position. That equals the light-curve depth, so the full dip comes from the target's PSF. | — |
| Gaia neighbours | The nearest star able to mimic the dip (G = 18.6, needing a 5% eclipse) is 26″ away (1.2 TESS px). Every capable neighbour is ≥ 26″ away. | `gaia_neighbors_120arcsec.csv` |
| **TRICERATOPS** (Giacalone et al. 2021) | **FPP = 0.159 ± 0.008, NFPP = 0.005** (5 runs). P(TP) ≈ 0.73. The remainder is mostly PTP/STP (an unresolved companion); EB scenarios total ≈ 0.04. The formal "likely planet" class needs NFPP < 0.001. The pixel-level on-target test above, which TRICERATOPS does not use, independently argues against nearby-star scenarios. | `triceratops_probs.csv` |

## Caveats
- The star is faint (T = 15.0). The detection relies on FFI photometry, and the per-sector SNR is only ~3–4.
- Per-sector mid-times scatter more than their formal errors (χ²/dof = 5.1). Ephemeris uncertainties are inflated accordingly.
- An unresolved bound companion within ~1″ cannot be excluded without high-resolution imaging. This is the main source of the remaining FPP.
- Stellar parameters come from TIC v8.2 (photometric relations) and have not been spectroscopically confirmed.

## Why it matters, if confirmed
An Earth-sized (1.0 R⊕), temperate-to-warm (~355 K) planet transiting a 0.12 R☉ ultracool dwarf at 31.6 pc would be a TRAPPIST-1-like system. The ~0.7% transit depth is deep enough for transmission spectroscopy with JWST.

## Suggested follow-up
1. **Ground-based transit photometry in the red/near-IR** (e.g. SPECULOOS, TRAPPIST, MEarth-South, or a 1–2 m telescope in I/z) at a predicted time (`upcoming_transits.md`, ±14 min). A ~0.7% dip on an I ≈ 14–15 star is readily detectable from the ground. It would confirm the source at arcsecond resolution and refine the ephemeris.
2. High-resolution imaging (AO/speckle) to rule out companions within ~1″.
3. NIR spectroscopy to confirm the spectral type, then NIR RVs (K ≈ 2–3 m/s).

## Reporting route
As with the other candidates: ExoFOP accepts only published community candidates (a refereed journal or an AAS Research Note citing published methods), or the target can be raised with the TESS TOI working group / exofop-support@ipac.caltech.edu.

## Files

| File | Contents |
|---|---|
| `fit_results.json`, `fit_results_freerho.json`, `chain*.npy` | MCMC fits, with and without the density prior |
| `ephemeris.json`, `upcoming_transits.md` | Linear ephemeris (scaled errors); predicted transits |
| `fig1_transit_fit*.png` | Folded 200-s data with model |
| `fig2_detection_vs_independent_sectors.png` | Detection sectors vs the 13 independent sectors |
| `triceratops_probs.csv`, `vet_summary.json` | Statistical validation and on-target pixel test |
| `gaia_neighbors_120arcsec.csv` | Gaia DR3 neighbours |

Code: [`code/`](../../code/) (`fetch_ffi.py`, `search_core.py`, `ffi_deep.py`, `vet_198079261.py`, `fit_198079261.py`, `ephem_198079261.py`, `train_vetter3.py`).
