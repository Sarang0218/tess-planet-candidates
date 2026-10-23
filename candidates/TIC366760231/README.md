# TIC 366760231: a new candidate ~2 R⊕ planet on a 3.14-day orbit around a nearby M dwarf

**Status:** new planet candidate from an independent search of public TESS data. Not confirmed.
Checked on 2026-10-21. The signal is not a TOI, a CTOI or a confirmed planet. It is not in any of the 132 SPOC TCE catalogs (single- and multi-sector, through sector 96). It is not in the public candidate lists of the independent projects *ahrazzle/exoplanets* (~530 candidates) or *tess-transit-hunter*. A web search on its TIC, 2MASS and Gaia identifiers returns nothing. ExoFOP lists no follow-up observations for this star.

## Host star (TIC v8.2 / Gaia DR3)

| | |
|---|---|
| IDs | TIC 366760231, 2MASS J22202826+7427466, Gaia DR3 2279461636455560576, UCAC4 823-025985 |
| RA, Dec (J2000) | 335.11793, +74.46294; pm 129.6, 30.5 mas/yr |
| TESS / Gaia G mag | 12.39 / 13.56 |
| Teff, R★, M★ | 3456 ± 157 K, 0.389 ± 0.012 R☉, 0.379 ± 0.020 M☉ (early–mid M dwarf) |
| Distance | 55.3 pc |
| Gaia RUWE | 1.20 |
| TIC contamination ratio | 0.107 |

## Candidate parameters
Ephemeris: linear fit to per-sector mid-transit times from all 9 sectors (2019–2025).

| Parameter | Value |
|---|---|
| **Period** | **3.137321 ± 0.000003 d** |
| **Mid-transit epoch** | **2458962.04155 ± 0.00061 BJD_TDB** |
| Depth | 1640 ± 220 ppm (free-density fit); 1915 ± 340 ppm (catalog-density fit, at mid-transit) |
| Duration T14 | 0.54–0.60 h (± 0.04–0.07) |
| Impact parameter | High: 0.97 (+0.23/−0.04) with the catalog stellar density |
| Planet radius | About 1.9–2.5 R⊕. The best catalog-density model has Rp/R★ = 0.045 (1.9 R⊕) at b = 0.93. The posterior tail to larger radii comes from the near-grazing degeneracy. |
| a/R★ (catalog density) | 16.7 ± 0.6 |
| Teq (A = 0), insolation | About 600 K, about 21 S⊕ |
| Detection | SNR 9.4 (all sectors combined); transit Δχ² ≈ 70 relative to no transit |
| Expected RV semi-amplitude | About 3 m/s for ~4 M⊕ |

**Note on the short duration:** a central transit would last about 1.4 h. The observed ~0.55 h therefore needs a high impact parameter. A fit with free stellar density prefers a box-like shape with an implied ρ★ ≈ 144 g cm⁻³. That fit is only marginally better (Δχ² = 2.5), so the catalog-density, high-b solution is fully consistent with the data.

## Evidence

| Test | Result | Figure |
|---|---|---|
| Discovery | GPU box-least-squares search with rotation prewhitening, bad-time masking and flare-aware cleaning, on 2-min sectors 17, 18, 19, 24 and 25 | — |
| **Independent confirmation in data not used for detection** | TESS-SPOC FFI sectors 58, 59 and 85 (2022–2024) at the ephemeris from the detection sectors (±1 h timing window scanned). Depth 1202 ± 233 ppm. 0 of 300 identical scans at random phases reach this depth (random mean 569 ± 101 ppm, ~6σ). | fig2 |
| Sector 52 (2-min, also unused for detection) | 1350 ± 455 ppm | fig3 |
| Depth consistency | All 9 sectors positive and mutually consistent: 941–1883 ppm | fig3 |
| Transit timing | 9 per-sector mid-times over 6 years fit a linear ephemeris (χ² = 12.8/7) | — |
| Period aliases | Best alias SNR 8.4 vs 5.0 for the runner-up (detection data) | — |
| Secondary eclipse | −93 ± 222 ppm (none) | — |
| Pixel-level location | Too shallow for single-sector difference-image centroids (S/N 1–2). Fitting the stacked difference images with the target's PSF gives 855 ± 178 ppm at the target position (4.8σ). | fig4 |
| **TRICERATOPS** (Giacalone et al. 2021) | **FPP = 0.264 ± 0.017, NFPP = 0.0000** (5 runs). P(TP) ≈ 0.64; the remainder is mostly STP (planet on an unresolved companion) and PTP/EBx2P. All nearby-star scenarios are ≈ 0. This meets the "likely planet" criterion (FPP < 0.5, NFPP < 0.001) but not validation (FPP < 0.015). | `triceratops_probs.csv` |

## Caveats
- The SNR is modest (~9). Single-sector difference images cannot localize the signal.
- The duration is short for the host star and needs a high impact parameter. Near-grazing geometry leaves the radius poorly constrained and keeps unresolved-companion scenarios open.
- Nearest Gaia neighbour: G = 16.3 at 11.9″ (it would need a 2–3% eclipse). TRICERATOPS assigns it negligible probability, but direct confirmation needs higher angular resolution.
- QLP light curves for sectors 57, 77–79, 84 and 86 were missing or unreadable, so those sectors are untested.
- Stellar parameters come from TIC v8.2 and have not been spectroscopically confirmed.

## Suggested follow-up
1. **High-resolution imaging** (speckle or AO). Ruling out companions within ~1″ is what would most reduce the FPP, toward validation.
2. **Ground-based transit photometry** at a predicted time (`upcoming_transits.md`, ±3 min). 1–2″ seeing resolves the 11.9″ neighbour and tests on-target depth and chromaticity. The depth (~1.6–1.9 ppt) is within reach of 0.6–1 m telescopes for a G = 13.6 star.
3. Reconnaissance spectroscopy to rule out a spectroscopic binary, then precise near-infrared RVs (K ≈ 3 m/s).

## Reporting route
ExoFOP requires community candidates to be published, either in a refereed journal or in an AAS Research Note that cites previously published methods. The methods used here are all standard: BLS (Kovács et al. 2002); transit model (Mandel & Agol 2002; batman, Kreidberg 2015); emcee (Foreman-Mackey et al. 2013); TRICERATOPS (Giacalone et al. 2021); TESS-SPOC FFI light curves (Caldwell et al. 2020). Another option is to contact the TESS TOI working group or exofop-support@ipac.caltech.edu.

## Files

| File | Contents |
|---|---|
| `fit_results.json`, `fit_results_freerho.json`, `chain*.npy` | MCMC fits, with and without the stellar-density prior |
| `ephemeris.json`, `upcoming_transits.md` | Linear ephemeris with covariance; predicted transits |
| `fig1_transit_fit*.png` | Folded 2-min data with model |
| `fig2_independent_FFI_sectors.png` | Detection data vs the independent FFI sectors |
| `fig3_per_sector_folds.png`, `fig4_difference_images.png` | Per-sector folds; difference images |
| `triceratops*.json`, `triceratops_probs.csv` | Statistical validation runs |
| `gaia_neighbors_120arcsec.csv` | Gaia DR3 neighbours |

Code: [`code/`](../../code/) (`search_core.py`, `deep_vet.py`, `extra_sectors_366760231.py`, `fit_366760231.py`, `tric_366760231.py`, `ephem_366760231.py`).
