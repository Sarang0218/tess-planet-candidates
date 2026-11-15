# TIC 366760231.02: a second transit signal (1.1434 d, ~1 R⊕) on the host of the 3.137-d candidate

**Status:** planet *candidate*, found 2026-11-16 in a deeper all-sector re-search and first posted publicly the same day (v1.1.0). Not confirmed.
Checked 2026-11-16: TIC 366760231 is not in the TOI list (8,175 rows) or the CTOI list (5,139 rows) downloaded that day, not in any SPOC TCE catalog, and not in the Lastovka et al. 2026 candidate or TCE tables. A web search on its TIC, 2MASS and Gaia identifiers returns nothing.

## Parameters (both signals, one shared stellar density)

| | .01 (known candidate) | **.02 (new)** |
|---|---|---|
| Period (d) | 3.137316 ± 0.000005 | **1.143425 ± 0.000003** |
| Mid-transit (BJD_TDB) | 2459404.40233 ± 0.00084 | **2459403.40411 ± 0.00152** |
| Depth (ppm) | 1810 (+260/−210) | **585 ± 66** |
| Duration T14 (h) | 0.67 ± 0.04 | **1.01 (+0.05/−0.06)** |
| Rp/R★ | 0.048 (+0.008/−0.004) | 0.0226 ± 0.0014 |
| Impact parameter | 0.94 ± 0.02 | 0.23 (+0.19/−0.16) |
| Radius (R⊕) | 2.03 (+0.33/−0.17) | **0.96 ± 0.06** |
| a/R★, inclination | 16.7, 86.8° | 8.5, 88.4° |
| Teq (A = 0), insolation | 598 K, 21 S⊕ | 838 K, 82 S⊕ |

Stellar density: TIC prior 9.1 ± 0.9 g cm⁻³; from the two transit shapes alone (no prior) 5.9 (+5.3/−2.4). Both durations fit the same star.
The .01 ephemeris here supersedes the published one (3.137321 d); it predicts transits about 6 min earlier in Nov 2026, within the combined errors.

## Evidence for .02

| Test | Result |
|---|---|
| Detection | All-sector weighted BLS (9 sectors, 2019–2025): SNR 8.5, SDE 8.0, after removing the .01 transits. First pass (5 sectors) had it at SNR 6.3, below threshold. |
| Noise floor | 100 searches of the same light curve with each sector's flux cyclically shifted: strongest SNR 7.0 (median 6.3). Inverted light curve: 7.1. |
| Raw-file re-detrend | Depth 508–557 ppm (8.4–9.1σ) for detrending windows 0.5 / 0.75 / 1.5 d, with and without masking the signal. |
| Per sector | Positive in all 9 sectors (284–821 ppm each), χ² = 7–8 for 8 dof. |
| Odd / even | 527 ± 86 vs 545 ± 88 ppm. |
| Secondary eclipse | 79 ± 62 ppm (none). |
| Other phases | Same box at 181 other phases: scatter 67 ppm, largest 123 ppm. |
| Pixel level | Per-transit difference images fitted with the target's own image: 455 ± 60 ppm on target, positive in 9/9 sectors; 7.1σ against the same fit at 41 other phases. |
| TRICERATOPS (3 runs) | FPP 0.098 ± 0.004, NFPP 0.069. TP 0.77, PTP 0.13, STP 0.03; nearby-star scenarios 0.07 (TIC 366760229 at 24″, 366760220 at 29″). Background-star scenarios were skipped (TRILEGAL unavailable), so the FPP is a lower bound. |
| Not a TESS-orbit harmonic | 12 × P = 13.721 d vs the 13.66-d orbit: 7 cycles of drift over the baseline. |

## Caveats
- SNR 8.5 on a 585 ppm dip. It needs independent confirmation; TESS sectors after 85 would be the first test.
- NFPP 0.069 fails the "likely planet" cut (0.001). The pixel test argues for the target but cannot separate the 11.6″ neighbour.
- The folded dip looks slightly wider on the egress side than the model; the free-density fit allows T14 up to ~1.27 h.
- Too shallow for small ground-based telescopes. On 2026-11-17 KST both signals transit together (.01 mid 23:06 ± 5 min, .02 mid 23:11 ± 7 min).

Files (this folder): `ephemeris_final.json`, `fit2_results.json`, `fit2_results_freerho.json`, `upcoming_01.md`, `upcoming_02.md`, `triceratops_02_probs_*.csv`, `EVIDENCE_TIC366760231_02.png`, `fig_two_planets.png`.
Code: `code/deep2_*.py` (all-sector search), `vet_second_366760231.py`, `fit2_366760231.py`, `ephem2_366760231.py`, `pix2_366760231.py`, `tric2_366760231.py`, `evidence2_366760231.py`.
