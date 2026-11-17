# Code

Research scripts as used for this search. They are not a polished package: paths are hard-coded to the original working directories (`/home/sarang/tess_search/...`), so adjust them before running.

| Stage | Scripts |
|---|---|
| Download and pre-processing | `s3fetch.py`, `fetch_proc.py` (SPOC 2-min), `fetch_ffi.py` (TESS-SPOC FFI): flare masking, detrending, binning |
| Search | `gpu_bls.py` / `gpu_bls32.py` (GPU box-least-squares), `search_core.py` (cleaning, rotation prewhitening, iterative BLS, per-signal statistics), `search3.py` (batch runner), `badtimes.py` |
| Triage and vetting | `triage3.py` (flags and catalog cross-match), `ffi_tce.py` (TESS-SPOC FFI TCE lookup), `final_candidates.py` |
| ML vetter | `views.py`, `cnn_model.py`, `gen_train.py` (injections), `build_real.py` (TOI/EB labels), `build_domain_set.py`, `train_vetter3.py`, `vet3_gbm.py` |
| Deep vetting | `deep_vet.py` (2-min), `ffi_deep.py` (FFI independent-sector test), `vet_ffi_generic.py` (Gaia neighbours, difference-image PSF fit, TRICERATOPS), `single_event_check.py`, `extra_sectors_366760231.py`, `tric_366760231.py` |
| All-sector re-search (v1.1.0) | `deep2_avail.py`, `deep2_targets.py` (which stars have unused sectors), `deep2_fetch.py` (all sectors, per-sector noise), `deep2_core.py` (weighted multi-year BLS), `deep2_search.py`, `deep2_triage.py`, `deep2_null.py` (inverted and scrambled noise-floor test) |
| Vetting of all-sector signals (v1.2.0) | `deep2_mlset.py`, `deep2_vetml.py` (classifier trained on injections and scrambled light curves), `deep2_shortlist.py` and `deep2_vet.py` (per-sector depths, odd/even, secondary, other-phase scan, trapezoid shape), `deep2_nbcheck.py` (known TOIs/CTOIs on neighbouring stars), `deep2_tpfget.py` and `deep2_loc.py` (joint all-sector pixel localisation), `deep2_pix.py` (older per-sector difference images, kept for the figure), `deep2_fit.py` (batman + emcee, with and without a density prior), `deep2_tric.py` (TRICERATOPS on local pixel files), `deep2_evidence.py`, `deep2_writeup.py` |
| Extra checks (v1.3.0) | `deep2_indep.py` (independent re-reduction from raw aperture flux: depth, period, per-season transit times, depth by half and by cadence; `SRC=qlp` uses MIT QLP light curves instead), `deep2_archive.py` (old sky plates from SkyView at the star's present position, for background stars), `deep2_dvr.py` (reads the official SPOC/TESS-SPOC data-validation XML reports), `deep2_pix2.py` (which of two close stars dims). `deep2_loc.py` now has a `NOREG=1` option that takes star positions from the pixel-file coordinates instead of registering the image; use it when a much brighter star is near the target, whose glare can drag the registration (TIC 328936437). It only lists catalogue stars within about 60″ of the target. |
| TIC 366760231.02 | `vet_second_366760231.py`, `fit2_366760231.py` (two-signal fit, shared stellar density), `ephem2_366760231.py`, `pix2_366760231.py`, `tric2_366760231.py`, `evidence2_366760231.py` |
| Fits and figures | `fit_*.py` (batman + emcee), `ephem_*.py` (linear ephemerides), `evidence_198079261.py` |

Main dependencies: numpy, pandas, astropy, lightkurve, torch, batman-package, emcee, triceratops (a small compatibility shim was needed for pytransit ≥ 2.6), scikit-learn, astroquery.
