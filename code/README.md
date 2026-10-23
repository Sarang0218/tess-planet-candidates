# Code

Research scripts as used for this search. They are not a polished package: paths are hard-coded to the original working directories (`/home/sarang/tess_search/...`), so adjust them before running.

| Stage | Scripts |
|---|---|
| Download and pre-processing | `s3fetch.py`, `fetch_proc.py` (SPOC 2-min), `fetch_ffi.py` (TESS-SPOC FFI): flare masking, detrending, binning |
| Search | `gpu_bls.py` / `gpu_bls32.py` (GPU box-least-squares), `search_core.py` (cleaning, rotation prewhitening, iterative BLS, per-signal statistics), `search3.py` (batch runner), `badtimes.py` |
| Triage and vetting | `triage3.py` (flags and catalog cross-match), `ffi_tce.py` (TESS-SPOC FFI TCE lookup), `final_candidates.py` |
| ML vetter | `views.py`, `cnn_model.py`, `gen_train.py` (injections), `build_real.py` (TOI/EB labels), `build_domain_set.py`, `train_vetter3.py`, `vet3_gbm.py` |
| Deep vetting | `deep_vet.py` (2-min), `ffi_deep.py` (FFI independent-sector test), `vet_ffi_generic.py` (Gaia neighbours, difference-image PSF fit, TRICERATOPS), `single_event_check.py`, `extra_sectors_366760231.py`, `tric_366760231.py` |
| Fits and figures | `fit_*.py` (batman + emcee), `ephem_*.py` (linear ephemerides), `evidence_198079261.py` |

Main dependencies: numpy, pandas, astropy, lightkurve, torch, batman-package, emcee, triceratops (a small compatibility shim was needed for pytransit ≥ 2.6), scikit-learn, astroquery.
