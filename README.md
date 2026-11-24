# InSAR profiles distinguish localized deformation near a mapped lineament from broader subsidence in the Marvdasht Plain, southern Zagros, Iran

Derived observations and analysis tables for the revised manuscript. The original rasters, fault geometries, transects, and anomaly-center displacement series remain in their existing directories.

## Revised data products

The [revision directory](marvdasht_revision/) contains the reconstructed profiles, along-trace contrasts, model and coverage sensitivity tables, time series, comparisons with published velocities, and local spatial controls. SHA256SUMS records the checksum of each revised data file.

* Profiles and velocities use km for signed distance and mm/yr for velocity. Negative descending LOS velocity indicates motion away from the satellite. The primary descending frame is 137D_05960_131313.
* The ascending frame is 028A_06017_131313. Phase samples are patch medians in radians. Empty site values indicate insufficient valid pixels. The full and pilot tables overlap. The full table is a partial sample of the archive.
* T1_interval_components.csv reports relative displacement in mm, conditional on negligible differential northward motion. The comparison uses the same target and southwestern control throughout each interval. Positive up and east indicate uplift and eastward displacement.
* T1_local_reference_test.csv supports Appendix Table A1. The relaxed field is 0 for initial terrain matching, 1 for relaxed terrain matching, and 2 for the separate land-cover comparison. A translation selected by both methods appears twice. East and north translations are in km.
* endpoint_sensitivity.csv compares archived processing, unfiltered descending displacement, and GACOS for 2017-2019. short_interval.csv tests regional VMF3 sensitivity for 2023-2024. Height offset 0 is the baseline. Neither comparison gives a formal confidence interval.
* The published vertical and east-west velocity grids used for comparison are distributed by [Payne et al.](https://doi.org/10.5281/zenodo.13754200). Other inputs include COMET-LiCSAR observations, ESA WorldCover 2021 v200, and [TU Wien VMF3](https://vmf.geo.tuwien.ac.at/products.html). Consult the original providers for their data and terms.

## Reproduction

Install the pinned dependencies with `python -m pip install -r requirements.txt`. Run `python scripts/reproduce_profiles.py` from the repository root. The script reads the included georeferenced velocity rasters, mapped traces, archived transects, and displacement series. It writes profiles, along-trace contrasts, covariance model comparisons, temporal sensitivity estimates, and figures to `results/`.

Run `python scripts/reproduce_phase.py` to reconstruct both ascending closure-filtered phase series from the archived patch medians.

The additional geometry, atmospheric, and reference-site tables are derived products. Reprocessing these analyses from SAR observations requires the original LiCSAR and LiCSBAS archives and the atmospheric products described above.
