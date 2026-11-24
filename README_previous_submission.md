# Data supporting: Distinguishing Groundwater-Driven Subsidence from Fault Slip in InSAR Velocity Fields — A Decadal Sentinel-1 Test across the Kazerun Fault System, Zagros, Iran

Data repository for a manuscript submitted for review (author information withheld for double-blind review).

This repository holds the derived data products behind the results, figures, and tables of
the manuscript submitted to *Natural Hazards*. It does not include raw Sentinel-1 SAR data,
raw interferograms, or the processing scripts; these are publicly available or reproducible
as described below.

## Study area and InSAR processing

- Sentinel-1 frame: `137D_05960_131313` (descending, relative orbit 137).
- Interferograms: COMET-LiCS Sentinel-1 InSAR portal (Lazecky et al., 2020).
- Time series inversion: LiCSBAS v1.15.1 (Morishita et al., 2020), weighted least squares,
  212 acquisitions, 27 October 2014 to 31 May 2024, 886 interferograms after quality control.
- Corrections: first-order planar ramp (`-r 1`) and elevation-correlated component
  (`--hgt_linear`), estimated jointly; 2 km spatial / 49-day temporal filter (LiCSBAS step 16).
- Grid: 1532 x 1349 pixels, ~200 m posting, geographic coordinates (WGS84 / EPSG:4326).
  Upper-left corner (pixel center convention as stored): latitude 31.7305555°N,
  longitude 50.6744444°E; pixel spacing -0.0020000° (latitude, north to south) and
  +0.0020000° (longitude).
- Reference pixel: 31.202°N, 52.383°E (row 269, column 886 in the grid above).
- All velocities are line-of-sight (LOS); negative = motion away from the satellite.

Full reprocessing from the raw interferograms requires only public tools (LiCSBAS,
freely available at https://github.com/yumorishita/LiCSBAS) and the parameters above; no
custom code is required to reproduce the published rasters, and none is distributed here.

## Contents

```
rasters/
  vel_los_corrected_mm_yr.tif     Corrected, filtered, masked LOS velocity (Fig. 2b, Fig. 3)
  vel_los_uncorrected_mm_yr.tif   Uncorrected network-inversion LOS velocity (Fig. 2a)
  vstd_bootstrap_mm_yr.tif        Bootstrap velocity standard deviation (Fig. 4b)
  coh_avg.tif                     Mean interferometric coherence (Fig. 4a)
  mask.tif                        LiCSBAS step-15 quality mask (1 = valid)
  All GeoTIFF, float32, EPSG:4326, nodata = NaN.

transects/
  FaultA_T1..T4_swath12km.csv     Fault A fault-perpendicular transects (Fig. 5)
  FaultB_T1..T4_swath12km.csv     Fault B fault-perpendicular transects (Fig. 6b, Table 3)
  Columns: distance_km (positive east), lon, lat, los_vel_mm_yr, bootstrap_sd_mm_yr.
  12-km swath mean (5 parallel lines at -6,-3,0,+3,+6 km offset), 1-km sampling.

faults/
  FaultA_FaultB.geojson                  The two profiled basement-fracture lineaments
  GSI_fault_segments_frame_subset.geojson 596 fault segments within the study frame,
                                          subset of the Geological Survey of Iran digital
                                          fault dataset (1:250,000 map series), for
                                          orientation reference only; attribute "type"
                                          gives the source classification.

tables/
  ISC_reviewed_events_in_frame_1925_2024.csv   All Reviewed ISC Bulletin events within the
                                                study frame (Table 1 area), 1925-2024.
  ISC_events_M5.8plus_Table2.csv               Table 2 of the manuscript.
  faultB_transect_stats.csv                    Per-transect summary statistics (Table 3).
  faultB_models_v7.csv                         Fault-step / dislocation / bowl model
                                                comparison statistics (Table 3, Section 4.3).
  v7_summary.json                              Noise-window, time-series, and seismicity
                                                summary statistics quoted in the text.

timeseries/
  anomaly_center_29.896N_52.873E_5x5mean.csv   Cumulative LOS displacement at the Fault B
                                                anomaly center, corrected and uncorrected,
                                                and at the reference pixel (Fig. 8).
```

## Earthquake catalogue and magnitude convention

Event locations, origin times, and depths are from the Reviewed ISC Bulletin
(International Seismological Centre, https://doi.org/10.31905/D808B830; query: 29.0-31.8°N,
50.6-53.8°E, 1900-2024; events after 1 September 2024 not yet reviewed by the ISC and
excluded).

Each event in the ISC bulletin can carry several magnitudes reported by different networks
and of different types (mb, MS, ML, Mw, ...). A single, documented magnitude is assigned to
each event as follows, in order of preference:

1. Moment magnitude (Mw), preferring the Global CMT catalogue, then NEIC/USGS, then GFZ,
   then any other network reporting Mw.
2. If no Mw is reported, the ISC's own mb or MS magnitude.
3. If neither is available, the first magnitude listed for that event in the bulletin.

This rule is applied uniformly to every event in `ISC_reviewed_events_in_frame_1925_2024.csv`
and to the magnitude counts and thresholds quoted in the manuscript text (Section 4.1,
Section 5.3) and in Table 2. It is not the ISC's own "preferred magnitude" field, which is
not exposed in this bulletin export; it is a documented substitute chosen to avoid
mixing magnitude scales and to avoid systematically selecting the largest of several
reported values for the same event.

## Fault dataset

The 596-segment subset (`GSI_fault_segments_frame_subset.geojson`) is a portion of the
Geological Survey of Iran's digital fault dataset (1:250,000 geological map series),
reprojected from UTM Zone 39N to geographic coordinates, retained here for the coordinates
of the segments actually used or shown in the manuscript. It carries no fault names or
slip-sense attributes beyond the source's type classification and is provided for
reproducing Fig. 1 and Fig. 3 only; it is not a substitute for the original dataset.

## Local seismicity catalogue

The supplementary local catalogue (Institute of Geophysics, University of Tehran / Iranian
Seismological Center, 1996-2017) is not redistributed here; it was obtained directly from
the network and can be requested from IRSC (http://irsc.ut.ac.ir).

## License

Data files: CC BY 4.0. Please cite the manuscript.

## References

- Lazecky, M., et al. (2020). Remote Sensing, 12(15), 2430.
- Morishita, Y., et al. (2020). Remote Sensing, 12(3), 424.
- International Seismological Centre (2026). ISC Bulletin, Reviewed. https://doi.org/10.31905/D808B830
