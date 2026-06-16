Wind effects on coral bleaching severity — analysis code
Code accompanying Lapenis & Jiang (Nature Climate Change, submitted): analysis of
how sustained background mean wind speed relates to coral-bleaching severity, the
supporting ocean-carbonate (BUD) box model, and the CMIP6 projection of future
mean-wind change over reefs.
This repository contains analysis and modelling code only. Input datasets are
obtained from the public sources listed in the manuscript's Data Availability
statement; derived datasets are deposited separately (see Data & code availability
below). No local file paths are embedded — every script reads from `data/input/` and
writes to `data/output/` relative to the repository root.
Repository structure
```
coral-wind-bleaching/
├── R/
│   # --- predictor & response construction ---
│   ├── 07_tcp_pdi_holland.R               Tropical-cyclone power (TCP/PDI) within 400 km,
│   │                                      Holland radial-wind profile, IBTrACS v04r01
│   ├── 08_blended_wind_monthly_extraction.R  Monthly mean wind (1991–2010) to
│   │                                      reef cells (NOAA blended winds) — helper module (functions only)
│   ├── 09_severity_binning_gcbd.R         GCBD bleaching severity → 3 ordered classes (5/30/75)
│   ├── 03_omega_extraction.R              Aragonite saturation (Ω) at 0/50/100 m onto reefs
│   # --- statistical analysis ---
│   ├── 01_glo_olr_wind_turbidity.R        GLO ordinal regression (FE & ME), wind × turbidity
│   ├── 06_glo_shapley_mundlak.R           GLO Shapley R² decomposition + Mundlak within-region
│   ├── 02_saf_mediation_residual_wind.R   SAF FE+ME, residual-wind test, mediation
│   ├── 10_saf_shapley.R                   SAF Shapley decomposition (RECONSTRUCTED; see header)
│   ├── 04_bud_box_model.R                 Box Upwelling–Diffusion carbonate model
│   ├── 11_glo_selection_robustness.R     Wind-OR robustness: co-cell cluster + one-per-cell
├── python/
│   └── 05_cmip6_wind_delU_per_reef.py     CMIP6 ΔU per reef + ensemble + odds maps
├── data/{input,output}/                   (place inputs here; outputs land here)
├── docs/glo_dataset_construction.md       How the gridded (GLO) dataset was built
├── run_all.R          master run-order script (sources the R scripts in order)
├── requirements.txt   LICENSE.txt   .gitignore
```
(File-number prefixes are stable IDs, not run order; the pipeline order is: predictor/response
construction → statistical analysis → projection.)
> One file needs a reader's attention before citing: `10_saf_shapley.R` reconstructs the
> SAF decomposition by adapting the GLO routine — run it against the deposited SAF dataset
> and confirm it returns the reported wind share (~18%); the script prints the reproduced
> share and flags whether it matches. `08_blended_wind_monthly_extraction.R` is a helper
> module — sourcing it only defines `interpolate_windspeed()` and
> `extract_monthly_mean_wind()`, which you call with a NOAA blended-wind field you have loaded.
How the scripts map to the manuscript
Script	Produces	Manuscript location
`07_tcp_pdi_holland.R`	Cumulative TCP/PDI per site (6/12-mo, 1993–2020)	Supplementary Note 2 (TCP)
`08_blended_wind_monthly_extraction.R` (helper module)	Monthly mean wind per reef site	Supplementary Note 1 (wind metrics)
`09_severity_binning_gcbd.R`	3-class GLO severity response (5/30/75)	Methods / `docs/glo_dataset_construction.md`
`03_omega_extraction.R`	Ω at 0/50/100 m for each reef cell	Supplementary Note 3; Suppl. Fig.
`01_glo_olr_wind_turbidity.R`	GLO ordinal-regression coefficients, FE/ME	Main GLO results figure
`06_glo_shapley_mundlak.R`	GLO Shapley shares (the ~18% wind share) + Mundlak	Results; Supplementary Note 4
`02_saf_mediation_residual_wind.R`	SAF model-improvement, residual-wind, mediation	Main SAF results figure
`10_saf_shapley.R`	SAF Shapley shares (Thermal/Depth/Wind)	Supplementary Note 4 (reconstructed)
`04_bud_box_model.R`	BUD Ω/pCO₂/pH vs wind across (m,n,q)	Supplementary Note 3; Suppl. Fig.
`11_glo_selection_robustness.R`	Wind OR: cluster-bootstrap + one-per-cell; 6/12-mo/long-term window comparison	Supplementary robustness table
`05_cmip6_wind_delU_per_reef.py`	Per-reef ΔU, ensemble means, odds maps	CMIP6 Note; main-text odds map
Running
R (≥ 4.2). Install the packages listed in each script header (consolidated in
`r_requirements.txt`), then run from the repository root, e.g.
```r
setwd("/path/to/coral-wind-bleaching")
source("R/01_glo_olr_wind_turbidity.R")
```
To run the full pipeline in the intended order, use the master script
`source("run_all.R")` from the repository root. The scripts no longer
auto-install packages: if a required package is missing they stop with an
informative message, so install `r_requirements.txt` first (or use `renv`).
Each script begins with a CONFIG / USER-SETTINGS block; place the corresponding input
file in `data/input/` (filenames are documented in the header) and outputs are written
to `data/output/`. For an exact-reproducibility record, pin package versions from your
`sessionInfo()` (or capture the environment with `renv`).
Python (≥ 3.9): `pip install -r requirements.txt`. The CMIP6 script streams CMIP6
fields from the Pangeo Google-Cloud Zarr catalogue (no bulk download). It is written for
Google Colab (it will prompt to upload the reef-coordinates file) but runs locally if
you set `IN_COORDS` at the top to your equal-area reef-coordinates file.
Data & code availability
Code: this repository, archived at a DOI-minting service (e.g. Zenodo) on
acceptance.
Bleaching observations: the event-resolved (SAF) and gridded (GLO) source datasets
are from the published sources cited in the manuscript. The procedure used to reduce the
raw gridded database to one ordered severity value per 5 × 5 km reef site (the 0.25° grid
is used only to attach the wind and SST predictors) is documented in
`docs/glo_dataset_construction.md`. The derived
analysis tables used by the scripts are deposited in the data repository cited in the
manuscript, so all results can be reproduced without re-running the (frequently updated)
source database. GLO inputs: `01`/`06` read the derived GLO table
`glo_bleaching_variables_PCA.xlsx`, which holds the same variables as sheet
"Table S2 – GLO" of the published `Supplementary_Data_S1_S2.xlsx` read by `11`
(a flat single-sheet export of the same data); both are provided in the deposit.
Environmental fields: tropical-cyclone tracks (IBTrACS v04r01), NOAA blended sea-surface
winds, SST climatologies, aragonite saturation, and CMIP6 `sfcWind` are public; provenance
(including the exact CMIP6 models/members/grids used) is written by `05_…` to the
`models_used` sheet of its output workbook. Because `05_…` streams the live Pangeo
catalogue, a re-run may select a slightly different model set over time; the deposited
`models_used` sheet is the authoritative record of the ensemble that produced the
published numbers.
Citation
Lapenis, A. Wind effects on coral bleaching severity. (code) — see the
manuscript for the full reference and the archived DOI.
