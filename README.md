# Wind effects on coral bleaching severity — analysis code

Code accompanying **Lapenis** (*Nature Communications*, submitted): analysis of how
sustained background mean wind speed relates to coral-bleaching severity, the supporting
ocean-carbonate (Box Upwelling–Diffusion) model, and the CMIP6 projection of future
mean-wind change over reefs.

This repository contains analysis and modelling code only. Input datasets are obtained
from the public sources listed in the manuscript's Data Availability statement; the
derived analysis tables are provided as Supplementary Data with the manuscript and in the
data deposit (see **Data & code availability**). No local file paths are embedded — every
script reads from `data/input/` and writes to `data/output/` relative to the repository
root.

Archived on Zenodo: concept DOI [10.5281/zenodo.20720763](https://doi.org/10.5281/zenodo.20720763)
(always resolves to the latest version). The reef-scale CMIP6 wind-change dataset is
archived separately at [10.5281/zenodo.20721850](https://doi.org/10.5281/zenodo.20721850).

## Repository structure

```
corals_wind/
# --- predictor & response construction ---
├── 07_tcp_pdi_holland.R               Tropical-cyclone power (TCP/PDI) within 400 km,
│                                      Holland radial-wind profile, IBTrACS v04r01
├── 08_blended_wind_monthly_extraction.R  Monthly mean wind (1993–2020) to reef cells
│                                      (NOAA blended winds) — helper module (functions only)
├── 09_severity_binning_gcbd.R         GCBD bleaching severity → 3 ordered classes (5/30/75)
├── 03_omega_extraction.R              Aragonite saturation (Ω) at 0/50/100 m onto reefs
# --- statistical analysis ---
├── 01_glo_olr_wind_turbidity.R        GLO ordinal regression (FE & ME), wind × turbidity
├── 06_glo_shapley_mundlak.R           GLO Shapley R² decomposition + Mundlak within-region
├── 02_saf_mediation_residual_wind.R   SAF FE+ME, residual-wind test, mediation
├── 10_saf_shapley.R                   SAF Shapley decomposition
├── 04_bud_box_model.R                 Box Upwelling–Diffusion carbonate model
├── 11_glo_selection_robustness.R      Wind-OR robustness: co-cell cluster + one-per-cell
# --- projection ---
├── 05_cmip6_wind_delU_per_reef.py     CMIP6 ΔU per reef + ensemble + odds maps
├── data/{input,output}/               (place inputs here; outputs land here)
├── docs/glo_dataset_construction.md   How the gridded (GLO) dataset was built
├── run_all.R          master run-order script (sources the scripts in order)
├── requirements.txt (Python)   r_requirements.txt (R)   LICENSE   .gitignore
```

(File-number prefixes are stable IDs, not run order; the pipeline order is: predictor/response
construction → statistical analysis → projection.)

> One file needs a reader's attention before citing: `10_saf_shapley.R` reconstructs the
> SAF decomposition by adapting the GLO routine. It has been run against the deposited SAF
> dataset and verified to reproduce the reported wind share (~18%); the script prints the
> reproduced share and flags whether it matches, so a reader can confirm it independently.
> `08_blended_wind_monthly_extraction.R` is a helper module — sourcing it only defines
> `interpolate_windspeed()` and `extract_monthly_mean_wind()`, which you call with a NOAA
> blended-wind field you have loaded.

## How the scripts map to the manuscript

Script	Produces	Manuscript location
`07_tcp_pdi_holland.R`	Cumulative TCP/PDI per site (6/12-mo, 1993–2020)	Supplementary Note 2 (TCP)
`08_blended_wind_monthly_extraction.R` (helper module)	Monthly mean wind per reef site	Supplementary Note 1 (wind metrics)
`09_severity_binning_gcbd.R`	3-class GLO severity response (5/30/75)	Methods / `docs/glo_dataset_construction.md`
`03_omega_extraction.R`	Ω at 0/50/100 m for each reef cell	Supplementary Note 3; Supplementary Figure 1
`01_glo_olr_wind_turbidity.R`	GLO ordinal-regression coefficients, FE/ME	Main GLO results figure
`06_glo_shapley_mundlak.R`	GLO Shapley shares (the ~18% wind share) + Mundlak	Results; Supplementary Note 4
`02_saf_mediation_residual_wind.R`	SAF model-improvement, residual-wind, mediation	Main SAF results figure
`10_saf_shapley.R`	SAF Shapley shares (Thermal/Depth/Wind)	Supplementary Note 4
`04_bud_box_model.R`	BUD Ω/pCO₂/pH vs wind across (m,n,q)	Supplementary Note 3; Supplementary Figure 1; Supplementary Table 1
`11_glo_selection_robustness.R`	Wind OR: cluster-bootstrap + one-per-cell; 6/12-mo/long-term window comparison	Supplementary Note 7; Supplementary Table 3
`05_cmip6_wind_delU_per_reef.py`	Per-reef ΔU, ensemble means, odds maps	Supplementary Note 6; main-text Fig. 6; Supplementary Figure 2

## Running

R (≥ 4.2). Install the required packages (listed in `r_requirements.txt`), then run
from the repository root, e.g.

```r
setwd("/path/to/corals_wind")
source("01_glo_olr_wind_turbidity.R")
```

To run the full pipeline in the intended order, use the master script
`source("run_all.R")` from the repository root. The scripts do not auto-install packages:
if a required package is missing they stop with an informative message, so install the
dependencies first (or use `renv::restore()`).

Each script begins with a CONFIG / USER-SETTINGS block; place the corresponding input file
in `data/input/` (filenames are documented in the header) and outputs are written to
`data/output/`.

Python (≥ 3.9): `pip install -r requirements.txt`. The CMIP6 script streams CMIP6 fields
from the Pangeo Google-Cloud Zarr catalogue (no bulk download). It is written for Google
Colab (it will prompt to upload the reef-coordinates file) but runs locally if you set
`IN_COORDS` at the top to your equal-area reef-coordinates file.

### Demo inputs (the published Supplementary Data)

The two analysis datasets are provided with the manuscript as **Supplementary Data**
(`Supplementary_Data_S1_S2.xlsx`): sheet **Table S1 – SAF** and sheet **Table S2 – GLO**.
All GLO scripts (`01_glo_olr_wind_turbidity.R`, `06_glo_shapley_mundlak.R` and
`11_glo_selection_robustness.R`) read the same deposited file — sheet `Table S2 - GLO`,
with the header on the second row (`skip = 1`). Place `Supplementary_Data_S1_S2.xlsx` in
`data/input/` and the GLO analyses reproduce from the deposited table directly; no separate
or renamed input file is required.

## Tested with

The published numbers were produced with the environment below.

**Operating system:** Windows 11 (build 26100), x86_64.

**R 4.6.0** — package versions used:
`ordinal` 2025.12-29, `readxl` 1.5.0, `writexl` 1.5.4, `openxlsx` 4.2.8.1,
`seacarb` 3.3.4 (with `SolveSAPHE` 2.1.0, `oce` 1.8-3, `gsw` 1.2-0),
`dplyr` 1.2.1, `tidyr` 1.3.2, `purrr` 1.2.2, `stringr` 1.6.0, `forcats` 1.0.1,
`tibble` 3.3.1, `ggplot2` 4.0.3, `janitor` 2.2.1, `lubridate` 1.9.5, `MASS` 7.3-65.
The following are also required by individual scripts and should be pinned from the
machine on which those scripts were run: `lme4` (02), `data.table` and `geosphere` (07),
`ncdf4` (03), plus `readr` (02, 10) and `pbapply` (07).

**Python ≥ 3.9** (the CMIP6 step in `05_cmip6_wind_delU_per_reef.py` is written for and was
run in Google Colab). It requires `xarray`, `zarr`, `gcsfs`, `cftime`, `pandas`, `numpy`,
`matplotlib`, `openpyxl` and `pyarrow` (`cartopy` optional — maps fall back to no coastlines
if absent), installed via `requirements.txt`. Because the script streams the live Pangeo
CMIP6 catalogue and was run in a managed Colab environment, exact package versions were not
pinned; the deposited `models_used` sheet is the authoritative record of the ensemble that
produced the published numbers.

**Non-standard hardware:** none required for the statistical analyses or the demo (a normal
desktop suffices). The raw CMIP6 read in `05_cmip6_wind_delU_per_reef.py` streams ~84 model
fields from the Pangeo Google-Cloud catalogue and benefits from a fast connection and ample
memory; its output (per-reef ΔU) is deposited at
[10.5281/zenodo.20721850](https://doi.org/10.5281/zenodo.20721850), so reviewers need not
repeat it to reproduce the maps.

**Approximate run-times (normal desktop):** GLO ordinal ensemble (`01`) and Shapley/Mundlak
(`06`) are the heaviest R steps — the GLO ordinal ensemble (thousands of candidate OLR
models) takes ≈ 40 min on a laptop; SAF and BUD scripts run in well under a minute; the
CMIP6 odds-map step runs in a few minutes once the per-reef ΔU table is available. Typical package install time on a normal desktop is a few minutes.

## Reproducing the headline results

| Result | Script | Expected value |
|---|---|---|
| GLO 6-month sustained-wind odds ratio | `01`, `06`, `11` | ≈ 0.70 per s.d. (95% CI ≈ 0.61–0.80) |
| Within-realm wind odds ratio (FE-within / Mundlak) | `06` | ≈ 0.66 / ≈ 0.70 per s.d. |
| SAF wind odds ratio | `02`, `10` | ≈ 0.45 per s.d. |
| GLO wind Shapley share (pooled / within-realm / within-basin) | `06` | ≈ 17.7% / 17.9% / 17.9% |
| SAF wind Shapley share (Thermal / Depth / Wind groups) | `10` | ≈ 25.4% (Thermal 67.1%, Depth 7.6%) |
| Full-model McFadden R² (GLO) | `01`, `06` | ≈ 0.08 |
| BUD modelled Ω_arag–wind slopes (m=n=1 / m=n=2, q=0) | `04` | ≈ +0.023 / +0.036 per m s⁻¹ |
| Reef area with projected wind weakening (SSP2-4.5 / SSP5-8.5) | `05` | ≈ 51% / 59% |

## Data & code availability

**Code:** this repository, archived on Zenodo at
https://doi.org/10.5281/zenodo.20720763 (concept DOI, resolves to the latest version).

**Bleaching observations:** the event-resolved (SAF) and gridded (GLO) source datasets are
from the published sources cited in the manuscript (Safaie et al. 2018; van Woesik &
Kratochwill 2022). The procedure used to reduce the raw gridded database to one ordered
severity value per 5 × 5 km reef site (the 0.25° grid is used only to attach the wind and
SST predictors) is documented in `docs/glo_dataset_construction.md`. All three GLO scripts
read the deposited `Supplementary_Data_S1_S2.xlsx` (sheet "Table S2 – GLO"); the SAF scripts
read sheet "Table S1 – SAF" of the same workbook.

**Environmental fields:** tropical-cyclone tracks (IBTrACS v04r01), NOAA blended sea-surface
winds, SST climatologies, aragonite saturation, and CMIP6 `sfcWind` are public; provenance
(including the exact CMIP6 models/members/grids used) is written by `05_…` to the
`models_used` sheet of its output workbook. Because `05_…` streams the live Pangeo
catalogue, a re-run may select a slightly different model set over time; the deposited
`models_used` sheet is the authoritative record of the ensemble that produced the published
numbers. The reef-scale CMIP6 wind-change dataset is archived at
https://doi.org/10.5281/zenodo.20721850.

## License

MIT License — see `LICENSE`.

## Citation

Lapenis, A. *Wind effects on coral bleaching severity* (code). Zenodo,
https://doi.org/10.5281/zenodo.20720763 — see the manuscript for the full reference.
