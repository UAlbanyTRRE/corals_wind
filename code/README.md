# Wind effects on coral bleaching severity — analysis code

Code accompanying:

**Lapenas, A. ***Wind effects on coral bleaching severity* (submitted).

This repository contains the analysis and modelling code used to evaluate relationships between sustained background wind speed, tropical-cyclone power, ocean carbonate chemistry, and coral bleaching severity. The repository also includes scripts used to project future changes in reef wind regimes using CMIP6 climate-model ensembles.

## Repository structure

```text
R/                 Statistical analyses and data-processing scripts
python/            CMIP6 projection analyses
data/input/        Input datasets
data/output/       Generated outputs
docs/              Supporting documentation
run_all.R          Master workflow script
```

## Main analyses

| Analysis                            | Script                                         |
| ----------------------------------- | ---------------------------------------------- |
| Tropical-cyclone power calculation  | `07_tcp_pdi_holland.R`                         |
| Wind extraction and processing      | `08_blended_wind_monthly_extraction.R`         |
| GLO bleaching severity construction | `09_severity_binning_gcbd.R`                   |
| Ordinal-regression models           | `01_glo_olr_wind_turbidity.R`                  |
| Shapley decomposition               | `06_glo_shapley_mundlak.R`, `10_saf_shapley.R` |
| Mediation analyses                  | `02_saf_mediation_residual_wind.R`             |
| BUD carbonate model                 | `04_bud_box_model.R`                           |
| CMIP6 wind projections              | `05_cmip6_wind_delU_per_reef.py`               |

## Reproducibility

All scripts use paths relative to the repository root. Input files should be placed in `data/input/`, and outputs are written to `data/output/`.

The recommended workflow is:

```r
source("run_all.R")
```

from the repository root directory.

## Data availability

The original bleaching observations and environmental datasets are obtained from the public sources cited in the manuscript.

Because some source databases are continuously updated, the derived analysis tables used in the manuscript are archived separately with the publication. The procedure used to construct the gridded bleaching dataset is documented in:

`docs/glo_dataset_construction.md`

## Code availability

The archived version of this repository associated with the publication will be deposited in Zenodo and assigned a DOI.

## Citation

If you use this code, please cite:

Lapenas, A. *Wind effects on coral bleaching severity*.
Analysis scripts
