# Wind effects on coral bleaching severity — analysis code

Code accompanying:

**Lapenas, A.** *Wind effects on coral bleaching severity* (submitted).

This repository contains the analysis and modelling code used to evaluate relationships between sustained background wind speed, tropical-cyclone power, ocean carbonate chemistry, and coral bleaching severity. The repository also includes scripts used to project future changes in reef wind regimes using CMIP6 climate-model ensembles.

## Repository structure

```text
corals_wind/
│
├── code/
│   ├── README.md
│   ├── 01_glo_olr_wind_turbidity.R
│   ├── 02_saf_mediation_residual_wind.R
│   ├── 03_omega_extraction.R
│   ├── 04_bud_box_model.R
│   ├── 05_cmip6_wind_delU_per_reef.py
│   ├── 06_glo_shapley_mundlak.R
│   ├── 07_tcp_pdi_holland.R
│   ├── 08_blended_wind_monthly_extraction.R
│   ├── 09_severity_binning_gcbd.R
│   ├── 10_saf_shapley.R
│   └── 11_glo_selection_robustness.R
│
├── metadata/
│   └── GLO_dataset_construction.md
│
├── LICENSE
└── README.md
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
| Selection and robustness analyses   | `11_glo_selection_robustness.R`                |

## Reproducibility

All scripts are intended to be run from the repository root. Input datasets should be obtained from the sources listed in the manuscript and placed in locations specified within individual script headers.

The scripts are organized into three stages:

1. Construction of predictor and response variables.
2. Statistical analyses of the SAF and GLO bleaching datasets.
3. CMIP6 projections and BUD-model simulations.

## Data availability

The original bleaching observations and environmental datasets are obtained from the public sources cited in the manuscript.

Because some source databases are continuously updated, the derived analysis tables used in the manuscript are archived separately with the publication. The procedure used to construct the gridded bleaching dataset is documented in:

`metadata/GLO_dataset_construction.md`

## Code availability

The archived version of this repository associated with the publication will be deposited in Zenodo and assigned a DOI upon publication.

## Citation

If you use this code, please cite:

Lapenas, A. *Wind effects on coral bleaching severity*.
