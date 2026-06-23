# Construction of the gridded (GLO) bleaching dataset

The raw global coral-bleaching database is updated continuously, so we do not distribute
the extraction script (its output is not stable over time). Instead we document here the
fixed, reproducible procedure used to reduce the raw records to the analysis dataset, so
that the selection rules — not a snapshot of a moving database — define the dataset. The
derived analysis table produced by this procedure is provided as Supplementary Data
(sheet "Table S2 – GLO" of `Supplementary_Data_S1_S2.xlsx`).

## Two distinct griddings

This dataset uses two different grids, kept separate throughout:

- **5 × 5 km equal-area reef grid — the unit of analysis.** Reef survey locations are
  projected onto a 5 × 5 km grid so that retained sites lie at least 5 km apart. One
  record is kept per occupied 5 km cell, yielding the 1,195 site-level observations used
  in the GLO models.
- **0.25° × 0.25° grid — predictor attachment only.** The coarser 0.25° SST and
  surface-wind fields are attached to the 5 km sites by spatial overlay. Because several
  5 km sites can fall within a single 0.25° predictor cell (the 1,195 sites occupy 835
  such cells; 275 of those cells contain ≥ 2 sites), co-located sites may share identical
  gridded predictor values. This shared structure is addressed by the robustness checks in
  Supplementary Note 7 (cluster bootstrap over 0.25° cells; collapsing to one record per
  0.25° cell); it does not define the analysis unit.

(The 5 km analysis grid here is also distinct from the ~5 km equal-area reef polygons used
for the separate CMIP6 wind-change projection; the two are kept clearly separated in the
Methods.)

## Procedure

1. **Site grid.** Reef survey locations were projected onto a 5 × 5 km equal-area grid, so
   that retained sites are at least 5 km apart. Each occupied 5 km cell defines one spatial
   unit of analysis.
2. **One record per site.** Where a 5 km cell contained more than one bleaching
   observation, the **earliest (oldest) record** of bleaching severity in that cell was
   retained, giving a single, non-duplicated observation per site anchored on the first
   documented event. This yields 1,195 site-level observations.
3. **Severity binning.** Heterogeneous bleaching-severity and bleaching-prevalence reports
   were harmonised into three ordered classes — **low**, **moderate**, and **severe** — and
   labelled **5, 30, and 75** (approximate mid-percent of the reef bleached in each class).
   Categorical reports were mapped by their text code (low → 5, moderate → 30, severe → 75);
   quantitative percent-bleached reports were binned **0.1–10% → 5 (low)**,
   **11–50% → 30 (moderate)**, **>50% → 75 (severe)**, taking the text code where present and
   otherwise the percent bin. These labels serve only to **order** the three classes for the
   ordinal model; the analysis uses the ordered categories, not the magnitudes themselves.
   (Implemented in `09_severity_binning_gcbd.R`.)
4. **Predictor attachment.** Gridded SST and surface-wind fields (0.25°) were attached to
   each 5 km site by spatial overlay (see "Two distinct griddings" above).

The derived table from this procedure is read by the GLO analysis scripts
(`01_glo_olr_wind_turbidity.R`, `06_glo_shapley_mundlak.R`,
`11_glo_selection_robustness.R`) from sheet "Table S2 – GLO" of
`Supplementary_Data_S1_S2.xlsx`.
