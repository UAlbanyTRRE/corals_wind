# Construction of the gridded (GLO) bleaching dataset

The raw global coral-bleaching database is updated continuously, so we do not distribute
the extraction script (its output is not stable over time). Instead we document here the
fixed, reproducible procedure used to reduce the raw records to the analysis dataset, so
that the selection rules — not a snapshot of a moving database — define the dataset.

## Procedure

1. **Grid projection.** The global coral-reef polygon layer was projected onto a regular
   0.25° × 0.25° latitude–longitude grid, matching the resolution of the SST and
   surface-wind fields used as predictors. Each reef-containing grid cell defines one
   spatial unit of analysis.

2. **One record per cell.** Where a cell contained more than one bleaching observation,
   the **earliest (oldest) record** of bleaching severity in that cell was retained, giving
   a single observation per cell.

3. **Severity binning.** Heterogeneous bleaching-severity and bleaching-prevalence reports
   were harmonised into three ordered classes — **light**, **medium**, and **severe** — and
   assigned the representative values **5, 30, and 75** (approximate mean percent of the
   reef bleached in each class). Categorical reports were mapped by their text code
   (Mild → 5, Moderate → 30, Severe → 75); quantitative percent-bleached reports were binned
   **0.1–10% → 5 (light)**, **11–50% → 30 (medium)**, **>50% → 75 (severe)**, taking the text
   code where present and otherwise the percent bin. These values order the classes for the
   ordinal model; the analysis uses the ordered categories, not the magnitudes themselves.
   (Implemented in `R/09_severity_binning_gcbd.R`.)

## Notes for the Methods section (recommended additions)

Before this goes into the manuscript, two short clarifications will pre-empt obvious
reviewer questions (the class thresholds above are now explicit):

- **Justify "oldest record".** One clause on *why* the earliest record is kept (e.g. to
  obtain a single, non-duplicated observation per cell and to anchor on the first
  documented event), so the choice does not read as arbitrary. If results are robust to
  the alternative (most-severe, or random record per cell), say so.
- **Clarify the 5/30/75 values.** State explicitly that these are representative
  midpoints used only to order the three classes for ordinal regression, so a reader does
  not mistake them for a continuous response.

The 0.25° analysis grid here is distinct from the ~5 km equal-area reef polygons used for
the CMIP6 projection; keep the two griddings clearly separated in the Methods.
