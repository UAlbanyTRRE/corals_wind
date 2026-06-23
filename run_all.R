# =============================================================================
# run_all.R — master run-order script
# Wind effects on coral bleaching severity (Lapenis)
#
# Sources the R analysis scripts in the intended pipeline order:
#   predictor / response construction  ->  statistical analysis  ->  (projection)
#
# Run from the repository root, e.g.:
#   setwd("/path/to/coral-wind-bleaching")
#   source("run_all.R")
#
# Requirements before running:
#   * Place the input files named in each script header into data/input/.
#   * Install the packages in r_requirements.txt (the scripts no longer
#     auto-install; they stop with a message if a package is missing).
#
# Notes:
#   * 08_blended_wind_monthly_extraction.R is a HELPER MODULE: sourcing it only
#     defines functions (interpolate_windspeed, extract_monthly_mean_wind). It is
#     sourced here for availability but performs no extraction until you call
#     extract_monthly_mean_wind(...) with a loaded blended-wind field.
#   * 05_cmip6_wind_delU_per_reef.py is PYTHON (CMIP6 projection) and is NOT run
#     from R. Run it separately:  python 05_cmip6_wind_delU_per_reef.py
# =============================================================================

# Fail early if obviously not at the repository root.
if (!file.exists("01_glo_olr_wind_turbidity.R"))
  stop("Run from the repository root (the folder containing the analysis scripts and data/).",
       call. = FALSE)

dir.create("data/output", showWarnings = FALSE, recursive = TRUE)

run_order <- c(
  # --- predictor & response construction ---
  "07_tcp_pdi_holland.R",                  # tropical-cyclone power per site
  "09_severity_binning_gcbd.R",            # 3-class GLO severity response
  "03_omega_extraction.R",                 # aragonite saturation onto reefs
  "08_blended_wind_monthly_extraction.R",  # helper module (defines functions only)
  # --- statistical analysis ---
  "01_glo_olr_wind_turbidity.R",           # GLO ordinal regression (FE & ME)
  "06_glo_shapley_mundlak.R",              # GLO Shapley + Mundlak (~18% wind share)
  "02_saf_mediation_residual_wind.R",      # SAF FE+ME, residual wind, mediation
  "10_saf_shapley.R",                      # SAF Shapley (reconstruction; see header)
  "04_bud_box_model.R",                    # BUD carbonate box model
  "11_glo_selection_robustness.R"          # GLO wind-OR robustness
)

for (f in run_order) {
  message("\n==== Running: ", f, " ====")
  source(f, echo = FALSE)
}

message("\nAll R scripts completed. For the CMIP6 projection, run:\n",
        "  python 05_cmip6_wind_delU_per_reef.py")
