# =============================================================================
# SAF Shapley decomposition of McFadden R-squared
# Thermal / Depth / Wind predictor groups
#
# Part of : Wind effects on coral bleaching severity (Lapenis & Jiang)
# Language: R
# Input   : data/input/saf_synced_bleaching_wind_TC_zscored.csv
# Output  : data/output/saf_shapley_results.xlsx
# Depends : readr, ordinal, writexl
#
# Purpose
# -------
# This script reproduces the SAF group-Shapley decomposition reported in
# Supplementary Note 4. It fits cumulative-link ordinal-logit models and
# decomposes McFadden R-squared among three predictor groups:
#   1. Thermal
#   2. Depth
#   3. Wind
#
# Reproducibility check
# ---------------------
# The script prints the reproduced Wind share of McFadden R-squared at the end.
# The expected value is approximately 18%. The script flags the run as PASS if
# the reproduced value falls within EXPECTED_TOL percentage points of 18%.
#
# Run from the repository root.
# =============================================================================

## ---- packages ----
pkgs <- c("readr", "ordinal", "writexl")
missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    "Missing required package(s): ", paste(missing_pkgs, collapse = ", "),
    "\nInstall them with install.packages() or restore the recorded environment with renv::restore().",
    call. = FALSE
  )
}
invisible(lapply(pkgs, require, character.only = TRUE))

## ---- configuration ----
SAF_FILE <- "data/input/saf_synced_bleaching_wind_TC_zscored.csv"
OUT_FILE <- "data/output/saf_shapley_results.xlsx"
OUTCOME  <- "Y_nanless"          # ordinal bleaching severity, expected 0..4
EXPECTED_WIND_PCT <- 18           # reported approximate Wind share, percent
EXPECTED_TOL      <- 3            # pass if within +/- 3 percentage points

## ---- small utilities ----
pick <- function(nm, cand, label) {
  hit <- cand[cand %in% nm]
  if (length(hit) == 0) {
    stop(
      "Could not resolve column for ", label, ". Tried: ",
      paste(cand, collapse = ", "),
      call. = FALSE
    )
  }
  hit[1]
}

safe_dir_create <- function(path) {
  out_dir <- dirname(path)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

fit_clm_ll <- function(rhs, data) {
  form <- stats::as.formula(paste("y ~", rhs))
  fit <- tryCatch(
    ordinal::clm(form, data = data, link = "logit"),
    error = function(e) {
      message("Model failed for RHS: ", rhs)
      message("Error: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(fit)) NA_real_ else as.numeric(stats::logLik(fit))
}

mcfadden_r2 <- function(rhs, data, ll_null) {
  ll <- fit_clm_ll(rhs, data)
  if (!is.finite(ll)) return(NA_real_)
  1 - ll / ll_null
}

make_rhs <- function(terms) {
  terms <- unique(terms[!is.na(terms) & nzchar(terms)])
  if (length(terms) == 0) "1" else paste(terms, collapse = " + ")
}

subset_key <- function(S) paste(sort(S), collapse = "|")

shapley_decompose <- function(groups, base_terms, data, ll_null) {
  keys <- names(groups)
  n_groups <- length(keys)

  # All subsets of predictor-group names.
  subsets <- unlist(
    lapply(0:n_groups, function(r) utils::combn(keys, r, simplify = FALSE)),
    recursive = FALSE
  )

  r2_by_subset <- setNames(numeric(length(subsets)), vapply(subsets, subset_key, character(1)))

  for (S in subsets) {
    rhs_terms <- c(base_terms, unlist(groups[S], use.names = FALSE))
    rhs <- make_rhs(rhs_terms)
    r2_by_subset[subset_key(S)] <- mcfadden_r2(rhs, data, ll_null)
  }

  if (any(!is.finite(r2_by_subset))) {
    bad <- names(r2_by_subset)[!is.finite(r2_by_subset)]
    stop("At least one subset model failed: ", paste(bad, collapse = ", "), call. = FALSE)
  }

  phi <- setNames(numeric(n_groups), keys)

  for (gn in keys) {
    others <- setdiff(keys, gn)
    for (r in 0:(n_groups - 1)) {
      for (S in utils::combn(others, r, simplify = FALSE)) {
        weight <- factorial(length(S)) * factorial(n_groups - length(S) - 1) / factorial(n_groups)
        phi[gn] <- phi[gn] + weight *
          (r2_by_subset[subset_key(c(S, gn))] - r2_by_subset[subset_key(S)])
      }
    }
  }

  full_key <- subset_key(keys)
  full_r2 <- r2_by_subset[full_key]

  data.frame(
    group = keys,
    R2_contribution = as.numeric(phi),
    pct_of_full_R2 = 100 * as.numeric(phi) / as.numeric(full_r2),
    full_McFadden_R2 = as.numeric(full_r2),
    row.names = NULL,
    check.names = FALSE
  )
}

## ---- read data ----
if (!file.exists(SAF_FILE)) {
  stop(
    "Input file not found: ", SAF_FILE,
    "\nRun this script from the repository root, or update SAF_FILE in the CONFIG block.",
    call. = FALSE
  )
}

d <- as.data.frame(readr::read_csv(SAF_FILE, show_col_types = FALSE))
nm <- names(d)

if (!(OUTCOME %in% nm)) {
  stop("Outcome column not found: ", OUTCOME, call. = FALSE)
}

## ---- resilient column resolution ----
# Candidate names mirror earlier SAF scripts and common sanitized variants.
DEPTH <- pick(nm, c("depths", "Depth", "RAW_depths"), "Depth")
ROTC  <- pick(nm, c("ROTC_.SS.", "ROTCSS", "RAW_ROTC_.SS."), "ROTC")
AC1   <- pick(nm, c("Acute1...10", "Acute1...59", "RAW_Acute1"), "Acute1")
DHW30 <- pick(nm, c("DHW_.l30.", "DHW30", "RAW_DHW_.l30."), "DHW30")
DTR30 <- pick(nm, c("DTR_.30.", "DTR30", "RAW_DTR_.30."), "DTR30")
TT    <- pick(nm, c("TT...12", "TT...84", "RAW_TT"), "TT")
WIND  <- pick(nm, c("Wind_12m_mean_z", "Wind_12m_mean_z."), "Wind")

resolved_columns <- data.frame(
  role = c("Outcome", "Depth", "ROTC", "Acute1", "DHW30", "DTR30", "TT", "Wind"),
  column = c(OUTCOME, DEPTH, ROTC, AC1, DHW30, DTR30, TT, WIND),
  stringsAsFactors = FALSE
)

## ---- predictor groups ----
GROUPS <- list(
  Thermal = c(DHW30, DTR30, ROTC, AC1, TT),
  Depth   = c(DEPTH),
  Wind    = c(WIND)
)
ALLP <- unique(unlist(GROUPS, use.names = FALSE))

## ---- prepare response and complete-case data ----
d$y <- factor(d[[OUTCOME]], levels = sort(unique(d[[OUTCOME]])), ordered = TRUE)

analysis_cols <- c("y", ALLP)
d0_n <- nrow(d)
d <- d[stats::complete.cases(d[, analysis_cols]), , drop = FALSE]

if (nrow(d) == 0) stop("No complete cases remain after filtering.", call. = FALSE)
if (length(unique(d$y)) < 2) stop("Outcome has fewer than two levels after filtering.", call. = FALSE)

# Re-scale all predictors in the analysis dataset. This is harmless if columns
# are already z-scored; McFadden R2 and Shapley shares are invariant to linear
# rescaling of predictors in these models.
for (cc in ALLP) d[[cc]] <- as.numeric(scale(d[[cc]]))

cat("\n--- SAF Shapley decomposition setup ---\n")
cat(sprintf("Input file: %s\n", SAF_FILE))
cat(sprintf("Rows read: %d | complete cases used: %d | rows removed: %d\n", d0_n, nrow(d), d0_n - nrow(d)))
cat("Resolved columns:\n")
print(resolved_columns, row.names = FALSE)
cat("Predictor groups:\n")
print(data.frame(
  group = names(GROUPS),
  n_predictors = vapply(GROUPS, length, integer(1)),
  predictors = vapply(GROUPS, paste, character(1), collapse = ", "),
  row.names = NULL
), row.names = FALSE)

## ---- null and full OLR models ----
ll_null <- fit_clm_ll("1", d)
if (!is.finite(ll_null)) stop("Null model failed.", call. = FALSE)

full_rhs <- make_rhs(ALLP)
olr <- ordinal::clm(stats::as.formula(paste("y ~", full_rhs)), data = d, link = "logit")
full_ll <- as.numeric(stats::logLik(olr))
full_mcf <- 1 - full_ll / ll_null

cat("\n--- OLR reference model ---\n")
print(summary(olr)$coefficients)
cat(sprintf("\nNull logLik = %.6f\n", ll_null))
cat(sprintf("Full logLik = %.6f\n", full_ll))
cat(sprintf("Full McFadden R2 = %.6f\n", full_mcf))

## ---- pooled Shapley decomposition ----
sh <- shapley_decompose(GROUPS, base_terms = character(0), data = d, ll_null = ll_null)

cat("\n--- Pooled group-Shapley decomposition ---\n")
print(sh, digits = 4, row.names = FALSE)

## ---- reproducibility check ----
wind_pct <- sh$pct_of_full_R2[sh$group == "Wind"]
if (length(wind_pct) != 1 || !is.finite(wind_pct)) {
  stop("Could not compute a finite Wind percentage.", call. = FALSE)
}

pass <- abs(wind_pct - EXPECTED_WIND_PCT) <= EXPECTED_TOL
status <- if (pass) "PASS" else "CHECK REQUIRED"

check_df <- data.frame(
  metric = "SAF Wind share of McFadden R2",
  reproduced_percent = wind_pct,
  expected_approx_percent = EXPECTED_WIND_PCT,
  tolerance_percent_points = EXPECTED_TOL,
  difference_percent_points = wind_pct - EXPECTED_WIND_PCT,
  status = status,
  stringsAsFactors = FALSE
)

cat("\n--- REPRODUCIBILITY CHECK ---\n")
cat(sprintf(
  "SAF Wind share of McFadden R2 = %.1f%%; reported value is approximately %.0f%%.\n",
  wind_pct, EXPECTED_WIND_PCT
))
if (pass) {
  cat(sprintf(
    "PASS: reproduced Wind share is within +/- %.0f percentage points of the reported value.\n",
    EXPECTED_TOL
  ))
} else {
  cat(sprintf(
    "CHECK REQUIRED: reproduced Wind share differs by %.1f percentage points. Check input data, filtering, and GROUPS.\n",
    wind_pct - EXPECTED_WIND_PCT
  ))
}

## ---- write workbook ----
safe_dir_create(OUT_FILE)

olr_coef <- as.data.frame(summary(olr)$coefficients)
olr_coef$term <- row.names(olr_coef)
row.names(olr_coef) <- NULL
olr_coef <- olr_coef[, c("term", setdiff(names(olr_coef), "term"))]

model_summary <- data.frame(
  n_rows_read = d0_n,
  n_complete_cases_used = nrow(d),
  n_rows_removed = d0_n - nrow(d),
  null_logLik = ll_null,
  full_logLik = full_ll,
  full_McFadden_R2 = full_mcf,
  stringsAsFactors = FALSE
)

writexl::write_xlsx(
  list(
    Model_summary = model_summary,
    Resolved_columns = resolved_columns,
    OLR_coefficients = olr_coef,
    Shapley_pooled = sh,
    Reproducibility_check = check_df
  ),
  OUT_FILE
)

cat(sprintf("\nWrote workbook: %s\n", OUT_FILE))
cat("\nDone.\n")
