# =============================================================================
# GLO robustness: co-cell non-independence and record aggregation
# Part of: Wind effects on coral bleaching severity (Lapenis & Jiang)
# Language: R
# Inputs : data/input/Supplementary_Data_S1_S2.xlsx  (sheet "Table S2 - GLO")
# Outputs: data/output/glo_robustness.xlsx
# Depends: readxl, dplyr, ordinal, writexl
#
# WHAT THIS DOES
#   The GLO analysis keeps one bleaching record per 5x5 km reef SITE (sites are
#   >= 5 km apart); the 0.25-deg grid is used only to attach the gridded SST/wind
#   predictors. Because several 5 km sites can fall in one 0.25-deg wind/SST cell,
#   co-cell records share identical gridded predictors. This script shows the
#   wind result is robust to that non-independence by:
#     (1) refitting the published pooled OLR (naive SEs) as the reference;
#     (2) a cluster bootstrap that resamples whole 0.25-deg cells; and
#     (3) collapsing to one record per 0.25-deg cell (earliest) and refitting.
#   The 0.25-deg cell is defined by ROUND-to-nearest (not floor): this is the
#   alignment under which the long-term wind is exactly constant within a cell.
#
# NOTE on record-selection (acclimatization) robustness:
#   Testing oldest vs most-severe vs random EVENT per site requires the 6-month
#   antecedent wind (wind_mean_6m) for the non-earliest events, which is only
#   computed for the first event per site. Substituting the long-term mean wind
#   is NOT a valid stand-in, because the protective effect is specific to the
#   6-month window (6-mo OR ~0.70, p~2e-7; long-term OR ~1.08, n.s.). That check
#   therefore needs antecedent winds derived for the alternative events.
# =============================================================================
pkgs <- c("readxl", "dplyr", "ordinal", "writexl")
missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0)
  stop("Missing required package(s): ", paste(missing_pkgs, collapse = ", "),
       "\n  Install them, or restore the recorded environment with renv::restore().",
       call. = FALSE)
invisible(lapply(pkgs, require, character.only = TRUE))

## ---- CONFIG ----
# Published Supplementary Data workbook. Sheet "Table S2 - GLO" holds the same
# GLO variables that 01/06 read from the derived extract glo_bleaching_variables_PCA.xlsx.
GLO_FILE <- "data/input/Supplementary_Data_S1_S2.xlsx"
GLO_SHEET <- "Table S2 - GLO"
GLO_SKIP  <- 1                       # header sits on the 2nd row of that sheet
OUT_FILE  <- "data/output/glo_robustness.xlsx"
GRID_DEG  <- 0.25
WIND_TERM <- "wind_mean_6m"          # the published GLO wind predictor
RHS <- c("tsa_dhw", "sst_pc1", "sst_pc2", "sst_pc3", WIND_TERM,
         "tcpower_1993_2020_400km", "distance_to_shore", "exposure")
N_BOOT <- 250
set.seed(1)

## ---- load & 0.25-deg cell (round-to-nearest) ----
g <- as.data.frame(readxl::read_excel(GLO_FILE, sheet = GLO_SHEET, skip = GLO_SKIP))
g <- g[!is.na(g$latitude) & !is.na(g$longitude), ]
g$cell <- paste(round(g$latitude / GRID_DEG), round(g$longitude / GRID_DEG), sep = "_")

## ---- helper: fit OLR, return wind OR (per s.d.) ----
wind_or <- function(d) {
  d <- d[stats::complete.cases(d[, c("bleaching_categorical", RHS)]), ]
  for (v in RHS) d[[v]] <- as.numeric(scale(d[[v]]))
  d$y <- factor(d$bleaching_categorical,
                levels = sort(unique(d$bleaching_categorical)), ordered = TRUE)
  m <- ordinal::clm(stats::as.formula(paste("y ~", paste(RHS, collapse = "+"))),
                    data = d, link = "logit")
  s <- summary(m)$coefficients[WIND_TERM, ]
  b <- s["Estimate"]; se <- s["Std. Error"]
  c(n = nrow(d), OR = exp(b), lo = exp(b - 1.96 * se), hi = exp(b + 1.96 * se))
}

## ---- (1) published pooled OLR (naive SEs) ----
ref <- wind_or(g)

## ---- (2) cluster bootstrap by 0.25-deg cell ----
cells <- unique(g$cell); ors <- numeric(0)
for (b in seq_len(N_BOOT)) {
  samp <- sample(cells, length(cells), replace = TRUE)
  db <- do.call(rbind, lapply(samp, function(c) g[g$cell == c, ]))
  ors <- c(ors, tryCatch(wind_or(db)["OR"], error = function(e) NA_real_))
}
ors <- ors[is.finite(ors)]
clus <- c(OR = unname(ref["OR"]), lo = quantile(ors, .025), hi = quantile(ors, .975))

## ---- (3) one record per 0.25-deg cell (earliest) ----
g$.d <- if ("EventDate" %in% names(g)) as.numeric(as.Date(g$EventDate)) else seq_len(nrow(g))
one <- g[order(g$.d), ]; one <- one[!duplicated(one$cell), ]
col <- wind_or(one)

## ---- assemble table ----
tab <- data.frame(
  Variant = c("Pooled OLR, all sites (published)",
              sprintf("Cluster bootstrap by 0.25-deg cell (%d reps)", length(ors)),
              "One record per 0.25-deg cell (earliest)"),
  n       = c(ref["n"], ref["n"], col["n"]),
  Wind_OR = round(c(ref["OR"], clus["OR"], col["OR"]), 3),
  CI_low  = round(c(ref["lo"], clus["lo"], col["lo"]), 3),
  CI_high = round(c(ref["hi"], clus["hi"], col["hi"]), 3),
  row.names = NULL)
cat("\n--- GLO wind OR robustness (predictor:", WIND_TERM, ") ---\n"); print(tab)

## ---- (optional) wind-window comparison: 6m vs 12m vs long-term ----
window_compare <- function() {
  base <- setdiff(RHS, WIND_TERM)
  do.call(rbind, lapply(c("wind_mean_6m", "wind_mean_12m", "wind_mean_1993_2020"), function(w) {
    if (!w %in% names(g)) return(NULL)
    d <- g[stats::complete.cases(g[, c("bleaching_categorical", base, w)]), ]
    P <- c(base, w); for (v in P) d[[v]] <- as.numeric(scale(d[[v]]))
    d$y <- factor(d$bleaching_categorical, levels = sort(unique(d$bleaching_categorical)), ordered = TRUE)
    m <- ordinal::clm(stats::as.formula(paste("y ~", paste(P, collapse = "+"))), data = d, link = "logit")
    s <- summary(m)$coefficients[w, ]
    data.frame(window = w, OR = round(exp(s["Estimate"]), 3),
               CI_low = round(exp(s["Estimate"] - 1.96 * s["Std. Error"]), 3),
               CI_high = round(exp(s["Estimate"] + 1.96 * s["Std. Error"]), 3),
               p = signif(s["Pr(>|z|)"], 3), row.names = NULL)
  }))
}
wc <- window_compare(); cat("\n--- wind-window comparison ---\n"); print(wc)

writexl::write_xlsx(list(robustness = tab, wind_windows = wc), OUT_FILE)
cat(sprintf("\nWrote: %s\n", OUT_FILE))
