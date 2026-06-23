# =============================================================================
# SAF (n=81): FE + ME ordinal models, residual-wind test, approximate mediation
# Part of: Wind effects on coral bleaching severity (Lapenis)
# Language: R
# Inputs : data/input/saf_synced_bleaching_wind_TC_zscored.csv
# Outputs: data/output/
# Depends: readr, dplyr, tibble, stringr, ordinal, openxlsx, lme4
# Notes  : Paths are set in the CONFIG/USER-SETTINGS block below; place input
#          files in data/input/ and run from the repository root.
# =============================================================================

# ============================================================
# FINAL SCRIPT (K4 only): FE + ME OLR + Residual Wind test + Approx. Mediation
# FIXED: reef_group guaranteed to exist in ALL modeling data frames

# ============================================================
# NOTE: rm(list = ls()) was removed here. Clearing the global environment is
# discouraged in shared/deposited code because it silently wipes the caller's
# session when the script is sourced. Run this script in a fresh R session if a
# clean environment is needed.

# --------------------------
# USER SETTINGS

# --------------------------
FILE_PATH <- "data/input/saf_synced_bleaching_wind_TC_zscored.csv"
OUT_DIR   <- "data/output"
RAW_OUTCOME_COL <- "Y_nanless"                # 0..4
REEF_ID_COLNAME_CANDIDATES <- c("Reef","TS_Code","row81")  # tries these in order
OUT_XLSX <- file.path(OUT_DIR, "K4_FE_ME_WindResidual_Mediation_FINAL.xlsx")
# Wind to use (your best wind)
WIND_COL_CANDIDATES <- c("Wind_12m_mean_z","Wind_12m_mean_z.")  # add variants if needed
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# --------------------------
# PACKAGES

# --------------------------
pkgs <- c("readr","dplyr","tibble","stringr","ordinal","openxlsx","lme4")
missing_pkgs <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(missing_pkgs) > 0) {
  stop("Missing required package(s): ", paste(missing_pkgs, collapse = ", "),
       "\n  Install with: install.packages(c(",
       paste(sprintf('"%s"', missing_pkgs), collapse = ", "), "))",
       "\n  or restore the recorded environment with renv::restore().",
       call. = FALSE)
}
invisible(lapply(pkgs, library, character.only=TRUE))

# --------------------------
# openxlsx temp-dir patch (Windows-safe)
#
# WHY: openxlsx builds the .xlsx by zipping files staged in a temporary
# directory. On some Windows setups the default temp path contains spaces or is
# redirected (e.g. OneDrive), which can break the zip step. Pointing the temp
# dir at a local folder under data/output avoids this. Harmless on macOS/Linux.
# --------------------------
SAFE_TMP <- file.path(OUT_DIR, "_tmp_openxlsx")
dir.create(SAFE_TMP, showWarnings = FALSE, recursive = TRUE)
Sys.setenv(TMPDIR = SAFE_TMP, TEMP = SAFE_TMP, TMP = SAFE_TMP)
options(tmpdir = SAFE_TMP)

# --------------------------
# HELPERS

# --------------------------
`%||%` <- function(x,y) if(!is.null(x) && length(x)>0 && !all(is.na(x))) x else y
make_safe_names <- function(df){
  names(df) <- make.names(trimws(names(df)), unique=TRUE)
  df
}
pick_first_present <- function(cols, candidates){
  for (cn in candidates){
    cn1 <- make.names(cn)
    hit <- cols[tolower(cols)==tolower(cn)] %||% cols[tolower(cols)==tolower(cn1)]
    if(length(hit)>0) return(hit[1])
  }
  NA_character_
}
make_formula <- function(response, terms){
  if(length(terms)==0) return(as.formula(paste0(response," ~ 1")))
  as.formula(paste0(response," ~ ", paste(terms, collapse=" + ")))
}
fit_clm_safe <- function(formula, data){
  tryCatch(ordinal::clm(formula, data=data, link="logit"), error=function(e) NULL)
}
fit_clmm_safe <- function(formula, data, group_var){
  f2 <- update(formula, paste0(". ~ . + (1|", group_var, ")"))
  tryCatch(ordinal::clmm(f2, data=data, link="logit", Hess=TRUE), error=function(e) NULL)
}
AICc <- function(model){
  if(is.null(model)) return(NA_real_)
  k <- attr(logLik(model), "df")
  n <- nobs(model)
  aic <- AIC(model)
  if(anyNA(c(k,n)) || (n-k-1)<=0) return(NA_real_)
  aic + 2*k*(k+1)/(n-k-1)
}
mcfadden_r2 <- function(full, null){
  if(is.null(full) || is.null(null)) return(NA_real_)
  1 - as.numeric(logLik(full) / logLik(null))
}
nagelkerke_r2 <- function(full, null){
  if(is.null(full) || is.null(null)) return(NA_real_)
  n <- nobs(full)
  llf <- as.numeric(logLik(full))
  lln <- as.numeric(logLik(null))
  r2_cs <- 1 - exp((2/n)*(lln - llf))
  r2_nk <- r2_cs / (1 - exp((2/n)*lln))
  r2_nk
}
extract_slopes <- function(mod, model_type){
  if(is.null(mod)) {
    return(tibble(term=character(), estimate=numeric(), se=numeric(), z=numeric(), p=numeric(),
                  model_type=model_type))
  }
  sm <- try(summary(mod), silent=TRUE)
  if(inherits(sm,"try-error") || is.null(sm$coefficients)) {
    return(tibble(term=character(), estimate=numeric(), se=numeric(), z=numeric(), p=numeric(),
                  model_type=model_type))
  }
  cf <- sm$coefficients
  n_thres <- length(mod$alpha) %||% 0
  if(n_thres==0) n_thres <- mod$nTheta %||% 0
  if(nrow(cf) <= n_thres) return(tibble(term=character(), estimate=numeric(), se=numeric(), z=numeric(), p=numeric(),
                                        model_type=model_type))
  slopes <- cf[(n_thres+1):nrow(cf), , drop=FALSE]
  tibble(
    term = rownames(slopes),
    estimate = as.numeric(slopes[, "Estimate"]),
    se = as.numeric(slopes[, "Std. Error"]),
    z = as.numeric(slopes[, "z value"]),
    p = as.numeric(slopes[, "Pr(>|z|)"]),
    model_type = model_type
  )
}
add_wald_ci <- function(df, level=0.95){
  zcrit <- qnorm(1 - (1-level)/2)
  df %>%
    mutate(
      ci_low  = estimate - zcrit*se,
      ci_high = estimate + zcrit*se,
      sig = case_when(
        is.na(p) ~ "",
        p < 0.001 ~ "***",
        p < 0.01  ~ "**",
        p < 0.05  ~ "*",
        TRUE ~ ""
      )
    )
}
prep_model_data_K4 <- function(dat, outcome, group_var, preds){
  needed <- unique(c(outcome, group_var, preds))
  needed <- needed[!is.na(needed) & nzchar(needed)]
  stopifnot(group_var %in% names(dat))  # hard guard
  d <- dat[, needed, drop=FALSE]
  d <- d[complete.cases(d), , drop=FALSE]
  if(nrow(d) < 20) return(NULL)
  if(!is.ordered(d[[outcome]]) || length(levels(d[[outcome]])) != 4){
    stop("Outcome must be ordered factor with 4 levels: ", outcome)
  }
  d
}
# Robust probability extraction for "highest category"
p_highest <- function(mod, newdata){
  pr <- try(ordinal::predict(mod, newdata=newdata, type="prob"), silent=TRUE)
  if(inherits(pr,"try-error") || is.null(pr)) return(NA_real_)
  if(is.list(pr) && !is.null(pr$fit)) pr <- pr$fit
  if(is.vector(pr)) return(as.numeric(pr[length(pr)]))
  if(is.matrix(pr) || is.data.frame(pr)) return(as.numeric(pr[, ncol(pr), drop=TRUE]))
  NA_real_
}
deltaP_highest <- function(mod, d, term, terms_all, group_var=NULL){
  base <- d[1, , drop=FALSE]
  for(nm in terms_all){
    v <- d[[nm]]
    base[[nm]] <- if(is.numeric(v)) median(v, na.rm=TRUE) else levels(v)[1]
  }
  if(!is.null(group_var) && group_var %in% names(d)){
    base[[group_var]] <- d[[group_var]][1]
  }
  vterm <- d[[term]]
  if(!is.numeric(vterm)) return(tibble(term=term, deltaP=NA_real_, p_lo=NA_real_, p_hi=NA_real_))
  qlo <- as.numeric(quantile(vterm, 0.05, na.rm=TRUE, type=7))
  qhi <- as.numeric(quantile(vterm, 0.95, na.rm=TRUE, type=7))
  lo <- base; hi <- base
  lo[[term]] <- qlo; hi[[term]] <- qhi
  plo <- p_highest(mod, lo)
  phi <- p_highest(mod, hi)
  tibble(term=term, p_lo=plo, p_hi=phi, deltaP=phi-plo)
}
# Extract one slope term from clm/clmm
extract_one_slope <- function(mod, term, label){
  if(is.null(mod)) return(tibble(which=label, term=term, estimate=NA_real_, se=NA_real_, z=NA_real_, p=NA_real_))
  sm <- summary(mod)$coefficients
  n_thres <- length(mod$alpha) %||% 0
  if(n_thres==0) n_thres <- mod$nTheta %||% 0
  slopes <- sm[(n_thres+1):nrow(sm), , drop=FALSE]
  if(!(term %in% rownames(slopes))) return(tibble(which=label, term=term, estimate=NA_real_, se=NA_real_, z=NA_real_, p=NA_real_))
  tibble(which=label, term=term,
         estimate=as.numeric(slopes[term,"Estimate"]),
         se=as.numeric(slopes[term,"Std. Error"]),
         z=as.numeric(slopes[term,"z value"]),
         p=as.numeric(slopes[term,"Pr(>|z|)"]))
}

# --------------------------
# LOAD DATA + BUILD K4 OUTCOME + GROUPING

# --------------------------
if(!file.exists(FILE_PATH)) stop("File not found: ", FILE_PATH)
dat0 <- readr::read_csv(FILE_PATH, show_col_types=FALSE) |> make_safe_names()
cols <- names(dat0)
RAW_OUT <- pick_first_present(cols, c(RAW_OUTCOME_COL))
if(is.na(RAW_OUT)) stop("Raw outcome column not found: ", RAW_OUTCOME_COL)
# K4 outcome
dat0$BleachResp4 <- dplyr::case_when(
  dat0[[RAW_OUT]] %in% c(0,1) ~ 1,
  dat0[[RAW_OUT]] == 2 ~ 2,
  dat0[[RAW_OUT]] == 3 ~ 3,
  dat0[[RAW_OUT]] == 4 ~ 4,
  TRUE ~ NA_real_
)
dat0$BleachResp4 <- factor(dat0$BleachResp4, levels=1:4, ordered=TRUE)
# Create reef_group robustly
REEF_COL <- pick_first_present(cols, REEF_ID_COLNAME_CANDIDATES)
if(is.na(REEF_COL)) stop("No reef id column found. Tried: ", paste(REEF_ID_COLNAME_CANDIDATES, collapse=", "))
dat0$reef_group <- factor(dat0[[REEF_COL]])
if(!("reef_group" %in% names(dat0))) stop("reef_group was not created. Check Reef column detection.")
cat("Using Reef ID column:", REEF_COL, " | n_groups=", length(levels(dat0$reef_group)), "\n")
OUTCOME <- "BleachResp4"

# --------------------------
# MAP PREDICTORS

# --------------------------
DEPTH <- pick_first_present(cols, c("depths","Depth","RAW_depths"))
ROTC  <- pick_first_present(cols, c("ROTC_.SS.","ROTCSS","RAW_ROTC_.SS.","ROTCSS"))
AC1   <- pick_first_present(cols, c("Acute1...10","Acute1...59","RAW_Acute1"))
DHW30 <- pick_first_present(cols, c("DHW_.l30.","DHW30","RAW_DHW_.l30."))
DTR30 <- pick_first_present(cols, c("DTR_.30.","DTR30","RAW_DTR_.30."))
TT    <- pick_first_present(cols, c("TT...12","TT...84","RAW_TT"))
WIND <- pick_first_present(cols, WIND_COL_CANDIDATES)
if(is.na(WIND)) stop("Wind column not found. Tried: ", paste(WIND_COL_CANDIDATES, collapse=", "))
base_terms <- c(DEPTH, ROTC, AC1, DHW30, DTR30, TT)
if(any(is.na(base_terms))) stop("Missing baseline predictors (DEPTH/ROTC/AC1/DHW30/DTR30/TT).")

# --------------------------
# MODEL SPECS (K4)

# --------------------------
model_specs <- list(
  Safaie_Base = base_terms,
  Sens_BasePlusWind = c(base_terms, WIND)  # includes DTR and temperature variables
)

# --------------------------
# FIT FE + ME

# --------------------------
model_rows <- list()
coef_rows  <- list()
dp_rows    <- list()
for(mname in names(model_specs)){
  terms <- model_specs[[mname]]
  d <- prep_model_data_K4(dat0, OUTCOME, "reef_group", preds=terms)
  if(is.null(d)) next
  f_full <- make_formula(OUTCOME, terms)
  f_null <- make_formula(OUTCOME, c())
  fe_full <- fit_clm_safe(f_full, d)
  fe_null <- fit_clm_safe(f_null, d)
  me_full <- fit_clmm_safe(f_full, d, "reef_group")
  me_null <- fit_clmm_safe(f_null, d, "reef_group")
  model_rows[[mname]] <- bind_rows(
    tibble(model=mname, model_type="FE", n=nrow(d), k=length(terms),
           AICc=AICc(fe_full),
           R2_mcfadden=mcfadden_r2(fe_full, fe_null),
           R2_nagelkerke=nagelkerke_r2(fe_full, fe_null),
           note=ifelse(is.null(fe_full),"FE fit failed","")),
    tibble(model=mname, model_type="ME", n=nrow(d), k=length(terms),
           AICc=AICc(me_full),
           R2_mcfadden=mcfadden_r2(me_full, me_null),
           R2_nagelkerke=nagelkerke_r2(me_full, me_null),
           note=ifelse(is.null(me_full),"ME fit failed",""))
  )
  fe_cf <- extract_slopes(fe_full, "FE") |> mutate(model=mname)
  me_cf <- extract_slopes(me_full, "ME") |> mutate(model=mname)
  cf <- bind_rows(fe_cf, me_cf) |>
    filter(term %in% terms) |>
    filter(!is.na(se) & is.finite(se)) |>
    add_wald_ci(0.95) |>
    mutate(effect_ci = sprintf("%.4f [%.4f, %.4f]%s", estimate, ci_low, ci_high,
                               ifelse(sig=="","", paste0(" ",sig))))
  coef_rows[[mname]] <- cf
  if(!is.null(fe_full)){
    for(trm in terms){
      dp_rows[[paste(mname,"FE",trm,sep="|")]] <-
        deltaP_highest(fe_full, d, trm, terms_all=terms) |>
        mutate(model=mname, model_type="FE")
    }
  }
  if(!is.null(me_full)){
    for(trm in terms){
      dp_rows[[paste(mname,"ME",trm,sep="|")]] <-
        deltaP_highest(me_full, d, trm, terms_all=terms, group_var="reef_group") |>
        mutate(model=mname, model_type="ME")
    }
  }
}
model_tbl <- bind_rows(model_rows) |>
  group_by(model_type) |>
  mutate(dAICc = if(all(is.na(AICc))) NA_real_ else AICc - min(AICc, na.rm=TRUE)) |>
  ungroup() |>
  arrange(model_type, dAICc)
coef_tbl <- bind_rows(coef_rows) |>
  left_join(model_tbl |> dplyr::select(model, model_type, n, k, AICc, dAICc, R2_mcfadden, R2_nagelkerke),
            by=c("model","model_type")) |>
  arrange(model_type, dAICc, model, term)
dp_tbl <- bind_rows(dp_rows) |>
  left_join(model_tbl |> dplyr::select(model, model_type, n, k, AICc, dAICc, R2_mcfadden, R2_nagelkerke),
            by=c("model","model_type")) |>
  arrange(model_type, dAICc, model, term)

# ============================================================
# RESIDUAL WIND TEST (FE + ME) for Sens_BasePlusWind

# ============================================================
terms_best <- model_specs$Sens_BasePlusWind
d_best <- prep_model_data_K4(dat0, OUTCOME, "reef_group", preds=terms_best)
controls_for_wind <- setdiff(terms_best, WIND)
# FE residuals
f_wind <- as.formula(paste0(WIND, " ~ ", paste(controls_for_wind, collapse=" + ")))
wind_lm <- lm(f_wind, data=d_best)
d_best$Wind_resid_FE <- resid(wind_lm)
# ME residuals (ensure reef_group is present in d_best; it is)
wind_lmer <- lme4::lmer(update(f_wind, . ~ . + (1|reef_group)), data=d_best, REML=TRUE)
d_best$Wind_resid_ME <- resid(wind_lmer)
# Ordinal with residual wind
f_ord_resid_fe <- as.formula(paste0(OUTCOME, " ~ ", paste(controls_for_wind, collapse=" + "), " + Wind_resid_FE"))
f_ord_resid_me <- as.formula(paste0(OUTCOME, " ~ ", paste(controls_for_wind, collapse=" + "), " + Wind_resid_ME"))
ord_resid_fe <- fit_clm_safe(f_ord_resid_fe, d_best)
ord_resid_me <- fit_clmm_safe(f_ord_resid_me, d_best, "reef_group")
resid_tbl <- bind_rows(
  extract_one_slope(ord_resid_fe, "Wind_resid_FE", "FE"),
  extract_one_slope(ord_resid_me, "Wind_resid_ME", "ME")
) |>
  mutate(
    test="Residual wind in OLR (controls include DHW30, DTR30, TT, Acute, depth, ROTC)",
    note=case_when(
      is.na(p) ~ "fit failed or term missing",
      p < 0.05 ~ "Residual wind significant (wind has signal beyond temperature-related predictors)",
      TRUE ~ "Residual wind not significant (wind mostly shares variance with temperature predictors)"
    )
  )
mediator_tbl <- bind_rows(
  tibble(model="Wind_lm_FE", metric=c("n","R2","Adj_R2","sigma"),
         value=c(nobs(wind_lm), summary(wind_lm)$r.squared, summary(wind_lm)$adj.r.squared, sigma(wind_lm))),
  tibble(model="Wind_lmer_ME", metric=c("n","RE_sd(reef)","sigma"),
         value=c(nobs(wind_lmer),
                 sqrt(as.numeric(lme4::VarCorr(wind_lmer)$reef_group[1,1])),
                 sigma(wind_lmer)))
)

# ============================================================
# APPROXIMATE MEDIATION VIA DHW30 (difference-in-coef diagnostic)

# ============================================================
controls_no_DHW <- c(DEPTH, ROTC, AC1, DTR30, TT)
need_med <- unique(c(OUTCOME,"reef_group",WIND,DHW30,controls_no_DHW))
d_med <- dat0[, need_med, drop=FALSE]
d_med <- d_med[complete.cases(d_med), , drop=FALSE]
# Total effect (exclude DHW30)
f_total <- as.formula(paste0(OUTCOME, " ~ ", paste(c(controls_no_DHW, WIND), collapse=" + ")))
f_null  <- as.formula(paste0(OUTCOME, " ~ 1"))
fe_total <- fit_clm_safe(f_total, d_med)
fe_total_null <- fit_clm_safe(f_null, d_med)
me_total <- fit_clmm_safe(f_total, d_med, "reef_group")
me_total_null <- fit_clmm_safe(f_null, d_med, "reef_group")
# Direct effect (include DHW30)
f_direct <- as.formula(paste0(OUTCOME, " ~ ", paste(c(controls_no_DHW, DHW30, WIND), collapse=" + ")))
fe_direct <- fit_clm_safe(f_direct, d_med)
fe_direct_null <- fit_clm_safe(f_null, d_med)
me_direct <- fit_clmm_safe(f_direct, d_med, "reef_group")
me_direct_null <- fit_clmm_safe(f_null, d_med, "reef_group")
extract_term_coef <- function(mod, term){
  if(is.null(mod)) return(c(estimate=NA_real_, p=NA_real_))
  sm <- summary(mod)$coefficients
  n_thres <- length(mod$alpha) %||% 0
  if(n_thres==0) n_thres <- mod$nTheta %||% 0
  slopes <- sm[(n_thres+1):nrow(sm), , drop=FALSE]
  if(!(term %in% rownames(slopes))) return(c(estimate=NA_real_, p=NA_real_))
  c(estimate=as.numeric(slopes[term,"Estimate"]),
    p=as.numeric(slopes[term,"Pr(>|z|)"]))
}
b_fe_total  <- extract_term_coef(fe_total, WIND)
b_fe_direct <- extract_term_coef(fe_direct, WIND)
b_me_total  <- extract_term_coef(me_total, WIND)
b_me_direct <- extract_term_coef(me_direct, WIND)
med_tbl <- tibble(
  model_type = c("FE","ME"),
  n = nrow(d_med),
  wind = WIND,
  beta_total  = c(b_fe_total["estimate"],  b_me_total["estimate"]),
  p_total     = c(b_fe_total["p"],         b_me_total["p"]),
  beta_direct = c(b_fe_direct["estimate"], b_me_direct["estimate"]),
  p_direct    = c(b_fe_direct["p"],        b_me_direct["p"])
) |>
  mutate(
    beta_indirect_approx = beta_total - beta_direct,
    prop_mediated_approx = ifelse(is.finite(beta_total) & beta_total!=0,
                                  (beta_total - beta_direct)/beta_total, NA_real_),
    R2_total_mcfadden = c(mcfadden_r2(fe_total, fe_total_null), mcfadden_r2(me_total, me_total_null)),
    R2_direct_mcfadden = c(mcfadden_r2(fe_direct, fe_direct_null), mcfadden_r2(me_direct, me_direct_null)),
    note = "Approx mediation via DHW30: (beta_total - beta_direct). Diagnostic only."
  )
# Mediator model: DHW30 ~ Wind + controls (FE/ME)
f_mediator <- as.formula(paste0(DHW30, " ~ ", paste(c(controls_no_DHW, WIND), collapse=" + ")))
med_lm <- lm(f_mediator, data=d_med)
med_lmer <- lme4::lmer(update(f_mediator, . ~ . + (1|reef_group)), data=d_med, REML=TRUE)
mediator2_tbl <- bind_rows(
  tibble(model="Mediator_lm_FE", term=WIND,
         estimate=coef(summary(med_lm))[WIND,"Estimate"],
         se=coef(summary(med_lm))[WIND,"Std. Error"],
         t=coef(summary(med_lm))[WIND,"t value"],
         p=coef(summary(med_lm))[WIND,"Pr(>|t|)"]),
  {
    fx <- lme4::fixef(med_lmer)
    se <- sqrt(diag(as.matrix(vcov(med_lmer))))
    z  <- fx[WIND]/se[WIND]
    tibble(model="Mediator_lmer_ME", term=WIND,
           estimate=fx[WIND], se=se[WIND], t=z, p=2*pnorm(-abs(z)))
  }
)

# --------------------------
# EXPORT TO EXCEL

# --------------------------
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "K4_Model_Table")
openxlsx::writeData(wb, "K4_Model_Table", model_tbl)
openxlsx::addWorksheet(wb, "K4_Coefficients")
openxlsx::writeData(wb, "K4_Coefficients", coef_tbl)
openxlsx::addWorksheet(wb, "K4_DeltaP_Highest")
openxlsx::writeData(wb, "K4_DeltaP_Highest", dp_tbl)
openxlsx::addWorksheet(wb, "Residual_Wind_Test")
openxlsx::writeData(wb, "Residual_Wind_Test", resid_tbl)
openxlsx::addWorksheet(wb, "Wind_Residualization_Models")
openxlsx::writeData(wb, "Wind_Residualization_Models", mediator_tbl)
openxlsx::addWorksheet(wb, "Mediation_Approx")
openxlsx::writeData(wb, "Mediation_Approx", med_tbl)
openxlsx::addWorksheet(wb, "Mediator_DHW30_Models")
openxlsx::writeData(wb, "Mediator_DHW30_Models", mediator2_tbl)
openxlsx::addWorksheet(wb, "Model_Specs")
openxlsx::writeData(wb, "Model_Specs",
                    tibble(model=names(model_specs),
                           k=vapply(model_specs, length, integer(1)),
                           terms=vapply(model_specs, function(v) paste(v, collapse=" + "), character(1))))
for(s in names(wb)){
  try(openxlsx::setColWidths(wb, s, cols=1:40, widths="auto"), silent=TRUE)
}
tryCatch({
  openxlsx::saveWorkbook(wb, OUT_XLSX, overwrite=TRUE)
  cat("\nDONE.\nExcel written to:\n", OUT_XLSX, "\n")
}, error=function(e){
  message("saveWorkbook failed: ", e$message)
  fallback <- file.path(OUT_DIR, paste0("K4_FE_ME_WindResidual_Mediation_FINAL_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".xlsx"))
  message("Trying fallback:\n", fallback)
  openxlsx::saveWorkbook(wb, fallback, overwrite=TRUE)
  cat("\nDONE.\nExcel written to fallback:\n", fallback, "\n")
})
