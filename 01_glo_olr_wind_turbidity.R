# =============================================================================
# GLO ordinal logistic regression (FE & ME); wind x turbidity interaction
# Part of: Wind effects on coral bleaching severity (Lapenis)
# Language: R
# Inputs : data/input/Supplementary_Data_S1_S2.xlsx (sheet "Table S2 - GLO", header on row 2)
# Outputs: data/output/OLR_FE_ME_windXturbidity/
# Depends: dplyr, forcats, ggplot2, janitor, openxlsx, ordinal, purrr, readxl, stringr, tibble, tidyr
# Notes  : Paths are set in the CONFIG/USER-SETTINGS block below; place input
#          files in data/input/ and run from the repository root.
# =============================================================================

############################################
# OLR SCRIPT: m = 6 ONLY
# wind × turbidity interaction only
#
# Features:
# - full environmental control set always included
# - one wind + one TCP + one TSA at a time
# - only wind:turbidity interaction
# - FE and ME models
# - circular longitude (lon_sin, lon_cos)

############################################

# -------------------------
# User settings
# -------------------------
# GLO analysis table = sheet "Table S2 - GLO" of the deposited
# Supplementary_Data_S1_S2.xlsx (the same file read by 11_glo_selection_robustness.R).
input_file <- "data/input/Supplementary_Data_S1_S2.xlsx"
sheet_name <- "Table S2 - GLO"
glo_skip   <- 1                 # header sits on the 2nd row of that sheet
out_dir    <- "data/output/OLR_FE_ME_windXturbidity"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_xlsx   <- file.path(out_dir, "OLR_FE_ME_m6_windXTurbidity_ONE_WORKBOOK.xlsx")
delta_cutoff <- 2
max_models_to_plot <- 3
# Mixed-effects settings
run_me_models <- TRUE
drop_sparse_realms_for_ME <- TRUE
min_obs_per_realm_ME <- 5
random_group <- "realm_name"

# -------------------------
# Packages
# -------------------------
pkgs <- c(
  "readxl","dplyr","tidyr","purrr","stringr","ordinal","ggplot2",
  "forcats","openxlsx","tibble","janitor"
)
# Deposited code does not auto-install packages. Check that dependencies are
# present and stop with guidance if not (install manually, or use renv::restore()
# / r_requirements.txt).
missing_pkgs <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(missing_pkgs) > 0) {
  stop("Missing required package(s): ", paste(missing_pkgs, collapse = ", "),
       "\n  Install with: install.packages(c(",
       paste(sprintf('"%s"', missing_pkgs), collapse = ", "), "))",
       "\n  or restore the recorded environment with renv::restore().",
       call. = FALSE)
}
library(readxl)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(ordinal)
library(ggplot2)
library(forcats)
library(openxlsx)
library(tibble)
library(janitor)

# -------------------------
# Helpers
# -------------------------
clean_sheet_base <- function(x) {
  x <- gsub("[:\\\\/\\?\\*\\[\\]]", "_", x)
  x <- gsub("\\s+", "_", x)
  x <- gsub("__+", "_", x)
  x
}
make_unique_sheet_name <- function(wb, base, suffix = "", max_len = 31) {
  existing <- tolower(openxlsx::sheets(wb))
  base <- clean_sheet_base(base)
  suffix <- clean_sheet_base(suffix)
  candidate <- paste0(base, suffix)
  candidate <- substr(candidate, 1, max_len)
  if (!(tolower(candidate) %in% existing)) return(candidate)
  i <- 1L
  repeat {
    tag <- sprintf("_%02d", i)
    cut_len <- max_len - nchar(tag)
    cand2 <- substr(candidate, 1, cut_len)
    cand2 <- paste0(cand2, tag)
    if (!(tolower(cand2) %in% existing)) return(cand2)
    i <- i + 1L
    if (i > 999L) stop("Could not create unique sheet name after 999 attempts.")
  }
}
formula_to_string <- function(f) paste(deparse(f, width.cutoff = 500L), collapse = " ")
safe_scale <- function(x) {
  if (all(is.na(x))) return(x)
  s <- stats::sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(x)
  as.numeric(scale(x))
}
aicc <- function(aic, k, n) {
  denom <- (n - k - 1)
  ifelse(denom > 0, aic + (2 * k * (k + 1)) / denom, NA_real_)
}
akaike_weights <- function(delta_metric) {
  w <- exp(-0.5 * delta_metric)
  w / sum(w, na.rm = TRUE)
}
mcfadden_r2 <- function(ll_model, ll_null) {
  if (is.na(ll_model) || is.na(ll_null) || ll_null == 0) return(NA_real_)
  1 - (as.numeric(ll_model) / as.numeric(ll_null))
}
extract_vars_from_formula_string <- function(formula_str) {
  f <- stats::as.formula(formula_str)
  tl <- attr(stats::terms(f), "term.labels")
  paste(tl, collapse = "; ")
}
fit_clm_diag <- function(formula, data) {
  tryCatch(
    {
      mod <- suppressWarnings(ordinal::clm(formula, data = data, link = "logit", Hess = TRUE))
      list(ok = TRUE, model = mod, error = NA_character_)
    },
    error = function(e) {
      list(ok = FALSE, model = NULL, error = conditionMessage(e))
    }
  )
}
fit_clmm_diag <- function(formula, data) {
  tryCatch(
    {
      mod <- suppressWarnings(
        ordinal::clmm(
          formula,
          data = data,
          link = "logit",
          Hess = TRUE,
          nAGQ = 1
        )
      )
      list(ok = TRUE, model = mod, error = NA_character_)
    },
    error = function(e) {
      message("clmm failed for formula: ", paste(deparse(formula), collapse = " "))
      message("Reason: ", conditionMessage(e))
      list(ok = FALSE, model = NULL, error = conditionMessage(e))
    }
  )
}
coef_table_ordinal <- function(model) {
  beta <- stats::coef(model)
  nm <- names(beta)
  is_threshold <- grepl("\\|", nm)
  vc <- tryCatch(as.matrix(stats::vcov(model)), error = function(e) NULL)
  if (is.null(vc)) return(NULL)
  se <- sqrt(diag(vc))
  tibble::tibble(
    term = nm[!is_threshold],
    estimate = as.numeric(beta[!is_threshold]),
    std_error = as.numeric(se[!is_threshold])
  ) %>%
    dplyr::mutate(
      OR = exp(estimate),
      lo = exp(estimate - 1.96 * std_error),
      hi = exp(estimate + 1.96 * std_error),
      z = estimate / std_error,
      p = 2 * (1 - stats::pnorm(abs(z)))
    ) %>%
    dplyr::arrange(dplyr::desc(abs(estimate)))
}

# -------------------------
# FE-only prediction helpers
# -------------------------
compute_eta_manual <- function(model, newdata) {
  tt <- tryCatch(stats::terms(model), error = function(e) NULL)
  if (is.null(tt)) stop("Could not extract terms(model).")
  contr <- NULL
  if (!is.null(model$contrasts)) contr <- model$contrasts
  X <- stats::model.matrix(tt, newdata, contrasts.arg = contr)
  cf <- stats::coef(model)
  slope <- cf[!grepl("\\|", names(cf))]
  slope <- slope[names(slope) != "(Intercept)"]
  common <- intersect(colnames(X), names(slope))
  if (length(common) == 0) stop("No overlap between model matrix columns and slope coefficient names.")
  as.numeric(X[, common, drop = FALSE] %*% slope[common])
}
predict_prob_manual_clm <- function(model, newdata) {
  eta <- compute_eta_manual(model, newdata)
  cf <- stats::coef(model)
  thr <- cf[grepl("\\|", names(cf))]
  if (length(thr) < 1) stop("Could not find thresholds.")
  thr <- as.numeric(thr)
  cum <- sapply(thr, function(tk) stats::plogis(tk - eta))
  cum <- as.matrix(cum)
  n <- length(eta)
  Jm1 <- ncol(cum)
  J <- Jm1 + 1
  probs <- matrix(NA_real_, nrow = n, ncol = J)
  probs[, 1] <- cum[, 1]
  if (Jm1 >= 2) for (k in 2:Jm1) probs[, k] <- cum[, k] - cum[, k - 1]
  probs[, J] <- 1 - cum[, Jm1]
  levs <- levels(model$y)
  if (!is.null(levs) && length(levs) == J) colnames(probs) <- levs
  probs
}
full_impact_top_model <- function(model, data) {
  terms_used <- attr(terms(formula(model)), "term.labels")
  main_terms <- terms_used[!grepl(":", terms_used)]
  main_terms <- main_terms[main_terms %in% names(data)]
  if (length(main_terms) == 0) return(tibble::tibble())
  num_terms <- main_terms[sapply(data[main_terms], is.numeric)]
  if (length(num_terms) == 0) return(tibble::tibble())
  get_p_high <- function(probs) probs[, ncol(probs)]
  out <- list()
  for (v in num_terms) {
    d_lo <- data
    d_hi <- data
    d_lo[[v]] <- -1
    d_hi[[v]] <- 1
    p_lo <- get_p_high(predict_prob_manual_clm(model, d_lo))
    p_hi <- get_p_high(predict_prob_manual_clm(model, d_hi))
    out[[length(out) + 1]] <- tibble::tibble(
      variable = v,
      contrast = "-1SD to +1SD",
      mean_delta_P_highest = mean(p_hi - p_lo, na.rm = TRUE),
      median_delta_P_highest = stats::median(p_hi - p_lo, na.rm = TRUE)
    )
  }
  dplyr::bind_rows(out) %>%
    dplyr::arrange(dplyr::desc(abs(mean_delta_P_highest)))
}
save_forest_plot <- function(model, title, out_path) {
  beta <- stats::coef(model)
  nm <- names(beta)
  is_threshold <- grepl("\\|", nm)
  betas <- beta[!is_threshold]
  if (length(betas) == 0) return(FALSE)
  vc <- tryCatch(stats::vcov(model), error = function(e) NULL)
  if (is.null(vc)) return(FALSE)
  se <- sqrt(diag(vc))[!is_threshold]
  dfp <- tibble::tibble(
    term = names(betas),
    estimate = as.numeric(betas),
    std.error = as.numeric(se)
  ) %>%
    dplyr::mutate(
      OR = exp(estimate),
      lo = exp(estimate - 1.96 * std.error),
      hi = exp(estimate + 1.96 * std.error)
    ) %>%
    dplyr::mutate(term = forcats::fct_reorder(term, abs(estimate))) %>%
    dplyr::arrange(dplyr::desc(abs(estimate)))
  p <- ggplot2::ggplot(dfp, ggplot2::aes(x = OR, y = term)) +
    ggplot2::geom_vline(xintercept = 1, linetype = 2) +
    ggplot2::geom_point() +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = lo, xmax = hi), height = 0.2) +
    ggplot2::scale_x_log10() +
    ggplot2::labs(title = title, x = "Odds Ratio (log scale) with 95% CI", y = NULL) +
    ggplot2::theme_bw()
  ggplot2::ggsave(out_path, p, width = 10, height = 7, dpi = 300)
  TRUE
}

# -------------------------
# Formula builder: m = 6 only, wind:turbidity only
# -------------------------
build_formulas_m6_windXTurbidity <- function(response, wind_set, tc_set, tsa_set, other_pool) {
  stopifnot(length(other_pool) == 6)
  stopifnot("turbidity" %in% other_pool)
  sst_block <- c("sst_pc1", "sst_pc2", "sst_pc3")
  oth <- other_pool
  formulas <- list()
  meta <- list()
  idx <- 1L
  for (w in wind_set) {
    for (tc in tc_set) {
      for (tsa in tsa_set) {
        preds <- c(sst_block, w, tc, tsa, oth)
        preds <- unique(preds)
        interaction_terms <- c(paste0(w, ":turbidity"))
        rhs_terms <- c(preds, interaction_terms)
        fml <- stats::as.formula(paste(response, "~", paste(rhs_terms, collapse = " + ")))
        formulas[[idx]] <- fml
        meta[[idx]] <- tibble::tibble(
          model_id = idx,
          selected_wind = w,
          selected_tcp = tc,
          selected_tsa = tsa,
          m_other = 6,
          interaction_mode = "wind_x_turbidity",
          includes_turbidity = TRUE,
          selected_optional = paste(oth, collapse = "; ")
        )
        idx <- idx + 1L
      }
    }
  }
  list(
    formulas = formulas,
    meta = dplyr::bind_rows(meta)
  )
}
add_random_intercept <- function(formula_obj, group_var = "realm_name") {
  rhs <- paste(deparse(formula_obj[[3]], width.cutoff = 500L), collapse = " ")
  stats::as.formula(paste(as.character(formula_obj[[2]]), "~", rhs, "+ (1|", group_var, ")"))
}

# -------------------------
# Load and prepare data
# -------------------------
raw <- readxl::read_excel(input_file, sheet = sheet_name, skip = glo_skip) %>% janitor::clean_names()
required_cols_base <- c(
  "bleaching_categorical",
  "realm_name",
  "latitude",
  "longitude",
  "wind_mean_6m","wind_mean_12m","wind_mean_1993_2020",
  "tcpower_6m_400km","tcpower_12m_400km","tcpower_1993_2020_400km",
  "tsa_dhw","tsa_dhw_standard_deviation","tsa_dhwmax","tsa_dhwmean",
  "sst_pc1","sst_pc2","sst_pc3",
  "distance_to_shore","exposure","turbidity"
)
missing_base <- setdiff(required_cols_base, names(raw))
if (length(missing_base) > 0) stop(paste0("Missing columns:\n- ", paste(missing_base, collapse = "\n- ")))
df <- raw %>% dplyr::select(dplyr::all_of(required_cols_base))
# circular longitude representation
df <- df %>%
  dplyr::mutate(
    lon_rad = longitude * pi / 180,
    lon_sin = sin(lon_rad),
    lon_cos = cos(lon_rad)
  )
# Response
if (is.numeric(df$bleaching_categorical)) {
  df <- df %>%
    dplyr::mutate(
      bleaching_categorical = dplyr::case_when(
        bleaching_categorical == 5  ~ "Low",
        bleaching_categorical == 30 ~ "Moderate",
        bleaching_categorical == 75 ~ "Severe",
        TRUE ~ as.character(bleaching_categorical)
      ),
      bleaching_categorical = factor(
        bleaching_categorical,
        levels = c("Low","Moderate","Severe"),
        ordered = TRUE
      )
    )
} else {
  df <- df %>%
    dplyr::mutate(
      bleaching_categorical = as.character(bleaching_categorical),
      bleaching_categorical = dplyr::case_when(
        bleaching_categorical %in% c("5","Low","low") ~ "Low",
        bleaching_categorical %in% c("30","Moderate","moderate","Medium","medium") ~ "Moderate",
        bleaching_categorical %in% c("75","Severe","severe","High","high") ~ "Severe",
        TRUE ~ bleaching_categorical
      ),
      bleaching_categorical = factor(
        bleaching_categorical,
        levels = c("Low","Moderate","Severe"),
        ordered = TRUE
      )
    )
}
df <- df %>% dplyr::filter(!is.na(bleaching_categorical))
# Exposure and grouping
df <- df %>%
  dplyr::mutate(
    exposure = as.character(exposure),
    exposure = dplyr::case_when(
      exposure %in% c("0","protected","Protected") ~ "0",
      exposure %in% c("1","somewhat protected","Somewhat protected","Somewhat_Protected") ~ "1",
      exposure %in% c("3","exposed","Exposed") ~ "3",
      TRUE ~ exposure
    ),
    exposure = as.numeric(exposure),
    realm_name = factor(realm_name)
  )
# Scale numeric predictors
is_num <- sapply(df, is.numeric)
num_cols <- names(df)[is_num]
df_scaled <- df
df_scaled[num_cols] <- lapply(df_scaled[num_cols], safe_scale)
analysis_cols <- c(
  "bleaching_categorical", "realm_name",
  "latitude", "lon_sin", "lon_cos",
  "wind_mean_6m","wind_mean_12m","wind_mean_1993_2020",
  "tcpower_6m_400km","tcpower_12m_400km","tcpower_1993_2020_400km",
  "tsa_dhw","tsa_dhw_standard_deviation","tsa_dhwmax","tsa_dhwmean",
  "sst_pc1","sst_pc2","sst_pc3",
  "distance_to_shore","exposure","turbidity"
)
df_scaled_base <- df_scaled %>%
  dplyr::select(dplyr::all_of(analysis_cols)) %>%
  dplyr::filter(stats::complete.cases(.))
cat("\nRows after complete-case filtering (FE): ", nrow(df_scaled_base), "\n", sep = "")
# Prepare ME dataset
realm_tab_fe <- as.data.frame(table(df_scaled_base[[random_group]]))
colnames(realm_tab_fe) <- c("realm_name", "n")
if (run_me_models) {
  if (drop_sparse_realms_for_ME) {
    keep_realms <- realm_tab_fe$realm_name[realm_tab_fe$n >= min_obs_per_realm_ME]
    df_scaled_me <- df_scaled_base %>%
      dplyr::filter(.data[[random_group]] %in% keep_realms) %>%
      droplevels()
  } else {
    df_scaled_me <- df_scaled_base
  }
  realm_tab_me <- as.data.frame(table(df_scaled_me[[random_group]]))
  colnames(realm_tab_me) <- c("realm_name", "n")
  cat("Rows after ME realm filtering: ", nrow(df_scaled_me), "\n", sep = "")
  cat("Realms retained for ME: ", nrow(realm_tab_me), "\n", sep = "")
} else {
  df_scaled_me <- NULL
  realm_tab_me <- NULL
}

# -------------------------
# Simple ME diagnostic model
# -------------------------
simple_me_result <- NULL
if (run_me_models && !is.null(df_scaled_me)) {
  test_formula <- stats::as.formula(paste0("bleaching_categorical ~ sst_pc1 + (1|", random_group, ")"))
  simple_me_result <- fit_clmm_diag(test_formula, df_scaled_me)
}

# -------------------------
# Predictor sets
# -------------------------
wind_set <- c("wind_mean_6m","wind_mean_12m","wind_mean_1993_2020")
tc_set   <- c("tcpower_6m_400km","tcpower_12m_400km","tcpower_1993_2020_400km")
tsa_set  <- c("tsa_dhw","tsa_dhw_standard_deviation","tsa_dhwmax","tsa_dhwmean")
# All 6 controls always included
other_pool <- c("latitude","lon_sin","lon_cos","distance_to_shore","exposure","turbidity")

# -------------------------
# Universe evaluator
# -------------------------
evaluate_model_universe <- function(formulas, formula_meta, data, scenario_label,
                                    model_type = c("FE", "ME"),
                                    delta_cutoff = 2,
                                    max_plots = 3,
                                    plot_dir = out_dir,
                                    random_group = "realm_name") {
  model_type <- match.arg(model_type)
  cat("\n====================================================\n")
  cat("Running:", scenario_label, "|", model_type, "\n")
  cat("Total candidate models:", length(formulas), "\n")
  cat("Rows in dataset:", nrow(data), "\n")
  cat("====================================================\n")
  if (model_type == "FE") {
    null_fit <- fit_clm_diag(bleaching_categorical ~ 1, data)
  } else {
    null_fit <- fit_clmm_diag(stats::as.formula(
      paste0("bleaching_categorical ~ 1 + (1|", random_group, ")")
    ), data)
  }
  ll_null <- if (isTRUE(null_fit$ok)) as.numeric(stats::logLik(null_fit$model)) else NA_real_
  pb <- utils::txtProgressBar(min = 0, max = length(formulas), style = 3)
  on.exit(close(pb), add = TRUE)
  success_rows <- list()
  fail_rows <- list()
  coef_tables <- list()
  impact_tables <- list()
  plot_paths <- character(0)
  for (i in seq_along(formulas)) {
    f <- formulas[[i]]
    fit_formula <- if (model_type == "FE") f else add_random_intercept(f, random_group)
    fit_obj <- if (model_type == "FE") fit_clm_diag(fit_formula, data) else fit_clmm_diag(fit_formula, data)
    if (isTRUE(fit_obj$ok)) {
      mod <- fit_obj$model
      ll <- tryCatch(as.numeric(stats::logLik(mod)), error = function(e) NA_real_)
      aic_val <- tryCatch(stats::AIC(mod), error = function(e) NA_real_)
      k <- tryCatch(length(stats::coef(mod)), error = function(e) NA_integer_)
      n <- nrow(data)
      row_tbl <- tibble::tibble(
        scenario = scenario_label,
        model_type = model_type,
        model_id = i,
        formula = formula_to_string(fit_formula),
        n = n,
        k = as.integer(k),
        logLik = ll,
        AIC = as.numeric(aic_val),
        AICc = if (model_type == "FE") as.numeric(aicc(as.numeric(aic_val), as.integer(k), n)) else NA_real_,
        McFadden_R2 = as.numeric(mcfadden_r2(ll, ll_null))
      ) %>%
        dplyr::left_join(formula_meta, by = "model_id")
      success_rows[[length(success_rows) + 1]] <- row_tbl
      ct <- coef_table_ordinal(mod)
      if (!is.null(ct)) coef_tables[[paste0("model_", i)]] <- ct
    } else {
      fail_rows[[length(fail_rows) + 1]] <- tibble::tibble(
        scenario = scenario_label,
        model_type = model_type,
        model_id = i,
        formula = formula_to_string(fit_formula),
        error_message = fit_obj$error
      ) %>%
        dplyr::left_join(formula_meta, by = "model_id")
    }
    if (i %% 25 == 0 || i == length(formulas)) utils::setTxtProgressBar(pb, i)
  }
  ranking <- dplyr::bind_rows(success_rows)
  if (nrow(ranking) == 0) {
    message("No successful ", model_type, " fits in scenario: ", scenario_label)
    return(list(
      ranking = tibble::tibble(),
      best = tibble::tibble(),
      var_summary = tibble::tibble(),
      coef_tables = coef_tables,
      impact_tables = impact_tables,
      plot_paths = plot_paths,
      fail_log = dplyr::bind_rows(fail_rows),
      summary_tbl = tibble::tibble(
        scenario = scenario_label,
        model_type = model_type,
        total_candidate_models = length(formulas),
        successful_models = 0,
        failed_models = nrow(dplyr::bind_rows(fail_rows)),
        success_fraction = 0
      )
    ))
  }
  if (model_type == "FE") {
    ranking <- ranking %>%
      dplyr::arrange(AICc, AIC) %>%
      dplyr::mutate(
        delta_metric = AICc - min(AICc, na.rm = TRUE),
        akaike_weight = akaike_weights(delta_metric),
        selected_vars = purrr::map_chr(formula, extract_vars_from_formula_string)
      ) %>%
      dplyr::rename(delta_AICc = delta_metric)
    best <- ranking %>%
      dplyr::filter(delta_AICc < delta_cutoff)
  } else {
    ranking <- ranking %>%
      dplyr::arrange(AIC) %>%
      dplyr::mutate(
        delta_metric = AIC - min(AIC, na.rm = TRUE),
        akaike_weight = akaike_weights(delta_metric),
        selected_vars = purrr::map_chr(formula, extract_vars_from_formula_string)
      ) %>%
      dplyr::rename(delta_AIC = delta_metric)
    best <- ranking %>%
      dplyr::filter(delta_AIC < delta_cutoff)
  }
  if (nrow(best) == 0) {
    best <- ranking %>% dplyr::slice(1)
  }
  var_summary <- NULL
  if (nrow(best) > 0) {
    var_summary <- best %>%
      dplyr::select(model_id, selected_vars, akaike_weight) %>%
      tidyr::separate_rows(selected_vars, sep = ";\\s*") %>%
      dplyr::group_by(selected_vars) %>%
      dplyr::summarise(
        n_models = dplyr::n(),
        frac_models = n_models / nrow(best),
        sum_akaike_weight = sum(akaike_weight, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::rename(variable = selected_vars) %>%
      dplyr::arrange(dplyr::desc(sum_akaike_weight), dplyr::desc(n_models))
  }
  plot_tbl <- best %>% dplyr::slice_head(n = min(max_plots, nrow(best)))
  for (j in seq_len(nrow(plot_tbl))) {
    f_j <- stats::as.formula(plot_tbl$formula[j])
    fit_obj_j <- if (model_type == "FE") fit_clm_diag(f_j, data) else fit_clmm_diag(f_j, data)
    if (!isTRUE(fit_obj_j$ok)) next
    mod_j <- fit_obj_j$model
    title <- if (model_type == "FE") {
      paste0(
        scenario_label, " | FE | rank ", j,
        " | ΔAICc=", round(plot_tbl$delta_AICc[j], 2),
        " | AICc=", round(plot_tbl$AICc[j], 2),
        " | McFadden R2=", round(plot_tbl$McFadden_R2[j], 3),
        " | wind_x_turbidity"
      )
    } else {
      paste0(
        scenario_label, " | ME | rank ", j,
        " | ΔAIC=", round(plot_tbl$delta_AIC[j], 2),
        " | AIC=", round(plot_tbl$AIC[j], 2),
        " | McFadden R2=", round(plot_tbl$McFadden_R2[j], 3),
        " | wind_x_turbidity"
      )
    }
    if (model_type == "FE") {
      img_path <- file.path(
        plot_dir,
        paste0("forest_", clean_sheet_base(scenario_label), "_FE_", j, ".png")
      )
      ok <- save_forest_plot(mod_j, title, img_path)
      if (ok) plot_paths <- c(plot_paths, img_path)
      impact_tables[[paste0("rank", j)]] <- full_impact_top_model(mod_j, data)
    }
  }
  summary_tbl <- tibble::tibble(
    scenario = scenario_label,
    model_type = model_type,
    total_candidate_models = length(formulas),
    successful_models = nrow(ranking),
    failed_models = nrow(dplyr::bind_rows(fail_rows)),
    success_fraction = nrow(ranking) / length(formulas)
  )
  list(
    ranking = ranking,
    best = best,
    var_summary = var_summary,
    coef_tables = coef_tables,
    impact_tables = impact_tables,
    plot_paths = plot_paths,
    fail_log = dplyr::bind_rows(fail_rows),
    summary_tbl = summary_tbl
  )
}

# -------------------------
# Workbook writer
# -------------------------
write_universe_to_workbook <- function(wb, scenario_label, res) {
  base <- clean_sheet_base(scenario_label)
  sh_rank <- make_unique_sheet_name(wb, base, "_rank")
  addWorksheet(wb, sh_rank)
  writeData(wb, sh_rank, res$ranking)
  sh_best <- make_unique_sheet_name(wb, base, "_best")
  addWorksheet(wb, sh_best)
  writeData(wb, sh_best, res$best)
  if (!is.null(res$var_summary) && nrow(res$var_summary) > 0) {
    sh_var <- make_unique_sheet_name(wb, base, "_vars")
    addWorksheet(wb, sh_var)
    writeData(wb, sh_var, res$var_summary)
  }
  if (length(res$coef_tables) > 0) {
    sh_coef <- make_unique_sheet_name(wb, base, "_coef")
    addWorksheet(wb, sh_coef)
    row <- 1
    for (nm in names(res$coef_tables)) {
      writeData(wb, sh_coef, paste("Coefficients:", nm), startRow = row, startCol = 1)
      row <- row + 1
      writeData(wb, sh_coef, res$coef_tables[[nm]], startRow = row, startCol = 1)
      row <- row + nrow(res$coef_tables[[nm]]) + 3
    }
  }
  if (length(res$impact_tables) > 0) {
    sh_imp <- make_unique_sheet_name(wb, base, "_impact")
    addWorksheet(wb, sh_imp)
    row <- 1
    for (nm in names(res$impact_tables)) {
      writeData(wb, sh_imp, paste("Full impact:", nm), startRow = row, startCol = 1)
      row <- row + 1
      writeData(wb, sh_imp, res$impact_tables[[nm]], startRow = row, startCol = 1)
      row <- row + nrow(res$impact_tables[[nm]]) + 3
    }
  }
  if (length(res$plot_paths) > 0) {
    sh_plot <- make_unique_sheet_name(wb, base, "_plots")
    addWorksheet(wb, sh_plot)
    writeData(wb, sh_plot, paste0("Forest plots for: ", scenario_label), startRow = 1, startCol = 1)
    r <- 3
    for (p in res$plot_paths) {
      insertImage(wb, sh_plot, file = p, startRow = r, startCol = 1,
                  width = 9.5, height = 6.5, units = "in")
      r <- r + 35
    }
  }
  sh_fail <- make_unique_sheet_name(wb, base, "_fail")
  addWorksheet(wb, sh_fail)
  writeData(wb, sh_fail, res$fail_log)
}

# -------------------------
# Run analysis
# -------------------------
wb <- openxlsx::createWorkbook()
# diagnostics sheet
addWorksheet(wb, "diagnostics")
writeData(wb, "diagnostics", "Realm sizes in FE dataset", startRow = 1, startCol = 1)
writeData(wb, "diagnostics", realm_tab_fe, startRow = 2, startCol = 1)
if (!is.null(realm_tab_me)) {
  writeData(wb, "diagnostics", "Realm sizes in ME dataset", startRow = 1, startCol = 5)
  writeData(wb, "diagnostics", realm_tab_me, startRow = 2, startCol = 5)
}
diag_meta <- tibble::tibble(
  metric = c(
    "rows_FE",
    "rows_ME",
    "n_realms_FE",
    "n_realms_ME",
    "simple_ME_model_success",
    "simple_ME_model_error",
    "interaction_mode"
  ),
  value = c(
    nrow(df_scaled_base),
    ifelse(is.null(df_scaled_me), NA, nrow(df_scaled_me)),
    nrow(realm_tab_fe),
    ifelse(is.null(realm_tab_me), NA, nrow(realm_tab_me)),
    ifelse(is.null(simple_me_result), NA, simple_me_result$ok),
    ifelse(is.null(simple_me_result), NA, ifelse(isTRUE(simple_me_result$ok), "", simple_me_result$error)),
    "wind_x_turbidity"
  )
)
writeData(wb, "diagnostics", diag_meta, startRow = 20, startCol = 1)
# build formulas only once (m = 6 only)
built <- build_formulas_m6_windXTurbidity(
  response = "bleaching_categorical",
  wind_set = wind_set,
  tc_set = tc_set,
  tsa_set = tsa_set,
  other_pool = other_pool
)
formulas <- built$formulas
formula_meta <- built$meta
all_summaries <- list()
# FE
res_fe <- evaluate_model_universe(
  formulas = formulas,
  formula_meta = formula_meta,
  data = df_scaled_base,
  scenario_label = "m6_FE",
  model_type = "FE",
  delta_cutoff = delta_cutoff,
  max_plots = max_models_to_plot,
  plot_dir = out_dir,
  random_group = random_group
)
write_universe_to_workbook(wb, "m6_FE", res_fe)
all_summaries[[length(all_summaries) + 1]] <- res_fe$summary_tbl
# ME
if (run_me_models && !is.null(df_scaled_me)) {
  res_me <- evaluate_model_universe(
    formulas = formulas,
    formula_meta = formula_meta,
    data = df_scaled_me,
    scenario_label = "m6_ME",
    model_type = "ME",
    delta_cutoff = delta_cutoff,
    max_plots = max_models_to_plot,
    plot_dir = out_dir,
    random_group = random_group
  )
  write_universe_to_workbook(wb, "m6_ME", res_me)
  all_summaries[[length(all_summaries) + 1]] <- res_me$summary_tbl
}
# summary sheet
addWorksheet(wb, "summary")
writeData(wb, "summary", dplyr::bind_rows(all_summaries))
openxlsx::saveWorkbook(wb, out_xlsx, overwrite = TRUE)
cat("\nWorkbook saved to:\n", out_xlsx, "\n", sep = "")
