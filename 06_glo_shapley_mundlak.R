# =============================================================================
# GLO Shapley decomposition of McFadden R-squared (4 groups) + Mundlak within-region wind
# Part of: Wind effects on coral bleaching severity (Lapenis & Jiang)
# Language: R
# Inputs : data/input/glo_bleaching_variables_PCA.xlsx
# Outputs: data/output/glo_shapley_mundlak_results.xlsx
# Depends: readxl, writexl, ordinal
# Notes  : Paths are set in the CONFIG/USER-SETTINGS block below; place input
#          files in data/input/ and run from the repository root.
# =============================================================================

###############################################################################
# GLO_shapley_mundlak.R
#
# GLO ordinal-bleaching analysis for the wind paper. Produces:
#   (1) brief OLR (proportional-odds) reference model  [spec follows Safaie 2018]
#   (2) Shapley variance decomposition of McFadden R^2 across 4 predictor GROUPS
#       - POOLED, and
#       - WITHIN-region (deviance beyond region fixed effects), for realm & basin
#   (3) within-region causal checks on the wind coefficient:
#       - fixed-effect-within estimator (region as fixed dummies)
#       - Mundlak split (within = wind - region mean; between = region mean)
#       repeated for realm (primary) and ocean basin (robustness)
#
# Response: bleaching_categorical {5,30,75} -> ordered {Low<Moderate<Severe}
# Groups:   Thermal = tsa_dhw + sst_pc1/2/3 ; Wind = wind_mean_6m ;
#           Cyclone = tcpower_1993_2020_400km ; Other = distance_to_shore + exposure
# Anchor wind for within/Mundlak = wind_mean_6m.
#
# Deps: readxl, writexl, ordinal

###############################################################################
pkgs <- c("readxl","writexl","ordinal")
missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0)
  stop("Missing required package(s): ", paste(missing_pkgs, collapse = ", "),
       "\n  Install them, or restore the recorded environment with renv::restore().",
       call. = FALSE)
invisible(lapply(pkgs,require,character.only=TRUE))
## ---- CONFIG ----
# GLO analysis table = same data as sheet "Table S2 - GLO" of the published
# Supplementary_Data_S1_S2.xlsx (read by 11), exported as a single flat sheet.
GLO_FILE <- "data/input/glo_bleaching_variables_PCA.xlsx"
OUT_FILE <- "data/output/glo_shapley_mundlak_results.xlsx"
MIN_REGION <- 10          # drop regions smaller than this for WITHIN-region fits
GROUPS <- list(
  Thermal = c("tsa_dhw","sst_pc1","sst_pc2","sst_pc3"),
  Wind    = c("wind_mean_6m"),
  Cyclone = c("tcpower_1993_2020_400km"),
  Other   = c("distance_to_shore","exposure"))
ALLP   <- unique(unlist(GROUPS))
CTRL   <- setdiff(ALLP, "wind_mean_6m")     # controls for the wind-coefficient fits
## ---- load, basin, response, scale ----
g <- as.data.frame(readxl::read_excel(GLO_FILE))
basin_of <- function(r) ifelse(grepl("Atlantic",r),"Atlantic",
                        ifelse(r=="Western Indo-Pacific","Indian","Pacific"))
g$basin <- basin_of(g$realm_name)
g$y <- factor(c(`5`="Low",`30`="Moderate",`75`="Severe")[as.character(g$bleaching_categorical)],
              levels=c("Low","Moderate","Severe"), ordered=TRUE)
g <- g[stats::complete.cases(g[,c("y",ALLP,"realm_name","basin")]),]
for(c in ALLP) g[[c]] <- as.numeric(scale(g[[c]]))      # z-score (R^2 invariant; coefs per-SD)
g$realm_name <- factor(g$realm_name); g$basin <- factor(g$basin)
cat(sprintf("n=%d | basins: %s\n", nrow(g),
            paste(names(table(g$basin)),table(g$basin),collapse=", ")))
## ---- helpers ----
ll_of <- function(rhs,data){
  m <- tryCatch(ordinal::clm(stats::as.formula(paste("y ~",rhs)),data=data,link="logit"),
                error=function(e) NULL)
  if(is.null(m)) NA_real_ else as.numeric(stats::logLik(m))
}
shapley <- function(groups, base_terms, data, ll_null){
  keys <- names(groups); m <- length(keys)
  subs <- unlist(lapply(0:m,function(r) combn(keys,r,simplify=FALSE)),recursive=FALSE)
  K <- function(S) paste(sort(S),collapse="|")
  R <- setNames(numeric(length(subs)), sapply(subs,K))
  for(S in subs){
    rhs <- paste(c(base_terms, unlist(groups[S])), collapse="+"); if(rhs=="") rhs<-"1"
    R[K(S)] <- 1 - ll_of(rhs,data)/ll_null
  }
  phi <- setNames(numeric(m),keys)
  for(gn in keys){ oth <- setdiff(keys,gn)
    for(r in 0:(m-1)) for(S in combn(oth,r,simplify=FALSE)){
      w <- factorial(length(S))*factorial(m-length(S)-1)/factorial(m)
      phi[gn] <- phi[gn] + w*(R[K(c(S,gn))]-R[K(S)]) } }
  full <- R[K(keys)]; data.frame(group=keys, R2=phi, pct=100*phi/full, row.names=NULL,
                                 check.names=FALSE)
}
wind_fit <- function(rhs,data){
  m <- ordinal::clm(stats::as.formula(paste("y ~",rhs)),data=data,link="logit")
  s <- summary(m)$coefficients
  rn <- intersect(c("wind_mean_6m","wdev","wbar"), rownames(s))
  out <- s[rn,c("Estimate","Std. Error","Pr(>|z|)"),drop=FALSE]
  data.frame(term=rownames(out), beta=out[,1], OR=exp(out[,1]), p=out[,3], row.names=NULL)
}
## ---- (1) brief OLR reference (spec per Safaie 2018) ----
olr <- ordinal::clm(y ~ tsa_dhw+sst_pc1+sst_pc2+sst_pc3+wind_mean_6m+
                      tcpower_1993_2020_400km+distance_to_shore+exposure, data=g, link="logit")
cat("\n--- (1) OLR reference model ---\n"); print(summary(olr)$coefficients)
mcf <- 1 - as.numeric(logLik(olr))/ll_of("1",g)
cat(sprintf("McFadden R2 = %.4f\n", mcf))
## ---- (2) POOLED Shapley ----
cat("\n--- (2a) POOLED Shapley ---\n")
sh_pool <- shapley(GROUPS, character(0), g, ll_of("1",g)); print(sh_pool, digits=3)
## ---- (2b) WITHIN-region Shapley (beyond region FE) ----
gr <- g[g$realm_name %in% names(which(table(g$realm_name)>=MIN_REGION)),]
gr$realm_name <- droplevels(gr$realm_name)
sh_realm <- shapley(GROUPS, "realm_name", gr, ll_of("realm_name",gr))
sh_basin <- shapley(GROUPS, "basin",      g,  ll_of("basin",g))
cat("\n--- (2b) WITHIN-realm Shapley ---\n");  print(sh_realm, digits=3)
cat("\n--- (2b) WITHIN-basin Shapley ---\n");  print(sh_basin, digits=3)
## ---- (3) within-region wind coefficient: FE-within & Mundlak ----
fe_mund <- function(data, gvar){
  data <- data[data[[gvar]] %in% names(which(table(data[[gvar]])>=MIN_REGION)),]
  data[[gvar]] <- droplevels(factor(data[[gvar]]))
  fe  <- wind_fit(paste(c("wind_mean_6m",CTRL,gvar),collapse="+"), data)
  mu  <- ave(data$wind_mean_6m, data[[gvar]])           # region-mean wind
  data$wbar <- mu; data$wdev <- data$wind_mean_6m - mu
  mund <- wind_fit(paste(c("wdev","wbar",CTRL),collapse="+"), data)
  list(fe=fe, mund=mund)
}
rr <- fe_mund(g,"realm_name"); bb <- fe_mund(g,"basin")
cat("\n--- (3) FE-within & Mundlak: REALM ---\n"); print(rr$fe,digits=3); print(rr$mund,digits=3)
cat("\n--- (3) FE-within & Mundlak: BASIN ---\n"); print(bb$fe,digits=3); print(bb$mund,digits=3)
## ---- write workbook ----
writexl::write_xlsx(list(
  OLR              = as.data.frame(summary(olr)$coefficients),
  Shapley_pooled   = sh_pool,
  Shapley_within_realm = sh_realm,
  Shapley_within_basin = sh_basin,
  FEwithin_realm   = rr$fe,  Mundlak_realm = rr$mund,
  FEwithin_basin   = bb$fe,  Mundlak_basin = bb$mund), OUT_FILE)
cat(sprintf("\nWrote: %s\n", OUT_FILE))
