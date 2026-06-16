# =============================================================================
# Box Upwelling-Diffusion (BUD) carbonate model; U6 normalization; full m/n/q grid
# Part of: Wind effects on coral bleaching severity (Lapenis & Jiang)
# Language: R
# Inputs : (none - self-contained)
# Outputs: data/output/box_diffusion_mnq_withFa_U6.xlsx
# Depends: openxlsx, seacarb
# Notes  : Paths are set in the CONFIG/USER-SETTINGS block below; place input
#          files in data/input/ and run from the repository root.
# =============================================================================

# ================================
# BUD Box Model - U6 normalization, full m/n/q grid
# Output: data/output/box_diffusion_mnq_withFa_U6.xlsx
# Deps: seacarb, openxlsx

# ================================
suppressPackageStartupMessages({
  .need <- c("seacarb", "openxlsx")
  .miss <- .need[!vapply(.need, requireNamespace, logical(1), quietly = TRUE)]
  if (length(.miss) > 0)
    stop("Missing required package(s): ", paste(.miss, collapse = ", "),
         "\n  Install them, or restore the recorded environment with renv::restore().",
         call. = FALSE)
})
library(seacarb)
library(openxlsx)

# -------- Constants --------
rho       <- 1025          # kg m^-3
h         <- 40            # m
S         <- 35
T_C       <- 25            # deg C
TA_DEEP   <- 2350e-6       # mol kg^-1
DIC_DEEP  <- 2100e-6       # mol kg^-1
Ca_DEEP   <- 10.3e-3       # mol kg^-1
# Reference wind speed changed from U5 to U6
U6        <- 6             # m s^-1
W0        <- 4.7e-7        # m s^-1 at U = 6 m/s
Az0       <- 1e-5          # m^2 s^-1 at U = 6 m/s
B0_mol_s  <-32.4e-9        # mol C m^-2 s^-1 at U = 6 m/s
rcaco3    <- 0.07          # fraction as CaCO3
rN        <- 16/106        # mol alk per mol C org production
alpha     <- 6.97e-7       # m s^-1 per U^2, Wanninkhof 2014, Sc=660
Pa_atm    <- 430e-6        # atm
sec_yr    <- 365*24*3600
gC        <- 12.0

# -------- CO2 solubility, Weiss 1974 --------
K0_CO2_Weiss <- function(T, S){
  TK <- T + 273.15
  A1 <- -58.0931; A2 <- 90.5069; A3 <- 22.2940
  B1 <-  0.027766; B2 <- -0.025888; B3 <- 0.0050578
  lnK0 <- A1 + A2*(100/TK) + A3*log(TK/100) +
    S*(B1 + B2*(TK/100) + B3*(TK/100)^2)
  exp(lnK0)
}
K0 <- K0_CO2_Weiss(T_C, S)

# -------- Parameter grid --------
U_vals <- 1:10
m_set  <- c(0, 1, 2)
n_set  <- c(0, 1, 2)
q_set  <- c(0, 1, 2)

# -------- Helper functions --------
mix_rate <- function(U, m, q){
  W   <- W0  * (U/U6)^m
  Az  <- Az0 * (U/U6)^q
  M   <- W + Az/h
  list(W = W, Az = Az, M = M)
}
bio_fluxes <- function(U, n){
  scale <- (U/U6)^n
  Btot  <- B0_mol_s * scale
  # signed surface tendencies: negative = removal from mixed ocean layer
  F_org   <- -(1 - rcaco3) * Btot
  F_caco3 <- -rcaco3 * Btot
  # positive organic production for nitrate alkalinity term
  F_org_prod <- -F_org
  list(
    Btot = Btot,
    F_org = F_org,
    F_caco3 = F_caco3,
    F_org_prod = F_org_prod
  )
}
Fa_from_Pw <- function(U, Pw_atm){
  k <- alpha * U^2
  rho * k * K0 * (Pa_atm - Pw_atm)
}
solve_Pw <- function(U, M, F_org, F_caco3, F_org_prod, tol = 1e-7){
  gfun <- function(Pw){
    Fa <- Fa_from_Pw(U, Pw)
    DIC_MOL <- DIC_DEEP + (F_org + F_caco3 + Fa) / (rho * M)
    Ca_MOL  <- Ca_DEEP  + F_caco3 / (rho * M)
    TA_MOL  <- TA_DEEP  + (2 * F_caco3 + rN * F_org_prod) / (rho * M)
    co <- carb(
      flag = 15,
      var1 = TA_MOL,
      var2 = DIC_MOL,
      S = S,
      T = T_C,
      Patm = 1,
      P = 0
    )
    co$pCO2 / 1e6 - Pw
  }
  lo <- 200e-6
  hi <- 800e-6
  glo <- tryCatch(gfun(lo), error = function(e) NA_real_)
  ghi <- tryCatch(gfun(hi), error = function(e) NA_real_)
  if (is.finite(glo) && is.finite(ghi) && glo * ghi < 0){
    uniroot(gfun, interval = c(lo, hi), tol = tol)$root
  } else {
    Pw <- 450e-6
    for (it in 1:200){
      Fa <- Fa_from_Pw(U, Pw)
      DIC_MOL <- DIC_DEEP + (F_org + F_caco3 + Fa) / (rho * M)
      Ca_MOL  <- Ca_DEEP  + F_caco3 / (rho * M)
      TA_MOL  <- TA_DEEP  + (2 * F_caco3 + rN * F_org_prod) / (rho * M)
      co <- carb(
        flag = 15,
        var1 = TA_MOL,
        var2 = DIC_MOL,
        S = S,
        T = T_C,
        Patm = 1,
        P = 0
      )
      Pw_new <- max(lo, min(hi, co$pCO2 / 1e6))
      Pw <- 0.35 * Pw_new + 0.65 * Pw
      if (abs(Pw_new - Pw) < tol) break
    }
    Pw
  }
}
run_case <- function(m, n, q){
  rows <- vector("list", length(U_vals))
  for (i in seq_along(U_vals)){
    U <- U_vals[i]
    mr <- mix_rate(U, m, q)
    W  <- mr$W
    Az <- mr$Az
    M  <- mr$M
    bf <- bio_fluxes(U, n)
    F_org       <- bf$F_org
    F_caco3     <- bf$F_caco3
    F_org_prod  <- bf$F_org_prod
    Pw <- tryCatch(
      solve_Pw(U, M, F_org, F_caco3, F_org_prod),
      error = function(e) 450e-6
    )
    Fa <- Fa_from_Pw(U, Pw)
    DIC_MOL <- DIC_DEEP + (F_org + F_caco3 + Fa) / (rho * M)
    Ca_MOL  <- Ca_DEEP  + F_caco3 / (rho * M)
    TA_MOL  <- TA_DEEP  + (2 * F_caco3 + rN * F_org_prod) / (rho * M)
    co <- try(
      carb(
        flag = 15,
        var1 = TA_MOL,
        var2 = DIC_MOL,
        S = S,
        T = T_C,
        Patm = 1,
        P = 0
      ),
      silent = TRUE
    )
    if (inherits(co, "try-error")) {
      co <- list(pCO2 = NA, pH = NA, OmegaAragonite = NA)
    }
    Fa_gC  <- Fa * gC * sec_yr
    Fup_gC <- rho * W * (DIC_DEEP - DIC_MOL) * gC * sec_yr
    Fd_gC  <- rho * (Az/h) * (DIC_DEEP - DIC_MOL) * gC * sec_yr
    B_gC   <- (-(F_org + F_caco3)) * gC * sec_yr
    Net_gC <- Fa_gC + Fup_gC + Fd_gC - B_gC
    rows[[i]] <- data.frame(
      U = U,
      Uref = U6,
      m = m,
      n = n,
      q = q,
      W = W,
      Az = Az,
      mix_velocity = M,
      DICMOL = DIC_MOL,
      CaMOL = Ca_MOL,
      TAMOL = TA_MOL,
      PwCO2_uatm = co$pCO2,
      pH = co$pH,
      OmegaArag = co$OmegaAragonite,
      Fa_gC = Fa_gC,
      Fup_gC = Fup_gC,
      Fd_gC = Fd_gC,
      B_gC = B_gC,
      Net_gC = Net_gC,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

# -------- Driver --------
out_dir <- "data/output"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
out_file <- file.path(out_dir, "box_diffusion_mnq_withFa_U6.xlsx")
wb <- createWorkbook()
for (m in m_set){
  for (n in n_set){
    for (q in q_set){
      cat(sprintf("Running m=%g n=%g q=%g with U6 normalization...\n", m, n, q))
      res <- run_case(m, n, q)
      sh <- sprintf("m%g_n%g_q%g", m, n, q)
      addWorksheet(wb, sh)
      writeData(wb, sh, res)
    }
  }
}
saveWorkbook(wb, out_file, overwrite = TRUE)
cat(sprintf("Done. Wrote '%s'.\n", out_file))
