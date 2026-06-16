# =============================================================================
# Transform GCBD bleaching severity (0.1-100% or text codes) into 3 ordered classes (5/30/75)
# Part of: Wind effects on coral bleaching severity (Lapenis & Jiang)
# Language: R
# Inputs : data/input/gcbd_database.xlsx
# Outputs: data/output/gcbd_with_aggregated_severity.csv
# Depends: dplyr, readxl
# Notes  : Paths are set in the CONFIG/USER-SETTINGS block below; place input
#          files in data/input/ and run from the repository root.
#
# Deterministic severity-mapping step of the GLO dataset (base-independent).
# The spatial gridding and oldest-record selection are described in
# docs/glo_dataset_construction.md rather than shipped as code.
# =============================================================================

library(dplyr)
library(readxl)

# --------------------------------------------------------------
# 1. Load your Excel dataset

# --------------------------------------------------------------
input_file <- "data/input/gcbd_database.xlsx"
# Read first sheet; change sheet="name" if needed
df <- read_excel(input_file)

# --------------------------------------------------------------
# 2. Functions to convert severity categories

# --------------------------------------------------------------
severity_from_code <- function(code) {
  case_when(
    grepl("Severe", code, ignore.case = TRUE) ~ 75,
    grepl("Moderate", code, ignore.case = TRUE) ~ 30,
    grepl("Mild", code, ignore.case = TRUE) ~ 5,
    TRUE ~ NA_real_
  )
}
severity_from_percent <- function(p) {
  case_when(
    is.na(p) ~ NA_real_,
    p >= 0.1 & p <= 10 ~ 5,
    p >= 11  & p <= 50 ~ 30,
    p > 50 ~ 75,
    TRUE ~ NA_real_
  )
}

# --------------------------------------------------------------
# 3. Create NEW Aggregated_Severity column

# --------------------------------------------------------------
df <- df %>%
  mutate(
    sev_from_code    = severity_from_code(Severity_Code),
    sev_from_percent = severity_from_percent(Percent_Bleached_Sum),
    Aggregated_Severity = coalesce(sev_from_code, sev_from_percent)
  )

# --------------------------------------------------------------
# 4. Save updated dataset to your directory

# --------------------------------------------------------------
output_file <- "data/output/gcbd_with_aggregated_severity.csv"
write.csv(df, output_file, row.names = FALSE)
cat("SUCCESS: Updated file saved to:\n", output_file, "\n")
