# =============================================================================
# Extract aragonite saturation (Omega) at 0/50/100 m onto reef cells (1-deg grid)
# Part of: Wind effects on coral bleaching severity (Lapenis)
# Language: R
# Inputs : data/input/longterm_climatology_PCA.xlsx, data/input/aragonite_saturation.nc
# Outputs: data/output/climatology_PCA_with_omega.xlsx
# Depends: ncdf4, readxl, writexl
# Notes  : Paths are set in the CONFIG/USER-SETTINGS block below; place input
#          files in data/input/ and run from the repository root.
# =============================================================================

library(ncdf4)
library(readxl)
library(writexl)
xlsx_in  <- "data/input/longterm_climatology_PCA.xlsx"
nc_path  <- "data/input/aragonite_saturation.nc"
xlsx_out <- "data/output/climatology_PCA_with_omega.xlsx"

# -----------------------
# Helpers

# -----------------------
closest_index <- function(grid_vec, x) which.min(abs(grid_vec - x))
# Convert lon to NetCDF convention: 20..380
to_nc_lon_20_380 <- function(lon) {
  lon <- as.numeric(lon)
  lon360 <- ifelse(lon < 0, lon + 360, lon)            # -75 -> 285
  lon380 <- ifelse(lon360 < 20, lon360 + 360, lon360)  # 10 -> 370 (grid starts at 20)
  lon380
}
# Nearest 1-degree index for lon grid 20..380 (361 pts)
lon_to_idx <- function(lon380) {
  idx <- round(lon380) - 20 + 1
  idx <- pmin(pmax(idx, 1), 361)
  as.integer(idx)
}
# Nearest 1-degree index for lat grid -90..90 (181 pts)
lat_to_idx <- function(lat) {
  idx <- round(lat) - (-90) + 1
  idx <- pmin(pmax(idx, 1), 181)
  as.integer(idx)
}

# -----------------------
# Read Excel targets

# -----------------------
df <- read_excel(xlsx_in)
if (!all(c("latitude", "longitude") %in% names(df))) {
  stop("Excel must contain columns named exactly: 'latitude' and 'longitude'.")
}
lat_pts <- as.numeric(df$latitude)
lon_pts <- as.numeric(df$longitude)
lon_pts_nc <- to_nc_lon_20_380(lon_pts)
lon_idx <- lon_to_idx(lon_pts_nc)
lat_idx <- lat_to_idx(lat_pts)
# QA: snapped grid coords
df$omega_grid_lat <- -90 + (lat_idx - 1)
df$omega_grid_lon <-  20 + (lon_idx - 1)

# -----------------------
# Open NetCDF

# -----------------------
nc <- nc_open(nc_path)
on.exit(nc_close(nc), add = TRUE)

# -----------------------
# Get available depth levels robustly
# Depth is stored as a 3-D field in this file, so take unique values.

# -----------------------
depth_all <- as.vector(ncvar_get(nc, "Depth"))
depth_levels <- sort(unique(round(depth_all, 6)))   # rounding avoids tiny float noise
depth_levels <- depth_levels[is.finite(depth_levels)]
# Remove any pathological repeats (sometimes 0 dominates) - keep unique sorted
message("Depth levels found in file: ", paste(depth_levels, collapse = ", "))
target_depths <- c(0, 50, 100)
depth_idx <- vapply(target_depths, function(z) closest_index(depth_levels, z), integer(1))
depth_used <- depth_levels[depth_idx]
message("Requested depths: ", paste(target_depths, collapse = ", "),
        " | Using closest: ", paste(depth_used, collapse = ", "))

# -----------------------
# Determine Aragonite dimension order
# We identify which dim corresponds to:
# - depth (length == length(depth_levels) OR == 9 typically)
# - lon   (length == 361)
# - lat   (length == 181)

# -----------------------
v <- nc$var[["Aragonite"]]
if (is.null(v)) stop("NetCDF variable 'Aragonite' not found.")
dim_lens  <- vapply(v$dim, function(d) d$len, integer(1))
dim_names <- vapply(v$dim, function(d) d$name, character(1))
# Identify lon/lat dims by length (this file is 1-degree global grid)
pos_lon <- which(dim_lens == 361)
pos_lat <- which(dim_lens == 181)
# Identify depth dim: usually 9 in this product
pos_dep <- which(dim_lens %in% c(length(depth_levels), 9))
if (length(pos_lon) != 1 || length(pos_lat) != 1 || length(pos_dep) < 1) {
  stop(paste0(
    "Could not uniquely identify Aragonite dims.\n",
    "Dim names: ", paste(dim_names, collapse = ", "), "\n",
    "Dim lengths: ", paste(dim_lens, collapse = ", ")
  ))
}
pos_dep <- pos_dep[1]  # if multiple candidates, take first
message("Aragonite dim order (name:length): ",
        paste0(dim_names, ":", dim_lens, collapse = " | "))
message("Using positions -> depth:", pos_dep, " lon:", pos_lon, " lat:", pos_lat)

# -----------------------
# Scalar extractor that respects the detected dim order

# -----------------------
get_arag_scalar <- function(depth_i, lon_i, lat_i) {
  start <- rep(1L, length(dim_lens))
  count <- rep(1L, length(dim_lens))
  start[pos_dep] <- as.integer(depth_i)
  start[pos_lon] <- as.integer(lon_i)
  start[pos_lat] <- as.integer(lat_i)
  as.numeric(ncvar_get(nc, "Aragonite",
                      start = as.integer(start),
                      count = as.integer(count)))[1]
}
# NOTE: depth_idx we computed is within depth_levels (unique depths),
# but the NetCDF depth index is 1..9. In this product, depth ordering is:
# 0, 50, 100, 200, 500, 1000, 2000, 3000, 4000 (typically).
# So we map requested depths to the nearest of that expected sequence by index.
#
# Robust approach: build the expected 9-level vector from the file itself:
# Take one vertical profile of Depth at any lon/lat and use that as the depth axis.
# We'll do that now using the Depth variable's Aragonite-like dimension order.
depth_var <- nc$var[["Depth"]]
if (is.null(depth_var)) stop("NetCDF variable 'Depth' not found.")
d_dim_lens <- vapply(depth_var$dim, function(d) d$len, integer(1))
# Create start/count to read a vertical profile: vary the depth dimension, fix lon/lat to 1
# Identify depth dimension in Depth var as the one with length 9
d_pos_dep <- which(d_dim_lens == 9)
if (length(d_pos_dep) != 1) {
  # fallback: use the dimension that matches one of Aragonite's depth candidates
  d_pos_dep <- which(d_dim_lens %in% c(9, length(depth_levels)))[1]
}
start_d <- rep(1L, length(d_dim_lens))
count_d <- rep(1L, length(d_dim_lens))
count_d[d_pos_dep] <- d_dim_lens[d_pos_dep]
depth_profile <- as.vector(ncvar_get(nc, "Depth", start = start_d, count = count_d))
depth_profile <- round(depth_profile, 6)
message("Depth profile (axis) read from file: ", paste(depth_profile, collapse = ", "))
# Now map target depths to indices of this depth_profile (this is the actual depth index 1..9)
depth_i0   <- closest_index(depth_profile, 0)
depth_i50  <- closest_index(depth_profile, 50)
depth_i100 <- closest_index(depth_profile, 100)
message("Depth indices used (in NetCDF): 0m=", depth_i0,
        " 50m=", depth_i50, " 100m=", depth_i100)

# -----------------------
# Extract omega at the three depths

# -----------------------
n <- nrow(df)
omega_0m   <- numeric(n)
omega_50m  <- numeric(n)
omega_100m <- numeric(n)
for (i in seq_len(n)) {
  omega_0m[i]   <- get_arag_scalar(depth_i0,   lon_idx[i], lat_idx[i])
  omega_50m[i]  <- get_arag_scalar(depth_i50,  lon_idx[i], lat_idx[i])
  omega_100m[i] <- get_arag_scalar(depth_i100, lon_idx[i], lat_idx[i])
  if (i %% 200 == 0) message("Processed ", i, " / ", n)
}
df$omega_0m   <- omega_0m
df$omega_50m  <- omega_50m
df$omega_100m <- omega_100m
message("Sanity summaries:")
message("omega_0m:");   print(summary(df$omega_0m))
message("omega_50m:");  print(summary(df$omega_50m))
message("omega_100m:"); print(summary(df$omega_100m))
write_xlsx(df, xlsx_out)
message("Wrote: ", xlsx_out)
