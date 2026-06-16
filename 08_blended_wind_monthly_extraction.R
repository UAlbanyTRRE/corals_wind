# =============================================================================
# Monthly mean wind speed (1991-2010) interpolated to reef coordinates (NOAA blended winds)
# Part of: Wind effects on coral bleaching severity (Lapenis & Jiang)
# Language: R
# Inputs : NOAA blended-wind daily field + reef coordinates (supplied by the caller)
# Outputs: data/output/windspeed_1991_2010_monthly.xlsx
# Depends: openxlsx
# Notes  : Run from the repository root; output lands in data/output/.
#
# WHAT THIS FILE IS
#   A self-contained HELPER MODULE. Sourcing it only *defines* two functions and
#   runs nothing, so source("R/08_...R") never errors. The heavy extraction is
#   performed only when you call extract_monthly_mean_wind(...) with a daily
#   NOAA blended-wind field that you have loaded (the field itself is a public
#   NOAA product, not redistributed here; see the Data Availability statement).
#
#   This replaces the earlier partial fragment, which referenced undefined
#   objects (results_df, interpolate_windspeed, windspeed, latitudes,
#   longitudes, dayofyear) and could not run standalone.
#
# USAGE
#   source("R/08_blended_wind_monthly_extraction.R")
#   # load your daily blended-wind field and reef table, then:
#   extract_monthly_mean_wind(
#     results_df  = reef_table,        # data.frame with LATITUDE, LONGITUDE
#     windspeed   = ws_array,          # daily field, see interpolate_windspeed()
#     latitudes   = lat_vec,           # ascending grid latitudes
#     longitudes  = lon_vec,           # ascending grid longitudes
#     dayofyear   = doy_vec)           # day index matching dim 3 of windspeed
# =============================================================================

# interpolate_windspeed() below needs no extra packages; the workbook writer
# (extract_monthly_mean_wind) checks for openxlsx at call time, so simply
# sourcing this file only defines functions and never errors.

# -----------------------------------------------------------------------------
# Bilinear interpolation of a daily wind field to an arbitrary (lat, lon).
#
# ASSUMED LAYOUT (confirm against your blended-wind file and adjust if needed):
#   - `windspeed`  is a 3-D array indexed [longitude, latitude, day]
#   - `longitudes` and `latitudes` are the matching coordinate vectors, ascending
#   - returns a numeric vector of length dim(windspeed)[3] (one value per day)
#
# Longitudes are wrapped into the grid's range so that, e.g., a -75 deg reef can
# be matched against a 0..360 grid. Points outside the grid are clamped to the
# nearest edge cell.
# -----------------------------------------------------------------------------
interpolate_windspeed <- function(lat, lon, windspeed, latitudes, longitudes) {
  nlon <- length(longitudes); nlat <- length(latitudes)
  if (nlon < 2L || nlat < 2L)
    stop("interpolate_windspeed(): need at least 2 grid points per axis.", call. = FALSE)

  # wrap longitude into the grid's convention (handles 0..360 vs -180..180 grids)
  if (lon < min(longitudes)) lon <- lon + 360
  if (lon > max(longitudes)) lon <- lon - 360

  # bracketing indices (all.inside keeps jx+1 / jy+1 valid at the upper edge)
  jx <- findInterval(lon, longitudes, all.inside = TRUE)
  jy <- findInterval(lat, latitudes,  all.inside = TRUE)

  x1 <- longitudes[jx]; x2 <- longitudes[jx + 1L]
  y1 <- latitudes[jy];  y2 <- latitudes[jy + 1L]
  tx <- if (x2 > x1) (lon - x1) / (x2 - x1) else 0
  ty <- if (y2 > y1) (lat - y1) / (y2 - y1) else 0
  tx <- min(max(tx, 0), 1); ty <- min(max(ty, 0), 1)

  # four corner daily series
  q11 <- windspeed[jx,      jy,      ]
  q21 <- windspeed[jx + 1L, jy,      ]
  q12 <- windspeed[jx,      jy + 1L, ]
  q22 <- windspeed[jx + 1L, jy + 1L, ]

  (1 - tx) * (1 - ty) * q11 + tx * (1 - ty) * q21 +
    (1 - tx) * ty * q12 + tx * ty * q22
}

# -----------------------------------------------------------------------------
# Compute the monthly mean wind speed at every reef site, for each year, and
# write one sheet per year to an .xlsx workbook.
#
# Months are formed from fixed 30-day blocks of `dayofyear` (block m =
# days (m-1)*30+1 .. min(m*30, length(dayofyear))), matching the original
# workflow. Pass a different `years` vector to change the span.
# -----------------------------------------------------------------------------
extract_monthly_mean_wind <- function(results_df, windspeed, latitudes, longitudes,
                                      dayofyear, years = 1991:2010,
                                      out_file = "data/output/windspeed_1991_2010_monthly.xlsx") {
  stopifnot(is.data.frame(results_df),
            all(c("LATITUDE", "LONGITUDE") %in% names(results_df)))
  if (!requireNamespace("openxlsx", quietly = TRUE))
    stop("Package 'openxlsx' is required to write the workbook. Install it with ",
         "install.packages(\"openxlsx\"), or restore the environment with renv::restore().",
         call. = FALSE)
  dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)

  total_rows <- nrow(results_df)
  # Pre-compute the day ranges for each of the 12 months (fixed 30-day blocks).
  day_ranges <- lapply(1:12, function(month) {
    start_day <- (month - 1) * 30 + 1
    end_day   <- min(month * 30, length(dayofyear))
    start_day:end_day
  })

  monthly_mean_results <- list()
  for (year in years) {
    year_windspeed <- matrix(NA_real_, nrow = total_rows, ncol = 12)
    cat(sprintf("Processing Year: %d\n", year))

    interpolated_windspeed_list <- lapply(seq_len(total_rows), function(i) {
      interpolate_windspeed(results_df$LATITUDE[i], results_df$LONGITUDE[i],
                            windspeed, latitudes, longitudes)
    })

    for (i in seq_len(total_rows)) {
      interpolated_windspeed <- interpolated_windspeed_list[[i]]
      if (!all(is.na(interpolated_windspeed))) {
        for (month in 1:12) {
          monthly_data <- interpolated_windspeed[day_ranges[[month]]]
          year_windspeed[i, month] <- mean(monthly_data, na.rm = TRUE)
        }
      }
      if (i %% 100 == 0)
        cat(sprintf("Year %d - %.2f%% complete\n", year, (i / total_rows) * 100))
    }
    monthly_mean_results[[as.character(year)]] <- year_windspeed
  }

  wb <- openxlsx::createWorkbook()
  month_names <- c("January", "February", "March", "April", "May", "June",
                   "July", "August", "September", "October", "November", "December")
  for (year in years) {
    openxlsx::addWorksheet(wb, as.character(year))
    sheet_data <- cbind(results_df[, c("LATITUDE", "LONGITUDE")],
                        monthly_mean_results[[as.character(year)]])
    colnames(sheet_data) <- c("LATITUDE", "LONGITUDE", month_names)
    openxlsx::writeData(wb, sheet = as.character(year), sheet_data)
  }
  openxlsx::saveWorkbook(wb, out_file, overwrite = TRUE)
  cat(sprintf("Monthly mean wind speeds for %d-%d written to %s\n",
              min(years), max(years), out_file))
  invisible(out_file)
}
