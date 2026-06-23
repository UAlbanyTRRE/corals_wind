# =============================================================================
# Tropical-cyclone power (TCP/PDI) within 400 km via Holland radial-wind profile (IBTrACS v04r01)
# Part of: Wind effects on coral bleaching severity (Lapenis)
# Language: R
# Inputs : data/input/saf_event_coordinates_dates.xlsx, data/input/ibtracs.ALL.list.v04r01.csv
# Outputs: data/output/tcp_pdi_400km_6m_12m_1993_2020.xlsx (+ TCpower attenuated-Rmax workbook)
# Depends: data.table, readxl, writexl, lubridate, geosphere, parallel, pbapply
# Notes  : Paths are set in the CONFIG/USER-SETTINGS block below; place input
#          files in data/input/ and run from the repository root.
#          Radial wind profile follows Holland (1980) (B-parameter, Coriolis-free
#          proxy), matching Supplementary Note 2.
# =============================================================================

# ============================================================
# PDI within 400 km of each site
# - 6 months prior to event date
# - 12 months prior to event date
# - total over 1993-2020
# IBTrACS v04r01 CSV
# Optimized: data.table + time key + bbox prefilter + parallel + progress

# ============================================================
suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(writexl)
  library(lubridate)
  library(geosphere)
  library(parallel)
  library(pbapply)
})

# ---------------------------
# User inputs

# ---------------------------
events_xlsx <- "data/input/saf_event_coordinates_dates.xlsx"
ibtracs_csv <- "data/input/ibtracs.ALL.list.v04r01.csv"
out_xlsx <- "data/output/tcp_pdi_400km_6m_12m_1993_2020.xlsx"
R_km <- 400
R_m  <- R_km * 1000
start_year <- 1993
end_year   <- 2020
KT_TO_MS <- 0.514444
# Parallel: keep 1 core free
ncores <- max(1, detectCores() - 1)
pboptions(type = "timer")

# ---------------------------
# Helpers

# ---------------------------
pick_first_existing <- function(dt, candidates) {
  hit <- candidates[candidates %in% names(dt)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}
pick_wind_col <- function(names_vec) {
  # Prefer USA_WIND if available; otherwise common alternatives
  candidates <- c("USA_WIND", "WMO_WIND", "TOKYO_WIND", "REUNION_WIND",
                  "BOM_WIND", "NEWDELHI_WIND", "HKO_WIND", "CMA_WIND")
  hit <- candidates[candidates %in% names_vec]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}
to_Date_robust <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (inherits(x, "POSIXct") || inherits(x, "POSIXt")) return(as.Date(x))
  if (is.numeric(x)) return(as.Date(x, origin = "1899-12-30"))
  x_chr <- trimws(as.character(x))
  x_chr[x_chr == "" | tolower(x_chr) %in% c("na", "null")] <- NA_character_
  parsed <- suppressWarnings(parse_date_time(
    x_chr,
    orders = c(
      "Ymd", "Y-m-d", "Y/m/d",
      "mdY", "m/d/Y", "m-d-Y",
      "dmY", "d/m/Y", "d-m-Y",
      "Ymd HMS", "Y-m-d H:M:S", "Y/m/d H:M:S"
    ),
    tz = "UTC"
  ))
  as.Date(parsed)
}
to_posix_safe <- function(x) {
  if (inherits(x, "POSIXct")) return(x)
  if (inherits(x, "Date")) return(as.POSIXct(x, tz = "UTC"))
  suppressWarnings({
    p <- ymd_hms(x, tz = "UTC", quiet = TRUE)
    if (all(is.na(p))) p <- ymd_hm(x, tz = "UTC", quiet = TRUE)
    if (all(is.na(p))) p <- ymd(x, tz = "UTC", quiet = TRUE)
    p
  })
}
bbox_filter <- function(dt, lat0, lon0, km) {
  dlat <- km / 111.32
  dlon <- km / (111.32 * max(cos(lat0 * pi/180), 0.01))
  dt[lat >= (lat0 - dlat) & lat <= (lat0 + dlat) &
       lon >= (lon0 - dlon) & lon <= (lon0 + dlon)]
}
# Core worker for one site & one time window
cum_pdi_site_window <- function(tracks, lat0, lon0, t_start, t_end, R_m, R_km) {
  sub <- tracks[time_utc >= t_start & time_utc <= t_end]
  if (nrow(sub) == 0) return(0)
  # fast bbox prefilter
  sub <- bbox_filter(sub, lat0, lon0, R_km)
  if (nrow(sub) == 0) return(0)
  d <- distHaversine(matrix(c(sub$lon, sub$lat), ncol = 2), c(lon0, lat0))
  idx <- which(d <= R_m)
  if (length(idx) == 0) return(0)
  # PDI sum: V^3 * dt_seconds
  sum((sub$wind_ms[idx]^3) * sub$dt_s[idx], na.rm = TRUE)
}

# ---------------------------
# Read events

# ---------------------------
message("Reading events: ", events_xlsx)
events <- as.data.table(read_excel(events_xlsx))
# Explicit mapping (your column names)
events[, lat := as.numeric(Latitude_Degree)]
events[, lon := as.numeric(Longitude_Degree)]
events[, date := to_Date_robust(Date)]
events[, event_id := .I]
if (anyNA(events$lat) || anyNA(events$lon) || anyNA(events$date)) {
  bad <- events[is.na(lat) | is.na(lon) | is.na(date)]
  message("Some rows have missing/invalid lat/lon/date. First few:")
  print(head(bad, 10))
  stop("Fix missing/invalid values in the events file, then rerun.")
}

# ---------------------------
# Read IBTrACS

# ---------------------------
message("Reading IBTrACS (can be large): ", ibtracs_csv)
ib <- fread(ibtracs_csv, showProgress = TRUE)
col_time <- pick_first_existing(ib, c("ISO_TIME", "iso_time", "time", "TIME"))
col_sid  <- pick_first_existing(ib, c("SID", "sid", "STORMID", "storm_id"))
col_lat  <- pick_first_existing(ib, c("LAT", "lat", "Latitude"))
col_lon  <- pick_first_existing(ib, c("LON", "lon", "Longitude"))
wind_col <- pick_wind_col(names(ib))
need <- c(col_time, col_sid, col_lat, col_lon, wind_col)
if (any(is.na(need))) {
  stop("Could not find required columns in IBTrACS. Needed: ISO_TIME, SID, LAT, LON, and a WIND column.\n",
       "Found columns include: ", paste(names(ib), collapse = ", "))
}
tracks <- ib[, .(
  SID      = get(col_sid),
  time_raw = get(col_time),
  lat      = suppressWarnings(as.numeric(get(col_lat))),
  lon      = suppressWarnings(as.numeric(get(col_lon))),
  wind_kt  = suppressWarnings(as.numeric(get(wind_col)))
)]
tracks[, time_utc := to_posix_safe(time_raw)]
tracks <- tracks[!is.na(time_utc) & !is.na(lat) & !is.na(lon) & !is.na(wind_kt)]
# Limit to years of interest
tracks <- tracks[year(time_utc) >= start_year & year(time_utc) <= end_year]
# Wind -> m/s
tracks[, wind_ms := wind_kt * KT_TO_MS]
tracks <- tracks[is.finite(wind_ms) & wind_ms > 0]
if (nrow(tracks) == 0) stop("No usable track points after filtering. Check wind/time columns.")
# Compute dt_s per trackpoint within storm
setorder(tracks, SID, time_utc)
tracks[, dt_s := as.numeric(difftime(shift(time_utc, type = "lead"),
                                    time_utc, units = "secs")), by = SID]
# Clean dt: if missing/invalid/too large, default to 6 hours
tracks[is.na(dt_s) | dt_s <= 0 | dt_s > 24*3600, dt_s := 6*3600]
# Key by time for fast slicing
setkey(tracks, time_utc)

# ---------------------------
# Precompute FULL-PERIOD PDI for unique sites

# ---------------------------
sites <- unique(events[, .(lat, lon)])
sites[, site_id := .I]
full_start <- as.POSIXct(sprintf("%d-01-01 00:00:00", start_year), tz = "UTC")
full_end   <- as.POSIXct(sprintf("%d-12-31 23:59:59", end_year), tz = "UTC")
message(sprintf("Computing FULL-PERIOD PDI (%d-%d) within %dkm for %d unique sites using %d cores...",
                start_year, end_year, R_km, nrow(sites), ncores))
site_chunks <- split(seq_len(nrow(sites)), cut(seq_len(nrow(sites)), breaks = ncores, labels = FALSE))
full_list <- pblapply(site_chunks, function(idx) {
  s <- sites[idx]
  out <- numeric(nrow(s))
  for (i in seq_len(nrow(s))) {
    out[i] <- cum_pdi_site_window(tracks, s$lat[i], s$lon[i], full_start, full_end, R_m, R_km)
  }
  out
}, cl = ncores)
sites[, PDI_1993_2020_400km := unlist(full_list)]

# ---------------------------
# Compute 6m and 12m windows per event (parallel across events)

# ---------------------------
events <- merge(events, sites, by = c("lat", "lon"), all.x = TRUE)
events_idx <- split(seq_len(nrow(events)), cut(seq_len(nrow(events)), breaks = ncores, labels = FALSE))
message(sprintf("Computing 6-month & 12-month prior-to-event PDI within %dkm for %d events using %d cores...",
                R_km, nrow(events), ncores))
event_results <- pblapply(events_idx, function(idx) {
  sub <- events[idx, .(event_id, lat, lon, date)]
  sub[, `:=`(
    start_6m  = as.POSIXct(date %m-% months(6),  tz = "UTC"),
    start_12m = as.POSIXct(date %m-% months(12), tz = "UTC"),
    end_t     = as.POSIXct(paste0(date, " 23:59:59"), tz = "UTC")
  )]
  out6  <- numeric(nrow(sub))
  out12 <- numeric(nrow(sub))
  for (i in seq_len(nrow(sub))) {
    out6[i] <- cum_pdi_site_window(tracks, sub$lat[i], sub$lon[i],
                                   sub$start_6m[i], sub$end_t[i], R_m, R_km)
    out12[i] <- cum_pdi_site_window(tracks, sub$lat[i], sub$lon[i],
                                    sub$start_12m[i], sub$end_t[i], R_m, R_km)
  }
  data.table(event_id = sub$event_id,
             PDI_6m_400km  = out6,
             PDI_12m_400km = out12)
}, cl = ncores)
event_results_dt <- rbindlist(event_results)
setkey(event_results_dt, event_id)
setkey(events, event_id)
events <- event_results_dt[events]

# ---------------------------
# Output

# ---------------------------
out <- events[, .(
  Latitude_Degree  = lat,
  Longitude_Degree = lon,
  Date             = date,
  PDI_6m_400km      = PDI_6m_400km,
  PDI_12m_400km     = PDI_12m_400km,
  PDI_1993_2020_400km = PDI_1993_2020_400km
)]
write_xlsx(out, out_xlsx)
message("DONE. Saved: ", out_xlsx)

# ============================================================
# TROPICAL CYCLONES BY HOLLAND MODEL AT GIVEN SITE

# ============================================================
# Tropical Cyclone Cumulative Power at Sites (attenuated winds)
#  - Radius of influence: 400 km
#  - 6 months before event date
#  - 12 months before event date
#  - Full period 1993-2020
#
# Wind at the site is attenuated with distance using a Holland-type
# radial wind profile, using Vmax at the track point and Rmax.
#
# Rmax handling:
#   - If IBTrACS RMW/RMAX is available and valid, use it.
#   - Otherwise Option B intensity-bin defaults:
#       TS/weak (<64 kt): 60 km
#       Cat 1-2 (64-95 kt): 40 km
#       Cat 3-5 (>=96 kt): 23 km
#
# Output "power" proxy: sum( V_site^3 * dt ), dt in seconds

# ============================================================
suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(writexl)
  library(lubridate)
  library(geosphere)
  library(parallel)
  library(pbapply)
})

# ---------------------------
# USER PATHS

# ---------------------------
events_xlsx <- "data/input/saf_event_coordinates_dates.xlsx"
ibtracs_csv <- "data/input/ibtracs.ALL.list.v04r01.csv"
out_xlsx    <- "data/output/tcp_power_400km_6m_12m_attenuatedRmaxBins.xlsx"

# ---------------------------
# SETTINGS

# ---------------------------
R_km <- 400
R_m  <- R_km * 1000
start_year <- 1993
end_year   <- 2020
KT_TO_MS <- 0.514444
# Holland B (kept constant here; you can later swap a variable B if you want)
B_default <- 1.5
# Option B defaults (km)
Rmax_TS_km   <- 60
Rmax_C12_km  <- 40
Rmax_C35_km  <- 23
# Parallel cores (leave 1 core free)
ncores <- max(1, detectCores() - 1)
pboptions(type = "timer")

# ---------------------------
# HELPERS

# ---------------------------
pick_first_existing <- function(dt, candidates) {
  hit <- candidates[candidates %in% names(dt)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}
to_date_safe <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (inherits(x, "POSIXct") || inherits(x, "POSIXt")) return(as.Date(x))
  if (is.numeric(x)) return(as.Date(x, origin = "1899-12-30"))  # Excel Windows origin
  x_chr <- trimws(as.character(x))
  x_chr[x_chr == "" | tolower(x_chr) %in% c("na", "null")] <- NA_character_
  suppressWarnings(as.Date(parse_date_time(
    x_chr,
    orders = c("Ymd", "Y-m-d", "Y/m/d", "mdY", "m/d/Y", "dmY", "d/m/Y", "Ymd HMS", "Y-m-d H:M:S"),
    tz = "UTC"
  )))
}
to_posix_safe <- function(x) {
  if (inherits(x, "POSIXct")) return(x)
  if (inherits(x, "Date")) return(as.POSIXct(x, tz = "UTC"))
  suppressWarnings({
    p <- ymd_hms(x, tz = "UTC", quiet = TRUE)
    if (all(is.na(p))) p <- ymd_hm(x, tz = "UTC", quiet = TRUE)
    if (all(is.na(p))) p <- ymd(x, tz = "UTC", quiet = TRUE)
    p
  })
}
# Holland-type radial wind profile (fast, coriolis-free proxy)
holland_wind_ms <- function(r_m, Rmax_m, Vmax_ms, B = B_default) {
  r_m    <- pmax(r_m, 1)       # avoid 0
  Rmax_m <- pmax(Rmax_m, 1000) # avoid too small
  x <- (Rmax_m / r_m)^B
  Vmax_ms * sqrt(x * exp(1 - x))
}
# Bounding box prefilter (fast)
bbox_filter <- function(dt, lat0, lon0, km) {
  dlat <- km / 111.32
  dlon <- km / (111.32 * max(cos(lat0 * pi/180), 0.01))
  dt[lat >= (lat0 - dlat) & lat <= (lat0 + dlat) &
       lon >= (lon0 - dlon) & lon <= (lon0 + dlon)]
}
# Compute site power in a window for many sites (used for full-period unique sites)
compute_power_for_window_manysites <- function(tracks_dt, sites_dt, win_start, win_end, R_m, R_km) {
  tr <- tracks_dt[time_utc >= win_start & time_utc <= win_end]
  if (nrow(tr) == 0) return(rep(0, nrow(sites_dt)))
  out <- numeric(nrow(sites_dt))
  for (i in seq_len(nrow(sites_dt))) {
    lat0 <- sites_dt$lat[i]
    lon0 <- sites_dt$lon[i]
    sub <- bbox_filter(tr, lat0, lon0, R_km)
    if (nrow(sub) == 0) { out[i] <- 0; next }
    d <- distHaversine(matrix(c(sub$lon, sub$lat), ncol = 2), c(lon0, lat0))
    idx <- which(d <= R_m)
    if (length(idx) == 0) { out[i] <- 0; next }
    w <- holland_wind_ms(r_m = d[idx], Rmax_m = sub$Rmax_m[idx], Vmax_ms = sub$Vmax_ms[idx], B = sub$B[idx])
    out[i] <- sum((w^3) * sub$dt_s[idx], na.rm = TRUE)
  }
  out
}
# Single-site convenience for event windows
compute_power_for_window_onesite <- function(tracks_dt, lat0, lon0, win_start, win_end, R_m, R_km) {
  tr <- tracks_dt[time_utc >= win_start & time_utc <= win_end]
  if (nrow(tr) == 0) return(0)
  sub <- bbox_filter(tr, lat0, lon0, R_km)
  if (nrow(sub) == 0) return(0)
  d <- distHaversine(matrix(c(sub$lon, sub$lat), ncol = 2), c(lon0, lat0))
  idx <- which(d <= R_m)
  if (length(idx) == 0) return(0)
  w <- holland_wind_ms(r_m = d[idx], Rmax_m = sub$Rmax_m[idx], Vmax_ms = sub$Vmax_ms[idx], B = sub$B[idx])
  sum((w^3) * sub$dt_s[idx], na.rm = TRUE)
}

# ---------------------------
# READ EVENTS (explicit mapping)

# ---------------------------
message("Reading events: ", events_xlsx)
events <- as.data.table(read_excel(events_xlsx))
# Your known column names:
events[, lat  := as.numeric(Latitude_Degree)]
events[, lon  := as.numeric(Longitude_Degree)]
events[, date := to_date_safe(Date)]
if (anyNA(events$lat) || anyNA(events$lon) || anyNA(events$date)) {
  bad <- events[is.na(lat) | is.na(lon) | is.na(date)]
  message("Some rows have missing/invalid lat/lon/date. First few shown:")
  print(head(bad, 10))
  stop("Fix missing/invalid values in the events file, then rerun.")
}
events[, event_id := .I]

# ---------------------------
# READ IBTRACS

# ---------------------------
message("Reading IBTrACS (can be large): ", ibtracs_csv)
ib <- fread(ibtracs_csv, showProgress = TRUE)
col_time <- pick_first_existing(ib, c("ISO_TIME", "iso_time", "TIME", "time"))
col_sid  <- pick_first_existing(ib, c("SID", "sid"))
col_lat  <- pick_first_existing(ib, c("LAT", "lat", "Latitude"))
col_lon  <- pick_first_existing(ib, c("LON", "lon", "Longitude"))
# Prefer USA_WIND then WMO_WIND (or vice versa if you want)
col_wind <- pick_first_existing(ib, c("USA_WIND", "WMO_WIND", "TOKYO_WIND", "REUNION_WIND",
                                      "BOM_WIND", "NEWDELHI_WIND", "HKO_WIND", "CMA_WIND"))
# Radius of max winds candidates (often incomplete)
col_rmax <- pick_first_existing(ib, c("USA_RMW", "WMO_RMW", "RMW", "RMW_KM", "RMW_NM", "USA_RMW_KM"))
need <- c(col_time, col_sid, col_lat, col_lon, col_wind)
if (any(is.na(need))) {
  stop("Could not find required IBTrACS columns. Required: time, SID, LAT, LON, WIND.\n",
       "Columns present: ", paste(names(ib), collapse = ", "))
}
tracks <- ib[, .(
  SID      = get(col_sid),
  time_raw = get(col_time),
  lat      = suppressWarnings(as.numeric(get(col_lat))),
  lon      = suppressWarnings(as.numeric(get(col_lon))),
  wind_kt  = suppressWarnings(as.numeric(get(col_wind))),
  rmax_raw = if (!is.na(col_rmax)) suppressWarnings(as.numeric(get(col_rmax))) else NA_real_
)]
tracks[, time_utc := to_posix_safe(time_raw)]
tracks <- tracks[!is.na(time_utc) & !is.na(lat) & !is.na(lon) & !is.na(wind_kt)]
tracks <- tracks[year(time_utc) >= start_year & year(time_utc) <= end_year]
# Convert Vmax to m/s
tracks[, Vmax_ms := wind_kt * KT_TO_MS]
tracks <- tracks[is.finite(Vmax_ms) & Vmax_ms > 0]
if (nrow(tracks) == 0) stop("No usable track points after filtering (check wind/time parsing).")

# ---------------------------
# Build dt_s (seconds between fixes) per storm globally
# Use LEAD difference; last point fallback to 6h

# ---------------------------
setorder(tracks, SID, time_utc)
tracks[, dt_s := as.numeric(difftime(shift(time_utc, type="lead"), time_utc, units="secs")), by = SID]
tracks[is.na(dt_s) | dt_s <= 0 | dt_s > 24*3600, dt_s := 6*3600]

# ---------------------------
# Rmax_m:
#   - use IBTrACS if available and looks valid
#   - else apply Option B intensity bins using wind_kt
# Heuristic for rmax_raw units if provided:
#   - if column name includes "NM" => nautical miles
#   - else if includes "KM" => km
#   - else assume km (common enough, and your Option B is km anyway)

# ---------------------------
rmax_name <- if (!is.na(col_rmax)) toupper(col_rmax) else ""
tracks[, Rmax_m := {
  # default from bins
  rbin_km <- fifelse(wind_kt < 64, Rmax_TS_km,
                    fifelse(wind_kt < 96, Rmax_C12_km, Rmax_C35_km))
  # if we have raw rmax and it's usable, prefer it
  use_raw <- !is.na(rmax_raw) & is.finite(rmax_raw) & rmax_raw > 0
  r_km <- rbin_km
  if (any(use_raw)) {
    if (grepl("NM", rmax_name)) {
      r_km[use_raw] <- rmax_raw[use_raw] * 1.852
    } else if (grepl("KM", rmax_name)) {
      r_km[use_raw] <- rmax_raw[use_raw]
    } else {
      # unknown units -> assume km
      r_km[use_raw] <- rmax_raw[use_raw]
    }
  }
  pmax(r_km, 5) * 1000  # enforce >=5 km
}]
# Holland B (constant)
tracks[, B := B_default]
setkey(tracks, time_utc)

# ---------------------------
# UNIQUE SITES for full-period caching

# ---------------------------
sites <- unique(events[, .(lat, lon)])
sites[, site_id := .I]
full_start <- as.POSIXct(sprintf("%d-01-01 00:00:00", start_year), tz="UTC")
full_end   <- as.POSIXct(sprintf("%d-12-31 23:59:59", end_year), tz="UTC")
message(sprintf("Computing FULL-PERIOD (%d-%d) power within %dkm for %d unique sites using %d cores...",
                start_year, end_year, R_km, nrow(sites), ncores))
site_idx <- split(seq_len(nrow(sites)), cut(seq_len(nrow(sites)), breaks = ncores, labels = FALSE))
full_list <- pblapply(site_idx, function(idx) {
  compute_power_for_window_manysites(
    tracks_dt = tracks,
    sites_dt  = sites[idx, .(lat, lon)],
    win_start = full_start,
    win_end   = full_end,
    R_m       = R_m,
    R_km      = R_km
  )
}, cl = ncores)
sites[, power_1993_2020 := unlist(full_list)]

# ---------------------------
# Merge site_id + full period power back into events

# ---------------------------
events <- merge(events, sites, by = c("lat", "lon"), all.x = TRUE)

# ---------------------------
# EVENT WINDOWS (6m, 12m) in parallel

# ---------------------------
message(sprintf("Computing 6-month and 12-month prior-to-event power within %dkm for %d events using %d cores...",
                R_km, nrow(events), ncores))
events_idx <- split(seq_len(nrow(events)), cut(seq_len(nrow(events)), breaks = ncores, labels = FALSE))
event_results <- pblapply(events_idx, function(idx) {
  sub <- events[idx, .(event_id, lat, lon, date)]
  # Use calendar months (your request) with %m-%
  sub[, `:=`(
    start_6m  = as.POSIXct(date %m-% months(6),  tz="UTC"),
    start_12m = as.POSIXct(date %m-% months(12), tz="UTC"),
    end_t     = as.POSIXct(date, tz="UTC")
  )]
  out6  <- numeric(nrow(sub))
  out12 <- numeric(nrow(sub))
  for (i in seq_len(nrow(sub))) {
    out6[i] <- compute_power_for_window_onesite(
      tracks_dt = tracks,
      lat0 = sub$lat[i], lon0 = sub$lon[i],
      win_start = sub$start_6m[i], win_end = sub$end_t[i],
      R_m = R_m, R_km = R_km
    )
    out12[i] <- compute_power_for_window_onesite(
      tracks_dt = tracks,
      lat0 = sub$lat[i], lon0 = sub$lon[i],
      win_start = sub$start_12m[i], win_end = sub$end_t[i],
      R_m = R_m, R_km = R_km
    )
  }
  data.table(event_id = sub$event_id, power_6m = out6, power_12m = out12)
}, cl = ncores)
event_results_dt <- rbindlist(event_results)
setkey(event_results_dt, event_id)
setkey(events, event_id)
events <- event_results_dt[events]

# ---------------------------
# OUTPUT

# ---------------------------
out <- events[, .(
  Latitude_Degree  = lat,
  Longitude_Degree = lon,
  Date             = date,
  TCpower_6m_400km  = power_6m,
  TCpower_12m_400km = power_12m,
  TCpower_1993_2020_400km = power_1993_2020
)]
write_xlsx(out, out_xlsx)
message("DONE. Saved: ", out_xlsx)
