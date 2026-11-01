# ============================================================
# Build the screened GPS data set from the raw AWT collar files (revision, 2026-10-30)
#
# Replaces all_elephants_combined.csv, whose upstream rule kept a fix only when it
# came 2-4 h after the previous raw fix. That rule discarded every hourly fix
# (3906, 3907, 4501), every genuine 6-hourly period (3906/3907 2023, 4142 2021)
# and the first fix after any gap > 4 h.
#
# Rule applied here, per animal, in order (counts logged to Screening_log_v6.csv):
#   1. parse; drop records with unparseable time or coordinates
#   2. drop duplicate timestamps (keep the lower collar-reported accuracy value)
#   3. pre-deployment: drop the initial 6-hourly run (collar in transit; spread
#      45-180 km); deployment starts at the first fix on another schedule
#   4. accuracy: drop fixes with collar-reported accuracy ("DOP" field, m) > 30 m
#   5. stationary collar: a stretch >= 14 days in which the 7-day rolling spread of
#      fixes is < 100 m marks a dropped collar; the track is truncated at its start
#      (less half the window), which also removes post-recovery transit fixes
#   6. regularise to a 3-hour minimum spacing: keep a fix only if >= 2.5 h after
#      the last kept fix (hourly schedules -> ~3-hourly; 3- and 6-hourly unchanged)
# Speed spikes, distance outliers and thinning to <= 1,500 fixes are applied later,
# per estimator, in combined2.R (Section B filter + prep_sub()).
#
# Output: data/processed/all_elephants_screened_v6.csv (same columns as all_elephants_combined.csv)
#         final/data_quality/tables/Screening_log_v6.csv
#
#   Rscript build_screened_v6.R
# ============================================================

Sys.unsetenv(c("PROJ_LIB", "PROJ_DATA", "GDAL_DATA"))
suppressPackageStartupMessages({ library(sf); library(dplyr) })

RAW_DIR  <- "data/raw"   # one AWT export per collar: "Elephant <ID>.csv"
DATA_DIR <- "."          # run from the repository root
setwd(DATA_DIR)

MAX_DOP_M     <- 30
STAT_WIN_D    <- 7      # rolling window for the stationary test (days)
STAT_SD_M     <- 100    # spread below which the collar is treated as stationary
STAT_MIN_D    <- 14     # minimum length of a stationary stretch (days)
MIN_SPACING_H <- 2.5

parse_coord <- function(s) {
  v <- as.numeric(sub("^[[:space:]]*([0-9.]+).*$", "\\1", s))
  ifelse(grepl("[SW]", s), -v, v)
}

# TRUE for fixes inside a stationary stretch
flag_stationary <- function(d) {
  n <- nrow(d); t <- as.numeric(d$ts) / 86400
  sd_win <- vapply(seq_len(n), function(i) {
    j <- which(t >= t[i] - STAT_WIN_D / 2 & t <= t[i] + STAT_WIN_D / 2)
    if (length(j) < 10) return(NA_real_)
    sqrt(var(d$x[j]) + var(d$y[j]))
  }, numeric(1))
  low <- !is.na(sd_win) & sd_win < STAT_SD_M
  r <- rle(low); ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1
  out <- rep(FALSE, n)
  for (k in which(r$values)) {
    if (t[ends[k]] - t[starts[k]] >= STAT_MIN_D) out[starts[k]:ends[k]] <- TRUE
  }
  out
}

thin_spacing <- function(ts, min_h) {
  keep <- logical(length(ts)); last <- -Inf; tt <- as.numeric(ts) / 3600
  for (i in seq_along(tt)) if (tt[i] - last >= min_h) { keep[i] <- TRUE; last <- tt[i] }
  keep
}

files <- sort(list.files(RAW_DIR, "^Elephant [0-9]+[.]csv$", full.names = TRUE))
log_rows <- list(); out <- list()
for (f in files) {
  id <- as.integer(sub("^Elephant ([0-9]+)[.]csv$", "\\1", basename(f)))
  r <- read.csv(f, fileEncoding = "latin1", strip.white = TRUE, check.names = FALSE)
  names(r) <- trimws(names(r))
  d <- data.frame(ts  = as.POSIXct(r[[4]], tz = "UTC"),
                  lat = parse_coord(r$Latitude), lon = parse_coord(r$Longitude),
                  dop = as.numeric(sub("^([0-9.]+).*$", "\\1", r$DOP)),
                  iv  = r$`Log Interval`)
  cnt <- c(raw = nrow(d))
  d <- d[!is.na(d$ts) & !is.na(d$lat) & !is.na(d$lon), ];       cnt["parsed"] <- nrow(d)
  d <- d[order(d$ts, d$dop), ]; d <- d[!duplicated(d$ts), ];     cnt["no_duplicates"] <- nrow(d)
  first_dep <- which(d$iv != "Every 6 hours")[1]
  deploy_start <- d$ts[first_dep]
  d <- d[seq(first_dep, nrow(d)), ];                             cnt["deployed"] <- nrow(d)
  d <- d[is.na(d$dop) | d$dop <= MAX_DOP_M, ];                   cnt["accuracy_le_30m"] <- nrow(d)
  p <- st_coordinates(st_transform(st_as_sf(d, coords = c("lon", "lat"), crs = 4326), 32736))
  d$x <- p[, 1]; d$y <- p[, 2]
  # A dropped collar does not resume: truncate at the first stationary stretch,
  # moved back by half the rolling window (the centred window flags it late).
  st <- flag_stationary(d)
  stat_from <- NA
  if (any(st)) {
    cut_ts <- min(d$ts[st]) - STAT_WIN_D / 2 * 86400
    stat_from <- format(cut_ts)
    d <- d[d$ts < cut_ts, ]
  }
  cnt["before_stationary"] <- nrow(d)
  d <- d[thin_spacing(d$ts, MIN_SPACING_H), ];                   cnt["spacing_ge_2.5h"] <- nrow(d)

  d$gap_hours  <- c(NA, diff(as.numeric(d$ts)) / 3600)
  step_km      <- c(NA, sqrt(diff(d$x)^2 + diff(d$y)^2) / 1000)
  d$speed_km_h <- round(ifelse(is.na(d$gap_hours), 0, step_km / d$gap_hours), 2)
  out[[as.character(id)]] <- data.frame(
    elephant_id = id, timestamp = format(d$ts, "%m/%d/%Y %H:%M"),
    latitude = d$lat, longitude = d$lon, gps_error_m = d$dop,
    speed_km_h = d$speed_km_h, speed_m_hr = round(d$speed_km_h * 1000),
    x_utm = d$x, y_utm = d$y, gap_hours = d$gap_hours)
  log_rows[[as.character(id)]] <- data.frame(ID = id, t(cnt), deploy_start = format(deploy_start),
    first_fix = format(min(d$ts)), last_fix = format(max(d$ts)), stationary_from = stat_from,
    schedules = paste(names(table(d$iv)), table(d$iv), sep = ": ", collapse = "; "))
}
all <- do.call(rbind, out)
write.csv(all, "data/processed/all_elephants_screened_v6.csv", row.names = FALSE)
lg <- do.call(rbind, log_rows)
dir.create(file.path(DATA_DIR, "final", "data_quality", "tables"), recursive = TRUE, showWarnings = FALSE)
write.csv(lg, file.path(DATA_DIR, "final", "data_quality", "tables", "Screening_log_v6.csv"), row.names = FALSE)
print(lg, row.names = FALSE)
cat("\nTotal screened fixes (incl. 4557):", nrow(all), "\n")
