# Per-animal fix counts removed by each analysis-stage filter (revision, 2026-10-31).
Sys.unsetenv(c("PROJ_LIB", "PROJ_DATA", "GDAL_DATA"))
suppressPackageStartupMessages(library(dplyr))
df <- read.csv("data/processed/all_elephants_screened_v6.csv")
df$timestamp <- as.POSIXct(df$timestamp, format = "%m/%d/%Y %H:%M", tz = "UTC")
df <- df[df$elephant_id != 4557, ]
out <- list()
for (id in sort(unique(df$elephant_id))) {
  s0 <- df[df$elephant_id == id, ]; n0 <- nrow(s0)
  s  <- s0[s0$speed_km_h <= 5, ]; n_speed <- n0 - nrow(s)
  # combined2.R also drops 3905 from 2023-08-08 (dropped collar; overlaps the stationary cut)
  if (id == 3905) s <- s[s$timestamp < as.POSIXct("2023-08-08", tz = "UTC"), ]
  n_date <- n0 - n_speed - nrow(s)
  s <- s[!duplicated(s$timestamp), ]; s <- s[order(s$timestamp), ]
  s <- s[c(diff(as.numeric(s$timestamp)), 1) > 0, ]
  n1 <- nrow(s); thr_used <- NA; maxd <- NA
  for (.i in 1:10) {
    cx <- median(s$x_utm); cy <- median(s$y_utm)
    d <- sqrt((s$x_utm - cx)^2 + (s$y_utm - cy)^2) / 1000
    thr <- max(median(d) + 6 * mad(d), 45); thr_used <- thr
    if (.i == 1) maxd <- max(d)
    ok <- d <= thr; if (all(ok)) break; s <- s[ok, ]
  }
  n_dist <- n1 - nrow(s); n2 <- nrow(s)
  for (.j in 1:10) {
    dt <- as.numeric(diff(s$timestamp), units = "secs") / 3600
    st <- sqrt(diff(s$x_utm)^2 + diff(s$y_utm)^2) / 1000
    sp <- st / pmax(dt, 0.25); si <- c(NA, sp); so <- c(sp, NA)
    bad <- !is.na(si) & !is.na(so) & si > 7 & so > 7
    if (!any(bad)) break; s <- s[!bad, ]
  }
  n_spike <- n2 - nrow(s)
  # steps 4-5 km/h (relevant to TGD Vmax = 4 km/h)
  dt <- as.numeric(diff(s$timestamp), units = "secs") / 3600
  st <- sqrt(diff(s$x_utm)^2 + diff(s$y_utm)^2) / 1000; sp <- st / dt
  out[[length(out) + 1]] <- data.frame(ID = id, screened = n0, speed_gt5 = n_speed, date_3905 = n_date,
    after_speed = n0 - n_speed - n_date,
    dist_filter = n_dist, dist_thr_km = round(thr_used, 1), max_dist_km = round(maxd, 1),
    spike = n_spike, analysed = nrow(s),
    steps_gt4kmh = sum(sp > 4 & dt <= 24, na.rm = TRUE), steps_total = sum(dt <= 24, na.rm = TRUE))
}
r <- do.call(rbind, out); print(r); cat("totals:", colSums(r[, c(2:4, 7:10)]), "\n")
write.csv(r, "final/data_quality/tables/Filter_counts_v6.csv", row.names = FALSE)
