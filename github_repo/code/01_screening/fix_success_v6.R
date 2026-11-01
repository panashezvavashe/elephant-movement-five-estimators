# Fix-success rate per collar from the programmed schedule (Reviewer 1; revision 2026-10-31).
# Expected fixes = sum over consecutive records of round(gap / programmed interval);
# success = fixes logged / expected, over the deployed period up to the last fix
# (or the stationary cut for 3905 and 3907). Accuracy filter NOT applied: a fix
# with poor DOP is still an acquired fix.
Sys.unsetenv(c("PROJ_LIB", "PROJ_DATA", "GDAL_DATA"))
src <- readLines("code/01_screening/build_screened_v6.R", encoding = "UTF-8")
stop_at <- grep("^files <- ", src) - 1
eval(parse(text = src[1:stop_at]))
cuts <- read.csv("final/data_quality/tables/Screening_log_v6.csv")
out <- list()
for (f in sort(list.files(RAW_DIR, "^Elephant [0-9]+[.]csv$", full.names = TRUE))) {
  id <- as.integer(gsub("[^0-9]", "", basename(f)))
  if (id == 4557) next
  r <- read.csv(f, fileEncoding = "latin1", strip.white = TRUE, check.names = FALSE)
  names(r) <- trimws(names(r))
  d <- data.frame(ts = as.POSIXct(r[[4]], tz = "UTC"), iv = r$`Log Interval`,
                  lat = parse_coord(r$Latitude))
  d <- d[!is.na(d$ts) & !is.na(d$lat), ]; d <- d[order(d$ts), ]; d <- d[!duplicated(d$ts), ]
  d <- d[seq(which(d$iv != "Every 6 hours")[1], nrow(d)), ]
  sf_ <- cuts$stationary_from[cuts$ID == id]
  if (!is.na(sf_) && nzchar(sf_)) d <- d[d$ts < as.POSIXct(sf_, tz = "UTC"), ]
  ivh <- as.numeric(gsub("[^0-9]", "", d$iv))
  gap <- diff(as.numeric(d$ts)) / 3600
  slots <- pmax(1, round(gap / ivh[-nrow(d)]))
  exp_n <- sum(slots) + 1
  out[[length(out) + 1]] <- data.frame(ID = id, logged = nrow(d), expected = exp_n,
    success_pct = round(100 * nrow(d) / exp_n, 1),
    gaps_missed = sum(slots > 1), longest_gap_h = round(max(gap), 1))
}
r <- do.call(rbind, out); print(r)
write.csv(r, "final/data_quality/tables/Fix_success_v6.csv", row.names = FALSE)
