# ============================================================
# Phase 2 (revision, 2026-10-28): season-year AKDE + ctmm::meta()
#
# Replaces the pooled wet/dry AKDE (one fit per season across all years) with
# one fit per season-year, so between-year range shifts are not folded into a
# single stationary process. Wet = Nov(Y)-Apr(Y+1) labelled W<Y>; dry = May-Oct(Y)
# labelled D<Y>. Same filters, telemetry object and model selection as
# combined2.R Section D (helpers are taken verbatim from combined2.R).
#
# Outputs -> final/AKDE/seasonal_v6/
#   SeasonYear_AKDE.csv            one row per season-year (model, N_eff, area CI)
#   SeasonYear_excluded.csv        season-years too short to fit, with reason
#   meta_*.txt / meta_*.png        wet vs dry ratio of mean areas (ctmm::meta)
#   PerAnimal_meta.csv             per-animal wet:dry ratio with CI
#   FullVsSeasonal_containment.csv share of each season's fixes outside the
#                                  full-period 95% contour (all animals)
#   Overlay_4142.png               wet vs full-period contours for 4142
#   fits/                          cached fit + UD per season-year (re-runs skip these)
#
#   Rscript rerun_seasonal_v6.R
# ============================================================

Sys.unsetenv(c("PROJ_LIB", "PROJ_DATA", "GDAL_DATA"))  # PostGIS values break sf CRS here
suppressPackageStartupMessages({
  library(ctmm); library(sp); library(sf); library(dplyr)
})

DATA_DIR <- "."   # run from the repository root
setwd(DATA_DIR)
N_CORES  <- 4L          # laptop loses power with more workers
PNG_DPI  <- 300L
FIT_PTS  <- 1500L       # same per-fit cap as combined2.R (MAX_PTS)
MIN_FIX  <- 100L        # minimum fixes in a season-year
MIN_DAYS <- 90          # minimum span (days) of a season-year (season = ~181 d)

AS  <- file.path(DATA_DIR, "final", "AKDE", "spatial")
OUT <- file.path(DATA_DIR, "final", "AKDE", "seasonal_v6")
FIT_DIR <- file.path(OUT, "fits")
dir.create(FIT_DIR, recursive = TRUE, showWarnings = FALSE)

# ── Helper functions from combined2.R ─────────────────────────
src   <- readLines("code/02_estimators/combined2.R", encoding = "UTF-8")
for (e in parse(text = src, keep.source = FALSE)) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) &&
      is.call(e[[3]]) && identical(e[[3]][[1]], as.name("function")))
    eval(e, envir = globalenv())
}
MAX_PTS <- .Machine$integer.max   # prep_sub() thins to MAX_PTS; thin per season-year instead

season_year <- function(ts) {
  m <- as.integer(format(ts, "%m")); y <- as.integer(format(ts, "%Y"))
  ifelse(m >= 11, paste0("W", y), ifelse(m <= 4, paste0("W", y - 1), paste0("D", y)))
}

# ── Data (same as combined2.R Section B) ──────────────────────
df <- read.csv("data/processed/all_elephants_screened_v6.csv")   # build_screened_v6.R
df$timestamp <- as.POSIXct(df$timestamp, format = "%m/%d/%Y %H:%M", tz = "UTC")
df_clean <- df %>%
  filter(speed_km_h <= 5, !elephant_id %in% c(4557L)) %>%
  arrange(elephant_id, timestamp)
# 3905's collar is stationary (weekly SD ~40 m) from 2023-08-08 to the end of
# its record — dropped collar; those fixes are not movement.
df_clean <- df_clean %>%
  filter(!(elephant_id == 3905 & timestamp >= as.POSIXct("2023-08-08", tz = "UTC")))
IDS <- sort(unique(df_clean$elephant_id))

# ── 1. Fit every season-year ──────────────────────────────────
rows <- list(); excl <- list(); uds <- list()
for (id in IDS) {
  sub <- prep_sub(df_clean, id)
  sub$sy <- season_year(sub$timestamp)
  for (sy in unique(sub$sy)) {
    s <- sub[sub$sy == sy, ]
    span <- as.numeric(difftime(max(s$timestamp), min(s$timestamp), units = "days"))
    if (nrow(s) < MIN_FIX || span < MIN_DAYS) {
      excl[[length(excl) + 1]] <- data.frame(ID = id, season_year = sy, n_fixes = nrow(s),
        span_days = round(span, 1),
        reason = if (nrow(s) < MIN_FIX) "too few fixes" else "span too short")
      next
    }
    s <- thin_to(s, FIT_PTS)
    cache <- file.path(FIT_DIR, sprintf("AKDE_sy_%s_%s.rds", id, sy))
    res <- if (file.exists(cache)) readRDS(cache) else tryCatch({
      cat(sprintf("[%s] fitting %s %s (%d fixes, %.0f d)\n",
                  format(Sys.time(), "%H:%M"), id, sy, nrow(s), span))
      tel  <- make_tel(s)
      fits <- ctmm.select(tel, ctmm.guess(tel, interactive = FALSE),
                          cores = N_CORES, verbose = TRUE)
      tab  <- summary(fits)
      ud   <- akde(tel, fits[[1]])
      r <- list(fit = fits[[1]], ud = ud, model = names(fits)[1],
                runner_up = if (length(fits) > 1) names(fits)[2] else NA_character_,
                dAICc_2nd = if (nrow(tab) > 1) round(tab[2, grep("AIC", colnames(tab))[1]], 2) else NA_real_,
                n_fixes = nrow(s), span_days = span,
                start = min(s$timestamp), end = max(s$timestamp))
      saveRDS(r, cache); r
    }, error = function(e) { cat("  FAILED:", conditionMessage(e), "\n"); NULL })
    if (is.null(res)) {
      excl[[length(excl) + 1]] <- data.frame(ID = id, season_year = sy, n_fixes = nrow(s),
        span_days = round(span, 1), reason = "fit failed")
      next
    }
    ci  <- summary(res$ud, units = FALSE)$CI[1, ] / 1e6
    dof <- summary(res$ud, units = FALSE)$DOF["area"]
    rows[[length(rows) + 1]] <- data.frame(
      ID = id, season = ifelse(substr(sy, 1, 1) == "W", "wet", "dry"), season_year = sy,
      start = format(res$start, "%Y-%m-%d"), end = format(res$end, "%Y-%m-%d"),
      n_fixes = res$n_fixes, span_days = round(res$span_days, 1),
      model = res$model, runner_up = res$runner_up, dAICc_2nd = res$dAICc_2nd,
      N_eff = round(unname(dof), 1),
      area_low = ci[1], area_est = ci[2], area_high = ci[3])
    uds[[paste(id, sy, sep = "_")]] <- res$ud
  }
}
sy_tab <- do.call(rbind, rows)
write.csv(sy_tab, file.path(OUT, "SeasonYear_AKDE.csv"), row.names = FALSE)
if (length(excl)) write.csv(do.call(rbind, excl), file.path(OUT, "SeasonYear_excluded.csv"),
                            row.names = FALSE)
cat(sprintf("\n%d season-year fits; %d excluded\n", nrow(sy_tab), length(excl)))

# ── 2. Population wet vs dry: ctmm::meta() ────────────────────
# meta() on a named list of two UD lists returns a [2, 2, low/est/high] array of
# ratios of population mean areas; ["wet/", "/dry", ] is wet:dry.
meta_summary <- function(keys_wet, keys_dry, tag = NA, plot_file = NULL) {
  x <- list(wet = uds[keys_wet], dry = uds[keys_dry])
  grp <- lapply(x, function(u) if (length(u) > 1)
    tryCatch(meta(u, variable = "area", plot = FALSE), error = function(e) NULL) else NULL)
  if (!is.null(plot_file)) png(plot_file, width = 8, height = 7, units = "in", res = PNG_DPI)
  txt <- capture.output(
    r <- tryCatch(meta(x, variable = "area", level = 0.95, plot = !is.null(plot_file)),
                  error = function(e) { cat("meta failed:", conditionMessage(e), "\n"); NULL }))
  if (!is.null(plot_file)) dev.off()
  mean_row <- function(g) if (is.null(g)) rep(NA_real_, 3) else unname(g[grep("^mean", rownames(g))[1], ])
  data.frame(tag = tag, n_wet = length(keys_wet), n_dry = length(keys_dry),
    wet_mean = mean_row(grp$wet)[2], wet_low = mean_row(grp$wet)[1], wet_high = mean_row(grp$wet)[3],
    dry_mean = mean_row(grp$dry)[2], dry_low = mean_row(grp$dry)[1], dry_high = mean_row(grp$dry)[3],
    ratio_wet_dry = if (is.null(r)) NA else r["wet/", "/dry", "est"],
    ratio_low     = if (is.null(r)) NA else r["wet/", "/dry", "low"],
    ratio_high    = if (is.null(r)) NA else r["wet/", "/dry", "high"],
    log = paste(txt, collapse = " | "))
}
run_meta <- function(keys_wet, keys_dry, tag) {
  s <- meta_summary(keys_wet, keys_dry, tag, file.path(OUT, paste0("meta_", tag, ".png")))
  writeLines(c(sprintf("=== %s ===", tag), strsplit(s$log, " \\| ")[[1]]),
             file.path(OUT, paste0("meta_", tag, ".txt")))
  print(s[, names(s) != "log"])
  s
}
key <- paste(sy_tab$ID, sy_tab$season_year, sep = "_")
w <- sy_tab$season == "wet";  no4501 <- sy_tab$ID != 4501

# (a) all season-years, as the reviewers asked
pop <- list()
pop[[1]] <- run_meta(key[w], key[!w], "all_season_years")
pop[[2]] <- run_meta(key[w & no4501], key[!w & no4501], "all_season_years_no4501")

# (b) one season-year per animal per season (highest N_eff): avoids letting
#     animals with more years dominate, and treats animals as the replicate
best <- sy_tab %>% mutate(key = key) %>% group_by(ID, season) %>%
  slice_max(N_eff, n = 1, with_ties = FALSE) %>% ungroup()
pop[[3]] <- run_meta(best$key[best$season == "wet"], best$key[best$season == "dry"], "one_per_animal")
pop[[4]] <- run_meta(best$key[best$season == "wet" & best$ID != 4501],
                     best$key[best$season == "dry" & best$ID != 4501], "one_per_animal_no4501")
write.csv(do.call(rbind, pop)[, -ncol(pop[[1]])], file.path(OUT, "Meta_population.csv"), row.names = FALSE)

# ── 3. Per-animal wet:dry ratio across that animal's season-years ──
per_animal <- dplyr::bind_rows(lapply(IDS, function(id) {
  kw <- key[w & sy_tab$ID == id];  kd <- key[!w & sy_tab$ID == id]
  if (!length(kw) || !length(kd)) return(NULL)
  s <- meta_summary(kw, kd, tag = id)
  s$ID <- id
  s$CI_excludes_1 <- !is.na(s$ratio_low) & (s$ratio_low > 1 | s$ratio_high < 1)
  s[, c("ID", setdiff(names(s), c("ID", "tag", "log")))]
}))
write.csv(per_animal, file.path(OUT, "PerAnimal_meta.csv"), row.names = FALSE)
print(per_animal)

# ── 4. Seasonal vs full-period AKDE: containment diagnostic ───
contour95 <- function(ud) {
  p <- sf::st_as_sf(SpatialPolygonsDataFrame.UD(ud, level.UD = 0.95, level = 0.95))
  p[grepl("est", p$name), ]
}
contain <- dplyr::bind_rows(lapply(IDS, function(id) {
  f_full <- file.path(AS, paste0("AKDE_", id, ".rds"))
  if (!file.exists(f_full)) return(NULL)
  ud_full <- readRDS(f_full)
  poly <- contour95(ud_full)
  sub  <- prep_sub(df_clean, id)
  pts  <- sf::st_transform(sf::st_as_sf(sub, coords = c("longitude", "latitude"), crs = 4326),
                           sf::st_crs(poly))
  inside <- lengths(sf::st_intersects(pts, poly)) > 0
  wet <- is_wet(sub$timestamp)
  out <- data.frame(ID = id,
    full_95_km2 = as.numeric(sf::st_area(poly)) / 1e6,
    n_wet = sum(wet), pct_wet_outside_full95 = round(100 * mean(!inside[wet]), 1),
    n_dry = sum(!wet), pct_dry_outside_full95 = round(100 * mean(!inside[!wet]), 1),
    pct_all_outside_full95 = round(100 * mean(!inside), 1))
  if (id == 4142) {
    f_wet <- file.path(AS, "AKDE_wet_4142.rds")
    if (file.exists(f_wet)) {
      poly_w <- sf::st_transform(contour95(readRDS(f_wet)), sf::st_crs(poly))
      png(file.path(OUT, "Overlay_4142.png"), width = 8, height = 7, units = "in", res = PNG_DPI)
      plot(sf::st_geometry(poly_w), border = "#2166AC", lwd = 2,
           main = "Elephant 4142: wet-season vs full-period AKDE 95%")
      plot(sf::st_geometry(poly), border = "black", lwd = 2, add = TRUE)
      plot(sf::st_geometry(pts[wet, ]), pch = 20, cex = 0.3, col = "#2166AC80", add = TRUE)
      plot(sf::st_geometry(pts[!wet, ]), pch = 20, cex = 0.3, col = "#B2182B80", add = TRUE)
      legend("topright", bty = "n",
             legend = c("Full-period 95%", "Pooled wet 95%", "Wet fixes", "Dry fixes"),
             col = c("black", "#2166AC", "#2166AC", "#B2182B"),
             lty = c(1, 1, NA, NA), pch = c(NA, NA, 20, 20), lwd = 2)
      dev.off()
      out$wet_pooled_95_km2 <- as.numeric(sf::st_area(poly_w)) / 1e6
    }
  }
  out
}))
write.csv(contain, file.path(OUT, "FullVsSeasonal_containment.csv"), row.names = FALSE)
print(contain)
cat("\nDone. Outputs in", OUT, "\n")
