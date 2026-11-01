# ============================================================
# ELEPHANT MOVEMENT ANALYSIS — MASTERS LEVEL CORRECT AND Final script compiled 
# Mana Pools, Zimbabwe Panashe ZVAVASHE  MPHIL @2026
# Tabs included: Variogram & Model, Data Quality, Seasonal Tests,
#                Method Comparison, Space Overlap, Movement Behaviour

# =============================================================================
# ELEPHANT GPS TELEMETRY — COMPLETE UNIFIED PIPELINE
# Mana Pools NP, Zimbabwe
#
# WHAT THIS SCRIPT DOES (in order):
#   PART 1 — Installs all packages
#   PART 2 — Reads every CSV in your data folder, runs full analysis per elephant:
#             CTMM/AKDE, Seasonal AKDE, dBBMM, HMM, iSSF, k/r/a-LoCoH,
#             Classical KDE, MCP, RF/SVM corridor models
#   PART 3 — Exports all publication-quality figures and tables (300 dpi PNG)
#             for Objectives 1 & 2 into Results/Publication_Figures/
#
# HOW TO RUN — one line in your R console:
#   (run from the repository root: Rscript code/02_estimators/combined2.R)
#
# Or in RStudio: open this file → click the SOURCE button (top right of editor)
# DO NOT use Ctrl+Enter to run line by line — run the whole file at once.
# ============================================================

# ── 0. Packages ───────────────────────────────────────────
suppressPackageStartupMessages({
  library(ctmm)          # AKDE, variogram, model selection, overlap
  library(move)          # dBBMM
  library(adehabitatHR)  # aLOCOH, kerneloverlaphr
  library(adehabitatLT)  # first-passage time (FPT)
  library(sf)
  library(ggplot2)
  library(ggrepel)       # geom_text_repel for Bland-Altman labels
  library(dplyr)
  library(tidyr)
  library(raster)
  library(parallel)
  library(viridis)
  library(circular)      # turning-angle statistics
  library(cowplot)       # multi-panel figure assembly
  library(rstatix)       # wilcox_effsize
  library(reshape2)      # melt() for heatmap
  library(vegan)         # mantel test
})

# ── Data directory ────────────────────────────────────────
DATA_DIR <- "."   # run from the repository root
setwd(DATA_DIR)

# ── Global settings ───────────────────────────────────────
N_CORES  <- min(4L, max(1L, detectCores() - 1L))
MAX_PTS  <- 1500L
PNG_DPI  <- 300L
ETD_SIGMA <- 500L   # ETD kernel sigma (m), fixed for all animals — see Section G2
cat(sprintf("Cores: %d  |  Max pts/elephant: %d  |  DPI: %d\n",
            N_CORES, MAX_PTS, PNG_DPI))

utm_proj <- CRS("+proj=utm +zone=36 +south +datum=WGS84 +units=m +no_defs")

# ── Output directory structure ────────────────────────────
FINAL_DIR <- file.path(DATA_DIR, "final")

AG <- file.path(FINAL_DIR, "AKDE",         "graphs")
AT <- file.path(FINAL_DIR, "AKDE",         "tables")
AS <- file.path(FINAL_DIR, "AKDE",         "spatial")

DG <- file.path(FINAL_DIR, "dBBMM",        "graphs")
DS <- file.path(FINAL_DIR, "dBBMM",        "spatial")

LG <- file.path(FINAL_DIR, "aLOCOH",       "graphs")
LS <- file.path(FINAL_DIR, "aLOCOH",       "spatial")

TG <- file.path(FINAL_DIR, "TGD",          "graphs")
TS <- file.path(FINAL_DIR, "TGD",          "spatial")

EG <- file.path(FINAL_DIR, "ETD",          "graphs")
ET <- file.path(FINAL_DIR, "ETD",          "tables")
ES <- file.path(FINAL_DIR, "ETD",          "spatial")

CG <- file.path(FINAL_DIR, "comparison",   "graphs")
CT <- file.path(FINAL_DIR, "comparison",   "tables")

QG <- file.path(FINAL_DIR, "data_quality", "graphs")
QT <- file.path(FINAL_DIR, "data_quality", "tables")

MG <- file.path(FINAL_DIR, "movement",     "graphs")
MT <- file.path(FINAL_DIR, "movement",     "tables")

for (d in c(AG,AT,AS, DG,DS, LG,LS, TG,TS, EG,ET,ES, CG,CT, QG,QT, MG,MT))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("Output root: %s\n", FINAL_DIR))

# ── Migrate pre-existing results from cleaned/ root into final/ ──────────────
# Copies files produced by older script runs (flat layout) into the new
# method-specific subfolders.  Already-present files are not overwritten.
.migrate <- function(src, dst) {
  if (file.exists(src) && !file.exists(dst))
    tryCatch(file.copy(src, dst), error = function(e) NULL)
}
.copy_shp <- function(stem_src, stem_dst) {
  for (ext in c(".shp", ".dbf", ".prj", ".shx"))
    .migrate(paste0(stem_src, ext), paste0(stem_dst, ext))
}
IDS_ALL <- c(3905, 3906, 3907, 4141, 4142, 4501, 4556, 4558, 4559, 4560)
OLD     <- DATA_DIR
for (.id in IDS_ALL) {
  # AKDE
  .migrate(file.path(OLD, paste0("AKDE_",          .id, ".tif")),          file.path(AS, paste0("AKDE_",          .id, ".tif")))
  .migrate(file.path(OLD, paste0("AKDE_",          .id, "_full.png")),     file.path(AG, paste0("AKDE_",          .id, "_full.png")))
  .migrate(file.path(OLD, paste0("AKDE_",          .id, "_seasonal.png")), file.path(AG, paste0("AKDE_",          .id, "_seasonal.png")))
  .migrate(file.path(OLD, paste0("Variogram_",     .id, ".png")),          file.path(AG, paste0("Variogram_",     .id, ".png")))
  .migrate(file.path(OLD, paste0("ModelSelection_",.id, ".csv")),          file.path(AT, paste0("ModelSelection_",.id, ".csv")))
  # dBBMM
  .migrate(file.path(OLD, paste0("dBBMM_",  .id, ".tif")),          file.path(DS, paste0("dBBMM_",  .id, ".tif")))
  .migrate(file.path(OLD, paste0("dBBMM_",  .id, "_ud.png")),       file.path(DG, paste0("dBBMM_",  .id, "_ud.png")))
  .migrate(file.path(OLD, paste0("dBBMM_",  .id, "_isopleths.png")),file.path(DG, paste0("dBBMM_",  .id, "_isopleths.png")))
  for (.pct in c("50","80","90","95","99"))
    .copy_shp(file.path(OLD, paste0("dBBMM_", .pct, "_", .id)),
              file.path(DS,  paste0("dBBMM_", .pct, "_", .id)))
  # aLOCOH
  .migrate(file.path(OLD, paste0("aLOCOH_",    .id, ".png")), file.path(LG, paste0("aLOCOH_",    .id, ".png")))
  .copy_shp(file.path(OLD, paste0("aLOCOH_50_",.id)),         file.path(LS, paste0("aLOCOH_50_",.id)))
  .copy_shp(file.path(OLD, paste0("aLOCOH_95_",.id)),         file.path(LS, paste0("aLOCOH_95_",.id)))
  # TGD
  .migrate(file.path(OLD, paste0("TGD_",.id, ".png")), file.path(TG, paste0("TGD_",.id, ".png")))
  .copy_shp(file.path(OLD, paste0("TGD_",.id)),         file.path(TS, paste0("TGD_",.id)))
  # ETD — graphs
  for (.sfx in c("_full.png","_speed_stratified.png")) {
    fp <- file.path(OLD, paste0("ETD_", .id, .sfx))
    if (file.exists(fp)) .migrate(fp, file.path(EG, paste0("ETD_", .id, .sfx)))
  }
  .migrate(file.path(OLD, paste0("ETD_sigma_sensitivity_",.id,".png")),
           file.path(EG, paste0("ETD_sigma_sensitivity_",.id,".png")))
  # ETD — tables
  .migrate(file.path(OLD, paste0("ETD_sigma_sensitivity_",.id,".csv")),
           file.path(ET, paste0("ETD_sigma_sensitivity_",.id,".csv")))
  # ETD — spatial
  for (.sfx in c("_full.tif","_slow.tif","_fast.tif"))
    .migrate(file.path(OLD, paste0("ETD_", .id, .sfx)), file.path(ES, paste0("ETD_", .id, .sfx)))
  .copy_shp(file.path(OLD, paste0("ETD_50_",.id)), file.path(ES, paste0("ETD_50_",.id)))
  .copy_shp(file.path(OLD, paste0("ETD_95_",.id)), file.path(ES, paste0("ETD_95_",.id)))
  # data_quality
  .migrate(file.path(OLD, paste0("FixTimeline_",.id,".png")), file.path(QG, paste0("FixTimeline_",.id,".png")))
  # movement
  .migrate(file.path(OLD, paste0("MovementBehaviour_",.id,".png")), file.path(MG, paste0("MovementBehaviour_",.id,".png")))
}
# Shared outputs
.migrate(file.path(OLD, "AKDE_CI_Summary.png"),    file.path(AG, "AKDE_CI_Summary.png"))
.migrate(file.path(OLD, "AKDE_Parameters.csv"),    file.path(AT, "AKDE_Parameters.csv"))
.migrate(file.path(OLD, "Data_Quality_Summary.csv"),file.path(QT, "Data_Quality_Summary.csv"))
.migrate(file.path(OLD, "ETD_Comparison.csv"),     file.path(ET, "ETD_Comparison.csv"))
.migrate(file.path(OLD, "Diurnal_Pattern.png"),    file.path(MG, "Diurnal_Pattern.png"))
.migrate(file.path(OLD, "StepLength_Summary.csv"), file.path(MT, "StepLength_Summary.csv"))
cat("Migration of existing outputs into final/ complete\n")

# ============================================================
# SECTION A — HELPERS
# ============================================================

thin_to <- function(df, n = MAX_PTS) {
  if (nrow(df) <= n) return(df)
  df[round(seq(1, nrow(df), length.out = n)), ]
}

# max_pts: AKDE and aLoCoH use the thinned track (MAX_PTS); the path-based estimators
# (dBBMM, TGD, ETD) and movement metrics use the full 3-hourly track (max_pts = Inf).
prep_sub <- function(df_clean, id, max_step_km = 15, max_pts = MAX_PTS) {
  sub <- df_clean[df_clean$elephant_id == id, ]
  sub <- sub[!duplicated(sub$timestamp), ]
  sub <- sub[order(sub$timestamp), ]
  diffs <- c(diff(as.numeric(sub$timestamp)), 1)
  sub   <- sub[diffs > 0, ]
  n_raw <- nrow(sub)

  if (nrow(sub) > 10) {
    for (.i in seq_len(10)) {
      cx  <- median(sub$x_utm);  cy <- median(sub$y_utm)
      d   <- sqrt((sub$x_utm - cx)^2 + (sub$y_utm - cy)^2) / 1000
      thr <- max(median(d) + 6 * mad(d), max_step_km * 3)
      ok  <- d <= thr
      if (all(ok)) break
      sub <- sub[ok, ]
    }
  }

  if (nrow(sub) > 2) {
    for (.j in seq_len(10)) {
      dt      <- as.numeric(diff(sub$timestamp), units = "secs") / 3600
      step_km <- sqrt(diff(sub$x_utm)^2 + diff(sub$y_utm)^2) / 1000
      spd     <- step_km / pmax(dt, 0.25)
      spd_in  <- c(NA_real_, spd)
      spd_out <- c(spd, NA_real_)
      bad     <- !is.na(spd_in) & !is.na(spd_out) & spd_in > 7 & spd_out > 7
      if (!any(bad)) break
      sub <- sub[!bad, ]
    }
  }

  if (nrow(sub) < n_raw)
    cat(sprintf("  Elephant %s: filters removed %d of %d points\n",
                id, n_raw - nrow(sub), n_raw))
  thin_to(sub, max_pts)
}

save_png <- function(expr, path, width = 8, height = 7) {
  png(path, width = width, height = height, units = "in", res = PNG_DPI)
  tryCatch(expr, finally = dev.off())
}

# gap_hours = 24 matches the TGD cut-off; 6 would split every 6-hourly schedule
# into single-fix bursts (revision, 2026-10-30)
prep_bursts <- function(sub, min_burst_size, gap_hours = 24) {
  sub      <- sub[order(sub$timestamp), ]
  dt_h     <- c(Inf, as.numeric(diff(sub$timestamp), units = "secs") / 3600)
  burst_id <- cumsum(dt_h > gap_hours)
  keep     <- as.integer(names(which(table(burst_id) >= min_burst_size)))
  sub[burst_id %in% keep, ]
}

is_wet <- function(ts) { m <- as.integer(format(ts, "%m")); m >= 11 | m <= 4 }

make_tel <- function(sub) {
  as.telemetry(data.frame(
    timestamp = sub$timestamp,
    longitude = as.double(sub$longitude),
    latitude  = as.double(sub$latitude),
    HDOP      = as.double(sub$gps_error_m)
  ), timezone = "UTC")
}

# Helper: build estUD from a DBBMM object and extract an sf isopleth polygon
.dbbmm_isopleth <- function(dbbmm_obj, pct) {
  tryCatch({
    px  <- as(dbbmm_obj, "SpatialPixelsDataFrame")
    ud  <- new("estUD", px)
    ud@vol        <- FALSE
    ud@h$meth     <- "dBBMM"
    poly <- adehabitatHR::getverticeshr(ud, percent = pct, standardize = TRUE)
    out  <- sf::st_as_sf(poly)
    sf::st_crs(out) <- 32736
    out
  }, error = function(e) NULL)
}

# Helper: padded bounding extent of all cells at/under a cumulative-volume
# threshold (e.g. 0.99 for the 99% isopleth), used to zoom dBBMM plots to the
# relevant area instead of the full (much larger) raster extent.
.pct_extent <- function(r, threshold, pad_frac = 0.10) {
  vals <- raster::values(r)
  m    <- vals <= threshold
  if (!any(m, na.rm = TRUE)) return(raster::extent(r))
  pts  <- raster::xyFromCell(r, which(m))
  xr   <- range(pts[, 1]); yr <- range(pts[, 2])
  padx <- max(diff(xr) * pad_frac, raster::res(r)[1] * 2)
  pady <- max(diff(yr) * pad_frac, raster::res(r)[2] * 2)
  e    <- raster::extent(xr[1] - padx, xr[2] + padx, yr[1] - pady, yr[2] + pady)
  raster::intersect(e, raster::extent(r))
}

ud_poly <- function(r, pct) {
  thr  <- pct / 100
  vals <- raster::values(r)
  if (!any(vals <= thr, na.rm = TRUE)) return(NULL)
  r2   <- r
  raster::values(r2) <- ifelse(vals <= thr, 1L, NA_integer_)
  fact <- min(5L, max(1L, as.integer(max(nrow(r2), ncol(r2)) / 400)))
  if (fact > 1L) r2 <- raster::aggregate(r2, fact = fact, fun = max, na.rm = TRUE)
  out <- sf::st_as_sf(raster::rasterToPolygons(r2, dissolve = TRUE))
  sf::st_crs(out) <- 32736
  out
}

# ── ETD helpers ───────────────────────────────────────────────
compute_etd <- function(sub, sigma, res = 250) {
  sub <- sub[order(sub$timestamp), ]
  n   <- nrow(sub)
  if (n < 3) return(NULL)
  dt_hr  <- as.numeric(diff(sub$timestamp), units = "secs") / 3600
  dx     <- diff(sub$x_utm);  dy <- diff(sub$y_utm)
  dist_m <- sqrt(dx^2 + dy^2)
  bear   <- atan2(dy, dx)
  mx     <- (sub$x_utm[-1] + sub$x_utm[-n]) / 2
  my     <- (sub$y_utm[-1] + sub$y_utm[-n]) / 2
  a      <- dist_m / 2 + sigma
  b      <- rep(sigma, n - 1)
  xmin   <- min(sub$x_utm) - 2 * sigma;  xmax <- max(sub$x_utm) + 2 * sigma
  ymin   <- min(sub$y_utm) - 2 * sigma;  ymax <- max(sub$y_utm) + 2 * sigma
  r <- raster::raster(xmn = xmin, xmx = xmax, ymn = ymin, ymx = ymax,
                      res = res,
                      crs = CRS("+proj=utm +zone=36 +south +datum=WGS84 +units=m"))
  vals       <- numeric(raster::ncell(r))
  cells_xy   <- raster::coordinates(r)
  total_time <- sum(dt_hr, na.rm = TRUE)
  for (i in seq_len(n - 1)) {
    if (is.na(dt_hr[i]) || dt_hr[i] <= 0) next
    cx <- cells_xy[, 1] - mx[i];  cy <- cells_xy[, 2] - my[i]
    cos_b <- cos(bear[i]);  sin_b <- sin(bear[i])
    u <-  cos_b * cx + sin_b * cy
    v <- -sin_b * cx + cos_b * cy
    kern <- exp(-0.5 * ((u / a[i])^2 + (v / b[i])^2)) / (2 * pi * a[i] * b[i])
    vals <- vals + dt_hr[i] * kern
  }
  vals <- vals / total_time
  vals <- vals / (sum(vals) * res^2 / 1e6)
  raster::values(r) <- vals
  r
}

etd_volume_area <- function(ud, level = 0.95, res = 250) {
  # 95% volume-isopleth area (km²) from cell counts — same cumulative-density
  # threshold as etd_isopleth(), without the polygonisation.
  v <- raster::values(ud);  v[is.na(v)] <- 0
  v_sorted <- sort(v, decreasing = TRUE)
  k <- which(cumsum(v_sorted) * (res / 1000)^2 >= level)[1]
  if (is.na(k)) return(NA_real_)
  k * (res / 1000)^2
}

sigma_sensitivity <- function(sub, sigmas = seq(250, 2000, 250), res = 250) {
  # Fixed 2026-10-28: previously counted every cell above the 5th percentile of
  # positive density (≈ the whole kernel footprint), not the 95% volume isopleth
  # used for the reported ETD areas, so the curve and the reported areas disagreed.
  areas <- sapply(sigmas, function(s) {
    ud <- compute_etd(sub, sigma = s, res = res)
    if (is.null(ud)) return(NA_real_)
    etd_volume_area(ud, level = 0.95, res = res)
  })
  data.frame(sigma_m = sigmas, area_km2 = areas)
}

etd_isopleth <- function(ud, level = 0.95, res = 250) {
  v <- raster::values(ud);  v[is.na(v)] <- 0
  v_sorted  <- sort(v, decreasing = TRUE)
  cum_sum   <- cumsum(v_sorted) * (res / 1000)^2
  threshold <- v_sorted[which(cum_sum >= level)[1]]
  if (is.na(threshold)) return(NULL)
  mask <- ud >= threshold
  raster::values(mask) <- ifelse(raster::values(mask) > 0, 1L, NA_integer_)
  out <- sf::st_as_sf(raster::rasterToPolygons(mask, dissolve = TRUE))
  out <- dplyr::mutate(out, level = level)
  sf::st_crs(out) <- 32736
  out
}

# Direct first-passage-time variance calculation (bypasses adehabitatLT quirks).
# Returns a numeric vector of var(log(FPT)) per radius (NA where insufficient data).
.fpt_varlog <- function(sub, radii) {
  x  <- sub$x_utm;  y <- sub$y_utm
  ts <- as.numeric(sub$timestamp) / 3600   # hours since epoch
  n  <- length(x)
  sapply(radii, function(r) {
    vals <- vapply(seq_len(n - 1L), function(i) {
      for (j in seq.int(i + 1L, n)) {
        dj <- sqrt((x[j] - x[i])^2 + (y[j] - y[i])^2)
        if (dj >= r) {
          d_prev <- if (j == i + 1L) 0
                    else sqrt((x[j-1L]-x[i])^2 + (y[j-1L]-y[i])^2)
          step_d <- dj - d_prev
          frac   <- if (step_d > 0) (r - d_prev) / step_d else 0
          t_prev <- if (j > i + 1L) ts[j - 1L] - ts[i] else 0
          return(t_prev + frac * (ts[j] - ts[j - 1L]))
        }
      }
      NA_real_   # animal never exits — censored
    }, numeric(1))
    vals <- vals[!is.na(vals) & vals > 0]
    if (length(vals) < 2L) NA_real_ else var(log(vals))
  })
}

# ── Robust ctmm CI-table extraction ─────────────────────────────────────────
# ctmm's summary(fit)$CI matrix labels rows with Unicode Greek names —
# "τ[position]", "τ[velocity]" — not the literal text "tau position"/
# "tau velocity" that earlier code searched for, and the unit in parentheses
# is chosen automatically per fit (verified against this project's own saved
# CTMM_fit_*.rds objects: τ[position] is reported in days, τ[velocity] in
# hours, speed in km/day — none of which match a hard-coded "(hours)"/
# "(km/h)" suffix). A fixed-string lookup like par_ci["tau position (hours)",
# "est"] therefore raises a caught "subscript out of bounds" error and
# silently returns NA for every individual, regardless of unit. This helper
# matches by (ASCII-safe) pattern instead and converts to a fixed target unit
# so values are comparable across elephants. Returns NA when the parameter
# genuinely does not apply to the fitted model (e.g. OU models have no
# τ[velocity] or speed row at all — that is a correct NA, not a bug).
.TAU_TO_HOURS <- c(
  second = 1/3600, seconds = 1/3600,
  minute = 1/60,   minutes = 1/60,
  hour   = 1,      hours   = 1,
  day    = 24,     days    = 24,
  week   = 168,    weeks   = 168
)
.SPEED_TO_KMH <- c(
  "meter/second" = 3.6, "meters/second" = 3.6, "m/s" = 3.6,
  "meter/minute" = 0.06, "meters/minute" = 0.06,
  "meter/hour" = 0.001, "meters/hour" = 0.001,
  "meter/day" = 0.001/24, "meters/day" = 0.001/24,
  "kilometer/second" = 3600, "kilometers/second" = 3600,
  "kilometer/minute" = 60,   "kilometers/minute" = 60,
  "kilometer/hour" = 1, "kilometers/hour" = 1, "km/h" = 1, "km/hour" = 1,
  "kilometer/day" = 1/24, "kilometers/day" = 1/24, "km/day" = 1/24
)

extract_ctmm_ci <- function(ci_table, pattern, conversions) {
  rn  <- rownames(ci_table)
  idx <- grep(pattern, rn, ignore.case = TRUE)
  if (length(idx) == 0) return(NA_real_)
  row_name <- rn[idx[1]]
  val      <- unname(ci_table[idx[1], "est"])
  unit     <- sub(".*\\(([^)]+)\\).*", "\\1", row_name)
  factor   <- conversions[[unit]]
  if (is.null(factor)) {
    cat(sprintf("  [WARN] extract_ctmm_ci: unrecognised unit '%s' in row '%s' — left unconverted\n",
                unit, row_name))
    factor <- 1
  }
  val * factor
}


# ============================================================
# SECTION B — LOAD & CLEAN DATA
# ============================================================

cat("\n=== Loading and cleaning data ===\n")
# Screened from the raw collar files by build_screened_v6.R (revision, 2026-10-30);
# replaces all_elephants_combined.csv, whose 2-4 h gap rule dropped all hourly fixes.
df <- read.csv("data/processed/all_elephants_screened_v6.csv")
df$timestamp <- as.POSIXct(df$timestamp, format = "%m/%d/%Y %H:%M", tz = "UTC")

# IDs excluded from all analyses (wide seasonal migrant — methods not suited)
SKIP_IDS <- c(4557L)

df_clean <- df %>%
  filter(speed_km_h <= 5, !elephant_id %in% SKIP_IDS) %>%
  arrange(elephant_id, timestamp)

# 3905's collar is stationary (weekly SD ~40 m) from 2023-08-08 to the end of
# its record — dropped collar; those fixes are not movement (revision, 2026-10-29).
df_clean <- df_clean %>%
  filter(!(elephant_id == 3905 & timestamp >= as.POSIXct("2023-08-08", tz = "UTC")))

cat(sprintf("Excluded elephant IDs: %s\n", paste(SKIP_IDS, collapse = ", ")))

# Initialise global accumulation tables
model_comparison <- data.frame()
param_table      <- data.frame()
ci_df            <- data.frame()
akde_results     <- list()
akde_wet         <- list()
akde_dry         <- list()
ctmm_fits        <- list()   # best-fit ctmm models keyed by elephant ID string
etd_results      <- list()
etd_comparison   <- data.frame()


# ============================================================
# SECTION C — DATA QUALITY SUMMARY
# ============================================================

cat("\n=== Data Quality Summary ===\n")

qual <- do.call(rbind, lapply(unique(df_clean$elephant_id), function(id) {
  sub <- df_clean[df_clean$elephant_id == id, ]
  sub <- sub[order(sub$timestamp), ]
  dts <- as.numeric(diff(sub$timestamp), units = "secs") / 3600
  total_hrs <- as.numeric(diff(range(sub$timestamp)), units = "secs") / 3600
  data.frame(
    ID                 = id,
    n_fixes_raw        = nrow(sub),
    tracking_days      = round(total_hrs / 24, 1),
    expected_fixes     = round(total_hrs),
    fix_success_pct    = round(nrow(sub) / max(total_hrs, 1) * 100, 1),
    median_gap_hr      = round(median(dts, na.rm = TRUE), 2),
    max_gap_hr         = round(max(dts, na.rm = TRUE), 1),
    n_gaps_gt_6hr      = sum(dts > 6, na.rm = TRUE),
    n_gaps_gt_24hr     = sum(dts > 24, na.rm = TRUE),
    median_HDOP_m      = round(median(sub$gps_error_m, na.rm = TRUE), 1),
    pct_wet_season     = round(mean(is_wet(sub$timestamp)) * 100, 1),
    pct_dry_season     = round(mean(!is_wet(sub$timestamp)) * 100, 1)
  )
}))

tryCatch(
  write.csv(qual, file.path(QT, "Data_Quality_Summary.csv"), row.names = FALSE),
  error = function(e) {
    alt <- file.path(QT, paste0("Data_Quality_Summary_", format(Sys.time(), "%H%M%S"), ".csv"))
    write.csv(qual, alt, row.names = FALSE)
    cat("Saved as", basename(alt), "(original was locked)\n")
  }
)
cat("Data quality table saved → data_quality/tables/\n")

for (id in unique(df_clean$elephant_id)) {
  out_p <- file.path(QG, paste0("FixTimeline_", id, ".png"))
  sub <- df_clean[df_clean$elephant_id == id, ]
  sub$month <- as.Date(format(sub$timestamp, "%Y-%m-01"))
  monthly   <- sub %>% group_by(month) %>% summarise(n = n(), .groups = "drop")
  monthly$season <- ifelse(
    as.integer(format(monthly$month, "%m")) %in% c(11, 12, 1:4),
    "Wet (Nov-Apr)", "Dry (May-Oct)"
  )
  p <- ggplot(monthly, aes(x = month, y = n, fill = season)) +
    geom_col() +
    scale_fill_manual(values = c("Wet (Nov-Apr)" = "#1D9E75", "Dry (May-Oct)" = "#D85A30")) +
    labs(title = paste("Monthly GPS fix count — Elephant", id),
         x = "Month", y = "Number of fixes", fill = "Season") +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(out_p, p, width = 10, height = 4, dpi = PNG_DPI)
}
cat("Fix-rate timelines saved → data_quality/graphs/\n")


# ============================================================
# SECTION D — AKDE WITH VARIOGRAM & MODEL SELECTION
# ============================================================

cat("\n=== Running AKDE with variogram & model selection ===\n")

for (id in unique(df_clean$elephant_id)) {


  sub <- prep_sub(df_clean, id)
  if (nrow(sub) < 20) { cat("Skipping AKDE for", id, "(too few points)\n"); next }

  tel   <- make_tel(sub)
  guess <- ctmm.guess(tel, interactive = FALSE)
  vg    <- variogram(tel)

  fits <- ctmm.select(tel, guess, cores = N_CORES, verbose = TRUE)
  fit  <- fits[[1]]

  save_png({
    plot(vg, CTMM = fit,
         main     = paste("Empirical Variogram — Elephant", id),
         col.CTMM = c("#D85A30", "#888780"))
    legend("bottomright",
           legend = c(paste("Best model:", summary(fit)$name), "IID (null)"),
           col = c("#D85A30", "#888780"), lty = 1, lwd = 2, bty = "n", cex = 0.85)
  }, file.path(AG, paste0("Variogram_", id, ".png")))

  aic_vals    <- sapply(fits, function(f) tryCatch(f$AIC, error = function(e) NA_real_))
  model_names <- sapply(fits, function(f) tryCatch(summary(f)$name, error = function(e) class(f)[1]))
  delta_aic   <- aic_vals - min(aic_vals, na.rm = TRUE)
  aic_weights <- exp(-0.5 * delta_aic) / sum(exp(-0.5 * delta_aic), na.rm = TRUE)
  ms_table <- data.frame(
    Elephant = id, Model = model_names,
    AIC = round(aic_vals, 2), delta_AIC = round(delta_aic, 2),
    AIC_weight = round(aic_weights, 3), Best = seq_along(fits) == 1
  )
  write.csv(ms_table, file.path(AT, paste0("ModelSelection_", id, ".csv")), row.names = FALSE)

  # ── Asymptote & Isotropy diagnostic tests ──────────────────────────────────
  {
    # class(fit)[1] is always the generic "ctmm" S3 class for every model type
    # (OU, OUF, IID all share it) — the actual fitted model name ("OUF
    # anisotropic", "OU anisotropic", etc.) lives in summary(fit)$name, which
    # is what the model-selection table above (model_names, line ~476)
    # already correctly uses. Using class(fit)[1] here made the asymptote
    # test always evaluate against the literal string "ctmm", so is_finite
    # was FALSE for every individual regardless of which model actually won.
    best_class <- summary(fit)$name

    # Asymptote test — OUF/OU/OUm have finite τ_p → finite home range
    is_finite   <- grepl("^OU", best_class) && !grepl("BM|IID", best_class)
    asym_result <- if (is_finite) "Finite range (asymptote reached)" else "Non-finite range (IID/BM)"

    # τ_p, τ_v, σ² from the fitted model — see extract_ctmm_ci() in Section A
    # for why this can't be a hard-coded string lookup.
    par_ci_diag <- summary(fit)$CI
    tau_p_est <- extract_ctmm_ci(par_ci_diag, "position", .TAU_TO_HOURS)
    tau_v_est <- extract_ctmm_ci(par_ci_diag, "velocity", .TAU_TO_HOURS)
    sigma_est <- tryCatch(par_ci_diag["area (square kilometers)", "est"], error = function(e) NA_real_)

    # Isotropy test — does an anisotropic model beat the best isotropic one?
    m_classes  <- sapply(fits, function(f) tryCatch(class(f)[1], error = function(e) ""))
    aniso_idx  <- which(grepl("aniso", m_classes, ignore.case = TRUE))
    iso_idx    <- which(!grepl("aniso", m_classes, ignore.case = TRUE))
    if (length(aniso_idx) > 0 && length(iso_idx) > 0) {
      dAIC_iso <- aic_vals[aniso_idx[1]] - aic_vals[iso_idx[1]]
      iso_result <- if (dAIC_iso < -2)
        sprintf("Anisotropic preferred (ΔAIC = %.1f) — elongated home range", -dAIC_iso)
      else
        sprintf("Isotropic (anisotropic ΔAIC = +%.1f) — circular home range", dAIC_iso)
    } else {
      iso_result <- if (length(aniso_idx) > 0) "Only anisotropic models fitted"
                    else                        "Isotropic (no anisotropic model in candidate set)"
    }

    writeLines(c(
      paste("Elephant:", id),
      paste("Best model:      ", best_class),
      paste("Asymptote test:  ", asym_result),
      paste("Isotropy test:   ", iso_result),
      sprintf("tau_p (h):       %.3f", tau_p_est),
      sprintf("tau_v (h):       %.3f", tau_v_est),
      sprintf("sigma^2 (km2):   %.4f", sigma_est)
    ), file.path(AT, paste0("AKDE_Tests_", id, ".txt")))
  }

  tryCatch({
    par_ci  <- summary(fit)$CI
    ud_akde <- akde(tel, fit)
    Neff    <- summary(ud_akde)$DOF["area"]
    area_ci <- summary(ud_akde)$CI
    param_table <<- rbind(param_table, data.frame(
      ID        = id, Model = summary(fit)$name,
      tau_p_hr  = extract_ctmm_ci(par_ci, "position", .TAU_TO_HOURS),
      tau_v_hr  = extract_ctmm_ci(par_ci, "velocity", .TAU_TO_HOURS),
      speed_kmh = extract_ctmm_ci(par_ci, "^speed",   .SPEED_TO_KMH),
      sigma_km2 = area_ci[2], N_eff = round(Neff, 1),
      AKDE_low  = area_ci[1], AKDE_est = area_ci[2], AKDE_high = area_ci[3],
      CI_width  = area_ci[3] - area_ci[1], Reliable = Neff >= 30
    ))
  }, error = function(e) cat("  Parameter extraction failed for", id, ":", e$message, "\n"))

  ud_akde    <- akde(tel, fit)
  area_ci    <- summary(ud_akde)$CI
  area_50_ci <- summary(ud_akde, level.UD = 0.50)$CI
  tryCatch(saveRDS(ud_akde, file.path(AS, paste0("AKDE_", id, ".rds"))), error = function(e) NULL)
  tryCatch(saveRDS(fit,    file.path(AS, paste0("CTMM_fit_", id, ".rds"))), error = function(e) NULL)

  akde_results[[as.character(id)]] <- ud_akde
  ctmm_fits[[as.character(id)]]    <- fit
  model_comparison <- dplyr::bind_rows(model_comparison,
    data.frame(ID = id, Method = "AKDE", Area_km2 = area_ci[2], AIC = fit$AIC))
  ci_df <- rbind(ci_df, data.frame(
    ID = id, low = area_ci[1], est = area_ci[2], high = area_ci[3],
    CI_width = area_ci[3] - area_ci[1]
  ))

  ras      <- raster(ud_akde, DF = "CDF")
  tif_path <- file.path(AS, paste0("AKDE_", id, ".tif"))
  writeRaster(ras, filename = tif_path, overwrite = TRUE)

  # Extract 95% isopleth as shapefile in UTM 36S (needed for BA concordance)
  tryCatch({
    mask95 <- ras
    raster::values(mask95) <- ifelse(raster::values(ras) <= 0.95, 1L, NA_integer_)
    fact <- min(5L, max(1L, as.integer(max(nrow(mask95), ncol(mask95)) / 400)))
    if (fact > 1L) mask95 <- raster::aggregate(mask95, fact = fact, fun = max, na.rm = TRUE)
    akde_sf95 <- sf::st_as_sf(raster::rasterToPolygons(mask95, dissolve = TRUE))
    proj4_str <- raster::projection(ras)
    if (!is.null(proj4_str) && nchar(proj4_str) > 0)
      sf::st_crs(akde_sf95) <- proj4_str
    akde_sf95_utm <- sf::st_transform(akde_sf95, 32736)
    sf::st_write(akde_sf95_utm, file.path(AS, paste0("AKDE_95_", id, ".shp")),
                 delete_layer = TRUE, quiet = TRUE)
    cat(sprintf("  AKDE 95%% isopleth shapefile saved for %s\n", id))
  }, error = function(e) cat(sprintf("  AKDE 95%% shp failed for %s: %s\n", id, e$message)))

  save_png({
    plot(ud_akde,
         main     = paste("Elephant", id, "— AKDE Home Range (Full Dataset)"),
         level.UD = c(0.50, 0.95), level = 0.95)
    title(sub = sprintf("50%% area: %.1f km²  |  95%% area: %.1f km²",
                        area_50_ci[2], area_ci[2]), cex.sub = 0.8)
  }, file.path(AG, paste0("AKDE_", id, "_full.png")))

  sub_wet <- thin_to(sub[is_wet(sub$timestamp), ])
  sub_dry <- thin_to(sub[!is_wet(sub$timestamp), ])
  has_wet <- nrow(sub_wet) >= 20
  has_dry <- nrow(sub_dry) >= 20
  ud_wet <- ud_dry <- NULL
  area_wet <- area_dry <- NA

  if (has_wet) {
    tel_w  <- make_tel(sub_wet)
    fit_w  <- ctmm.select(tel_w, ctmm.guess(tel_w, interactive = FALSE),
                          cores = N_CORES, verbose = TRUE)
    ud_wet <- akde(tel_w, fit_w[[1]])
    area_wet <- summary(ud_wet)$CI[2]
    tryCatch(saveRDS(ud_wet, file.path(AS, paste0("AKDE_wet_", id, ".rds"))), error = function(e) NULL)
    akde_wet[[as.character(id)]] <- ud_wet
    model_comparison <- dplyr::bind_rows(model_comparison,
      data.frame(ID = id, Method = "AKDE_Wet", Area_km2 = area_wet, AIC = fit_w[[1]]$AIC))
  }

  if (has_dry) {
    tel_d  <- make_tel(sub_dry)
    fit_d  <- ctmm.select(tel_d, ctmm.guess(tel_d, interactive = FALSE),
                          cores = N_CORES, verbose = TRUE)
    ud_dry <- akde(tel_d, fit_d[[1]])
    area_dry <- summary(ud_dry)$CI[2]
    tryCatch(saveRDS(ud_dry, file.path(AS, paste0("AKDE_dry_", id, ".rds"))), error = function(e) NULL)
    akde_dry[[as.character(id)]] <- ud_dry
    model_comparison <- dplyr::bind_rows(model_comparison,
      data.frame(ID = id, Method = "AKDE_Dry", Area_km2 = area_dry, AIC = fit_d[[1]]$AIC))
  }

  save_png({
    par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
    if (has_wet) {
      plot(ud_wet, main = paste0("Wet Season AKDE — Elephant ", id, "\n(Nov-Apr)"),
           level.UD = c(0.50, 0.95), level = 0.95,
           col.UD = "steelblue", col.level = "navy")
      title(sub = sprintf("50%%/95%% areas: %.0f / %.0f km²",
                          summary(ud_wet, level.UD = 0.50)$CI[2], area_wet), cex.sub = 0.75)
    } else { plot.new(); title(main = paste0("Wet Season — Elephant ", id, "\n(insufficient data)")) }
    if (has_dry) {
      plot(ud_dry, main = paste0("Dry Season AKDE — Elephant ", id, "\n(May-Oct)"),
           level.UD = c(0.50, 0.95), level = 0.95,
           col.UD = "darkorange", col.level = "saddlebrown")
      title(sub = sprintf("50%%/95%% areas: %.0f / %.0f km²",
                          summary(ud_dry, level.UD = 0.50)$CI[2], area_dry), cex.sub = 0.75)
    } else { plot.new(); title(main = paste0("Dry Season — Elephant ", id, "\n(insufficient data)")) }
  }, file.path(AG, paste0("AKDE_", id, "_seasonal.png")), width = 14, height = 7)

  cat("Finished AKDE for Elephant", id, "\n")
}

if (nrow(param_table) > 0) {
  write.csv(param_table, file.path(AT, "AKDE_Parameters.csv"), row.names = FALSE)
  cat("AKDE parameters saved → AKDE/tables/\n")
}

if (nrow(ci_df) > 0) {
  p_ci <- ggplot(ci_df, aes(x = reorder(factor(ID), -est), y = est)) +
    geom_col(fill = "#1D9E75", alpha = 0.75, width = 0.6) +
    geom_errorbar(aes(ymin = low, ymax = high), width = 0.3, colour = "#D85A30", linewidth = 0.8) +
    annotate("text", x = Inf, y = Inf, label = "Error bars = 95% CI",
             hjust = 1.1, vjust = 1.5, size = 3, colour = "grey40") +
    labs(title = "AKDE 95% home range area with confidence intervals",
         x = "Elephant ID", y = expression("Area (km"^2*")")) +
    theme_bw(base_size = 12)
  ggsave(file.path(AG, "AKDE_CI_Summary.png"), p_ci, width = 10, height = 5, dpi = PNG_DPI)
  cat("CI summary plot saved → AKDE/graphs/\n")
}


# ============================================================
# SECTION E — dBBMM  (Manual of Applied Spatial Ecology approach)
# ============================================================

cat("\n=== Running dBBMM ===\n")

utm_proj_sp <- sp::CRS("+proj=utm +zone=36 +south +datum=WGS84 +units=m +no_defs")

for (id in unique(df_clean$elephant_id)) {


  sub <- prep_sub(df_clean, id, max_pts = Inf)
  # Restrict to 50-km core to remove extreme outliers
  if (nrow(sub) > 10) {
    cx <- median(sub$x_utm); cy <- median(sub$y_utm)
    d  <- sqrt((sub$x_utm - cx)^2 + (sub$y_utm - cy)^2) / 1000
    sub <- sub[d <= 50, ]
  }
  if (nrow(sub) < 20) { cat("Skipping dBBMM for", id, "(too few points)\n"); next }

  # Remove step-speed outliers (> 15 km/step)
  for (.k in seq_len(20)) {
    stp_km <- c(0, sqrt(diff(sub$x_utm)^2 + diff(sub$y_utm)^2) / 1000)
    bad <- stp_km > 15; bad[1] <- FALSE
    if (!any(bad)) break
    sub <- sub[!bad, ]
  }
  if (nrow(sub) < 20) { cat("Skipping dBBMM for", id, "(too few after speed filter)\n"); next }

  loc_err <- median(sub$gps_error_m, na.rm = TRUE)
  if (is.na(loc_err) || loc_err <= 0) loc_err <- 10

  # time.step must be < minimum fix interval (minutes)
  min_dt_min <- min(as.numeric(diff(sub$timestamp), units = "secs") / 60, na.rm = TRUE)
  time_step  <- max(5, min(60, floor(min_dt_min * 0.75)))

  # ext increases across retries too (not just margin/window.size). A fixed
  # ext = 1.5 is enough for most individuals, but fails with "Lower y grid not
  # large enough" for sparser, wider-spread point sets (e.g. elephant 4501,
  # whose sub-hourly fixes are disproportionately thinned by the speed_km_h <=
  # 5 filter in Section B) — confirmed empirically that ext = 3 resolves this
  # at every margin/window.size combination tried.
  dbbmm_params <- list(
    list(margin = 7, window.size = 19, ext = 1.5),
    list(margin = 5, window.size = 11, ext = 3),
    list(margin = 3, window.size = 7,  ext = 5)
  )
  dbbmm_obj         <- NULL
  dbbmm_mv          <- NULL
  dbbmm_params_used <- NULL   # which of the 3 tiers below actually converged

  for (.p in dbbmm_params) {
    sub_b <- prep_bursts(sub, min_burst_size = .p$window.size + 1)
    if (nrow(sub_b) < 20) next

    mv_try <- tryCatch(
      move::move(
        x      = sub_b$x_utm,
        y      = sub_b$y_utm,
        time   = sub_b$timestamp,
        animal = as.character(id),
        proj   = utm_proj_sp
      ),
      error = function(e) {
        cat(sprintf("  Elephant %s: move() failed — %s\n", id, e$message)); NULL
      })
    if (is.null(mv_try)) next

    cat(sprintf("  Elephant %s: min timeLag = %.1f min\n",
                id, min(move::timeLag(mv_try, units = "mins"), na.rm = TRUE)))

    obj_try <- tryCatch(
      move::brownian.bridge.dyn(
        object         = mv_try,
        location.error = loc_err,
        margin         = .p$margin,
        window.size    = .p$window.size,
        dimSize        = 200,
        time.step      = time_step,
        ext            = .p$ext
      ),
      error = function(e) {
        cat(sprintf("  Elephant %s: dBBMM failed (margin=%d, ws=%d) — %s\n",
                    id, .p$margin, .p$window.size, e$message))
        NULL
      })

    if (!is.null(obj_try)) {
      dbbmm_mv          <- mv_try
      dbbmm_obj         <- obj_try
      sub               <- sub_b
      dbbmm_params_used <- .p
      break
    }
  }

  if (is.null(dbbmm_obj)) {
    cat(sprintf("  Elephant %s: dBBMM failed all param sets — skipping\n", id)); next
  }

  # Volume UD (cumulative, 0–1) for raster output and area calculation
  cat_cont  <- move::getVolumeUD(dbbmm_obj)
  vud_vals  <- raster::values(cat_cont)
  if (all(is.na(vud_vals)) || diff(range(vud_vals, na.rm = TRUE)) < 0.01) {
    cat(sprintf("  Elephant %s: dBBMM UD degenerate — skipping\n", id)); next
  }
  vud_clean <- pmin(pmax(vud_vals, 0), 1)
  vud_clean[is.na(vud_clean)] <- 1
  raster::values(cat_cont) <- vud_clean

  cell_km2 <- prod(raster::res(cat_cont)) / 1e6
  area50   <- sum(vud_clean <= 0.50, na.rm = TRUE) * cell_km2
  area95   <- sum(vud_clean <= 0.95, na.rm = TRUE) * cell_km2

  # Extract isopleths via estUD (Manual of Applied Spatial Ecology method)
  shp50 <- .dbbmm_isopleth(dbbmm_obj, 50)
  shp80 <- .dbbmm_isopleth(dbbmm_obj, 80)
  shp90 <- .dbbmm_isopleth(dbbmm_obj, 90)
  shp95 <- .dbbmm_isopleth(dbbmm_obj, 95)
  shp99 <- .dbbmm_isopleth(dbbmm_obj, 99)

  # Recompute areas from shapefiles where available (more accurate than raster count)
  if (!is.null(shp95)) area95 <- sum(as.numeric(sf::st_area(shp95)), na.rm = TRUE) / 1e6
  if (!is.null(shp50)) area50 <- sum(as.numeric(sf::st_area(shp50)), na.rm = TRUE) / 1e6

  model_comparison <- dplyr::bind_rows(model_comparison,
    data.frame(ID = id, Method = "dBBMM", Area_km2 = area95, AIC = NA,
               dBBMM_margin      = dbbmm_params_used$margin,
               dBBMM_window_size = dbbmm_params_used$window.size,
               dBBMM_ext         = dbbmm_params_used$ext))

  # Save volume-UD raster
  tif_path <- file.path(DS, paste0("dBBMM_", id, ".tif"))
  raster::writeRaster(cat_cont, filename = tif_path, overwrite = TRUE)

  # Save isopleths as shapefiles
  for (pct in c("50","80","90","95","99")) {
    shp_obj <- get(paste0("shp", pct))
    if (!is.null(shp_obj))
      tryCatch(
        sf::st_write(shp_obj, file.path(DS, paste0("dBBMM_", pct, "_", id, ".shp")),
                     delete_layer = TRUE, quiet = TRUE),
        error = function(e) NULL)
  }

  # Zoom extent for the 99% isopleth (shared by Plots 1, 3, 4 below), with padding
  ext99 <- .pct_extent(cat_cont, 0.99)

  # Plot 1: volume-UD heatmap with contours — zoomed to 99%, masked outside it,
  # 95% line recoloured (blue) so it no longer blends into the red core of the heatmap
  cat_cont_99 <- cat_cont
  cat_cont_99[cat_cont_99 > 0.99] <- NA
  cat_cont_99 <- raster::crop(cat_cont_99, ext99)
  pts_in_99   <- sub$x_utm >= ext99@xmin & sub$x_utm <= ext99@xmax &
                 sub$y_utm >= ext99@ymin & sub$y_utm <= ext99@ymax

  save_png({
    raster::plot(cat_cont_99, col = rev(heat.colors(100)), colNA = NA,
                 main = paste("dBBMM Utilization Distribution — Elephant", id),
                 xlab = "Easting (m)", ylab = "Northing (m)")
    raster::contour(cat_cont_99,
                    levels = c(0.50, 0.90, 0.95, 0.99), add = TRUE,
                    col = c("yellow","white","blue","grey80"), lwd = c(2,1.5,1.5,1), labcex = 0.8)
    points(sub$x_utm[pts_in_99], sub$y_utm[pts_in_99], pch = 20, cex = 0.3, col = "black")
    legend("bottomright",
           legend = c("50%","90%","95%","99%",
                      sprintf("50%% area: %.2f km²", area50),
                      sprintf("95%% area: %.2f km²", area95)),
           lty = c(1,1,1,1,NA,NA), col = c("yellow","white","blue","grey80",NA,NA),
           lwd = c(2,1.5,1.5,1,NA,NA), bty = "n", cex = 0.75)
  }, file.path(DG, paste0("dBBMM_", id, "_ud.png")))

  # Plot 2: isopleth polygons (reference-style 2×2 combined plot)
  save_png({
    par(mfrow = c(1, 1))
    base <- if (!is.null(shp99)) shp99 else if (!is.null(shp95)) shp95 else shp90
    if (is.null(base)) { plot.new(); title(main = "No isopleths") } else {
      plot(sf::st_geometry(base),
           main = paste("dBBMM Isopleths — Elephant", id),
           xlab = "Easting (m)", ylab = "Northing (m)")
      if (!is.null(shp95)) plot(sf::st_geometry(shp95), add = TRUE, border = "black", lwd = 1.5)
      if (!is.null(shp90)) plot(sf::st_geometry(shp90), add = TRUE, border = "blue",  lwd = 1.5)
      if (!is.null(shp80)) plot(sf::st_geometry(shp80), add = TRUE, border = "red",   lwd = 1.5)
      if (!is.null(shp50)) plot(sf::st_geometry(shp50), add = TRUE, border = "green", lwd = 2)
      points(sub$x_utm, sub$y_utm, pch = 1, cex = 0.5)
      legend("bottomright",
             legend = c("99%","95%","90%","80%","50%",
                        sprintf("50%% area: %.2f km²", area50),
                        sprintf("95%% area: %.2f km²", area95)),
             lty = c(1,1,1,1,1,NA,NA),
             col = c("black","black","blue","red","green",NA,NA),
             bty = "n", cex = 0.8)
    }
  }, file.path(DG, paste0("dBBMM_", id, "_isopleths.png")))

  # Plot 3: native move dBBMM UD heatmap (raw density), zoomed to 99% and masked
  # outside it so no colour is drawn beyond the meaningful UD area
  dbbmm_obj_99 <- dbbmm_obj
  dbbmm_obj_99[cat_cont > 0.99] <- NA
  dbbmm_obj_99 <- raster::crop(dbbmm_obj_99, ext99)

  save_png({
    raster::plot(dbbmm_obj_99, colNA = NA,
         main = paste("dBBMM Utilization Distribution (native) — Elephant", id),
         xlab = "Easting (m)", ylab = "Northing (m)")
  }, file.path(DG, paste0("dBBMM_", id, "_native.png")))

  # Plot 4: contour lines at 50%, 90%, 95%, 99%, zoomed to the 99% extent.
  # NOTE: previously drawn on the raw density raster (dbbmm_obj), whose values
  # are on a tiny 0-0.04-ish PDF scale — the requested 0.5/0.9/0.95/0.99 levels
  # never occurred in that data, so this plot was silently rendering with no
  # contour lines at all. Switched to cat_cont (the 0-1 cumulative-volume
  # raster used correctly in Plot 1) so the isopleth levels are meaningful.
  cat_cont_crop <- raster::crop(cat_cont, ext99)
  save_png({
    contour(cat_cont_crop, levels = c(0.5, 0.9, 0.95, 0.99),
            main  = paste("dBBMM Contour Lines — Elephant", id),
            xlab  = "Easting (m)", ylab = "Northing (m)",
            col   = c("blue", "green", "orange", "red"),
            lwd   = c(2, 1.5, 1.5, 1),
            labcex = 0.85)
    legend("bottomright",
           legend = c("50%", "90%", "95%", "99%"),
           col    = c("blue", "green", "orange", "red"),
           lwd    = c(2, 1.5, 1.5, 1), bty = "n", cex = 0.8)
  }, file.path(DG, paste0("dBBMM_", id, "_contour.png")))

  # Plot 5: GPS movement path (move object, type="o" as in reference)
  save_png({
    plot(dbbmm_mv, type = "o",
         col  = "#1D9E75", lwd = 1.5, pch = 20, cex = 0.5,
         main = paste("GPS Movement Path — Elephant", id),
         xlab = "Easting (m)", ylab = "Northing (m)")
  }, file.path(DG, paste0("dBBMM_", id, "_track.png")))

  # ── Behavioral Variance (σ_m²) time series + breakpoint plot ──────────────
  {
    # Try native move accessor; fall back to step-speed proxy from GPS track
    bmv_raw <- tryCatch(move::brownian.motion.variance.dyn(dbbmm_obj), error = function(e) NULL)

    if (!is.null(bmv_raw) && length(bmv_raw) > 1) {
      bmv_vals <- as.numeric(bmv_raw)
      bmv_t    <- seq_along(bmv_vals)
      x_lab    <- "Time step (window)"
    } else {
      dx  <- diff(sub$x_utm);  dy <- diff(sub$y_utm)
      dt_h <- pmax(as.numeric(diff(sub$timestamp), units = "secs") / 3600, 1e-3)
      bmv_vals <- (sqrt(dx^2 + dy^2) / dt_h)^2     # (speed m/h)^2 ≈ σ_m² proxy
      bmv_t    <- sub$timestamp[-1]
      x_lab    <- "Date"
    }

    thr_val <- median(bmv_vals, na.rm = TRUE)
    high_mv <- bmv_vals > 2 * thr_val

    save_png({
      par(mar = c(5, 5, 4, 2))
      plot(bmv_t, bmv_vals / 1e6, type = "l", col = "#1D9E75", lwd = 1,
           main  = paste("Behavioral Variance σᵐ² — Elephant", id),
           xlab  = x_lab,
           ylab  = expression("Variance (m"^2*"/h) × 10"^6))
      abline(h = thr_val / 1e6, lty = 2, col = "#D85A30", lwd = 1.5)
      if (any(high_mv, na.rm = TRUE))
        points(bmv_t[high_mv], bmv_vals[high_mv] / 1e6,
               pch = 16, col = "#D85A30", cex = 0.45)
      legend("topright",
             legend = c("Variance", "Median (rest/travel boundary)", "High activity (>2× median)"),
             lty = c(1, 2, NA), pch = c(NA, NA, 16),
             col = c("#1D9E75", "#D85A30", "#D85A30"), bty = "n", cex = 0.75)
    }, file.path(DG, paste0("dBBMM_", id, "_variance.png")))
  }

  # ── Sensitivity analysis — how margin + window.size affect 95% area ────────
  {
    sens_sets <- list(
      list(margin = 3, window.size = 7),
      list(margin = 5, window.size = 11),
      list(margin = 7, window.size = 19),
      list(margin = 9, window.size = 29)
    )
    sens_rows <- lapply(sens_sets, function(.sp) {
      sub_s <- prep_bursts(sub, min_burst_size = .sp$window.size + 1)
      if (nrow(sub_s) < 20)
        return(data.frame(margin = .sp$margin, window_size = .sp$window.size, area95_km2 = NA_real_))
      mv_s <- tryCatch(
        move::move(x = sub_s$x_utm, y = sub_s$y_utm, time = sub_s$timestamp,
                   animal = as.character(id), proj = utm_proj_sp),
        error = function(e) NULL)
      # Same "Lower y grid not large enough" issue as the main computation can
      # hit any individual margin/window combo here too (e.g. sparser
      # individuals like 4501) — retry with progressively larger ext instead
      # of a single fixed value, so one bad combo doesn't just silently NA out.
      obj_s <- NULL
      if (!is.null(mv_s)) {
        for (.ext_try in c(1.5, 3, 5, 8)) {
          obj_s <- tryCatch(
            move::brownian.bridge.dyn(object = mv_s, location.error = loc_err,
                                      margin = .sp$margin, window.size = .sp$window.size,
                                      dimSize = 150, time.step = time_step, ext = .ext_try),
            error = function(e) NULL)
          if (!is.null(obj_s)) break
        }
      }
      a_s <- if (!is.null(obj_s)) {
        sh_s <- .dbbmm_isopleth(obj_s, 95)
        if (!is.null(sh_s)) sum(as.numeric(sf::st_area(sh_s)), na.rm = TRUE) / 1e6 else NA_real_
      } else NA_real_
      data.frame(margin = .sp$margin, window_size = .sp$window.size, area95_km2 = a_s)
    })
    sens_df <- do.call(rbind, sens_rows)
    write.csv(sens_df, file.path(DS, paste0("dBBMM_Sensitivity_", id, ".csv")), row.names = FALSE)

    save_png({
      lbl  <- paste0("m=", sens_df$margin, "\nws=", sens_df$window_size)
      cols <- c("#888780", "#1D9E75", "#D85A30", "#5A6BBA")
      ym   <- max(c(sens_df$area95_km2, area95), na.rm = TRUE) * 1.25
      bp   <- barplot(sens_df$area95_km2, names.arg = lbl, col = cols,
                      ylim = c(0, ym), xlab = "Parameter set (margin / window size)",
                      ylab = expression("95% Home Range (km"^2*")"),
                      main = paste("dBBMM Sensitivity Analysis — Elephant", id))
      text(bp, sens_df$area95_km2 + ym * 0.03,
           ifelse(is.na(sens_df$area95_km2), "NA", sprintf("%.1f", sens_df$area95_km2)),
           cex = 0.8)
      abline(h = area95, lty = 2, col = "red", lwd = 1.5)
      legend("topright", "Chosen set (95% area)", lty = 2, col = "red", bty = "n", cex = 0.8)
    }, file.path(DG, paste0("dBBMM_", id, "_sensitivity.png")))
  }

  cat(sprintf("Finished dBBMM for Elephant %s — 50%%: %.1f km²  95%%: %.1f km²\n",
              id, area50, area95))
}


# ============================================================
# SECTION F — aLOCOH
# ============================================================

cat("\n=== Running aLOCOH ===\n")

# getverticeshr() can fail on some adehabitatHR versions/inputs; this is the
# fallback used to pull the nearest-percent polygon directly. Defined once
# here (not inside the loop) so it exists before its first use below.
extract_pct <- function(spdf, target) {
  if (!inherits(spdf, "SpatialPolygonsDataFrame")) return(NULL)
  if (!"percent" %in% names(spdf)) return(NULL)
  pcts <- as.numeric(as.character(spdf$percent))
  ok   <- !is.na(pcts)
  if (!any(ok)) return(NULL)
  spdf[ok, ][which.min(abs(pcts[ok] - target)), ]
}

for (id in unique(df_clean$elephant_id)) {


  sub <- prep_sub(df_clean, id)
  if (nrow(sub) < 20) { cat("Skipping aLOCOH for", id, "(too few points)\n"); next }

  sp_coords <- SpatialPoints(sub[, c("x_utm", "y_utm")], proj4string = utm_proj)

  med_step <- if ("speed_m_hr" %in% names(sub))
    median(sub$speed_m_hr, na.rm = TRUE) * 3
  else
    median(sub$speed_km_h * 1000, na.rm = TRUE) * 3
  a_param <- ifelse(!is.na(med_step) & med_step > 0, med_step * 10, 5000)

  # ── Parameter selection plot: 95% area vs. a value ──────────────────────────
  # Use a thinned copy (≤500 pts) so the 10 trial fits stay fast.
  # getverticeshr can fail; fall back to extract_pct (same as main computation).
  {
    thin_n    <- min(500L, nrow(sub))
    sub_thin  <- sub[sort(sample(nrow(sub), thin_n)), ]
    sp_thin   <- SpatialPoints(sub_thin[, c("x_utm", "y_utm")], proj4string = utm_proj)

    a_seq   <- unique(round(exp(seq(log(max(100, a_param * 0.20)),
                                   log(a_param * 4), length.out = 10))))
    a_areas <- vapply(a_seq, function(a_try) {
      ud_try <- tryCatch(
        LoCoH.a(sp_thin, a = a_try, unin = "m", unout = "km2"),
        error = function(e) NULL)
      if (is.null(ud_try)) return(NA_real_)
      sp95 <- tryCatch(getverticeshr(ud_try, percent = 95),
                       error = function(e) extract_pct(ud_try, 95))
      if (is.null(sp95)) return(NA_real_)
      tryCatch(as.numeric(sum(sf::st_area(sf::st_as_sf(sp95)), na.rm = TRUE)) / 1e6,
               error = function(e) NA_real_)
    }, numeric(1))

    valid_a <- !is.na(a_areas)
    if (sum(valid_a) >= 2) {
      save_png({
        y_max <- max(a_areas[valid_a]) * 1.25
        plot(a_seq[valid_a] / 1000, a_areas[valid_a],
             type = "b", pch = 19, col = "#1D9E75", lwd = 2,
             main = paste("aLOCOH Parameter Selection — Elephant", id),
             xlab = "a value (km)", ylab = expression("95% HR Area (km"^2*")"),
             ylim = c(0, y_max))
        abline(v = a_param / 1000, lty = 2, col = "#D85A30", lwd = 1.5)
        legend("bottomright",
               legend = c(sprintf("Area (n=%d subsample)", thin_n),
                          sprintf("Chosen a = %.0f m", a_param)),
               lty = c(NA, 2), pch = c(19, NA),
               col = c("#1D9E75", "#D85A30"), bty = "n", cex = 0.85)
      }, file.path(LG, paste0("aLOCOH_", id, "_param_selection.png")))
    } else {
      save_png({
        plot.new()
        title(main = paste("aLOCOH Parameter Selection — Elephant", id),
              sub  = "Insufficient valid a-values to plot curve")
      }, file.path(LG, paste0("aLOCOH_", id, "_param_selection.png")))
      cat(sprintf("  Elephant %s: aLOCOH param-selection — only %d valid trial(s)\n",
                  id, sum(valid_a)))
    }
  }

  ud_alocoh <- tryCatch(
    LoCoH.a(sp_coords, a = a_param, unin = "m", unout = "km2"),
    error = function(e) {
      cat(sprintf("  Elephant %s: LoCoH.a failed — %s\n", id, conditionMessage(e)))
      NULL
    }
  )
  if (is.null(ud_alocoh)) next

  spoly50 <- tryCatch(getverticeshr(ud_alocoh, percent = 50),
                      error = function(e) extract_pct(ud_alocoh, 50))
  spoly95 <- tryCatch(getverticeshr(ud_alocoh, percent = 95),
                      error = function(e) extract_pct(ud_alocoh, 95))

  if (is.null(spoly50) && is.null(spoly95)) {
    cat(sprintf("  Elephant %s: aLOCOH — could not extract polygons\n", id)); next
  }

  to_sf <- function(sp) {
    if (is.null(sp)) return(NULL)
    out <- st_as_sf(sp); st_crs(out) <- 32736; out
  }
  sf50 <- to_sf(spoly50);  sf95 <- to_sf(spoly95)
  area50_al <- if (!is.null(sf50)) as.numeric(st_area(sf50)) / 1e6 else NA
  area95_al <- if (!is.null(sf95)) as.numeric(st_area(sf95)) / 1e6 else NA

  model_comparison <- dplyr::bind_rows(model_comparison,
    data.frame(ID = id, Method = "aLOCOH", Area_km2 = area95_al, AIC = NA))

  if (!is.null(sf50)) st_write(sf50, file.path(LS, paste0("aLOCOH_50_", id, ".shp")),
                               delete_layer = TRUE, quiet = TRUE)
  if (!is.null(sf95)) st_write(sf95, file.path(LS, paste0("aLOCOH_95_", id, ".shp")),
                               delete_layer = TRUE, quiet = TRUE)

  base_sf <- if (!is.null(sf95)) sf95 else sf50
  save_png({
    plot(st_geometry(base_sf),
         main = paste("aLOCOH Isopleths — Elephant", id),
         col = adjustcolor("forestgreen", alpha.f = 0.15),
         border = "darkgreen", lwd = 2,
         xlab = "Easting (m)", ylab = "Northing (m)")
    if (!is.null(sf50))
      plot(st_geometry(sf50), add = TRUE,
           col = adjustcolor("darkgreen", alpha.f = 0.4), border = "darkgreen", lwd = 2)
    points(sub$x_utm, sub$y_utm, pch = 20, cex = 0.3, col = "black")
    legend("bottomright",
           legend = c(sprintf("95%% hull: %.2f km²", area95_al),
                      sprintf("50%% core: %.2f km²", area50_al), "GPS fixes"),
           fill = c(adjustcolor("forestgreen", 0.15), adjustcolor("darkgreen", 0.4), NA),
           pch = c(NA, NA, 20), col = c(NA, NA, "black"), bty = "n")
  }, file.path(LG, paste0("aLOCOH_", id, ".png")))

  cat("Finished aLOCOH for Elephant", id, "\n")
}


# ============================================================
# SECTION G — TGD
# ============================================================

cat("\n=== Computing TGD ===\n")

for (id in unique(df_clean$elephant_id)) {


  sub <- prep_sub(df_clean, id, max_pts = Inf)
  if (nrow(sub) > 10) {
    cx <- median(sub$x_utm); cy <- median(sub$y_utm)
    d  <- sqrt((sub$x_utm - cx)^2 + (sub$y_utm - cy)^2) / 1000
    sub <- sub[d <= 50, ]
  }
  if (nrow(sub) < 20) { cat("Skipping TGD for", id, "(too few points)\n"); next }

  sub$dt <- c(NA, as.numeric(difftime(
    sub$timestamp[-1], sub$timestamp[-nrow(sub)], units = "hours")))
  sub <- sub[!is.na(sub$dt) & sub$dt > 0, ]
  v_max <- 4000  # m/h

  ellipses <- lapply(1:(nrow(sub) - 1), function(i) {
    p1 <- c(sub$x_utm[i], sub$y_utm[i])
    p2 <- c(sub$x_utm[i + 1], sub$y_utm[i + 1])
    dt <- sub$dt[i + 1]
    if (dt > 24) return(NULL)
    d     <- sqrt(sum((p2 - p1)^2))
    r_max <- v_max * dt
    if (d > r_max) return(NULL)
    line        <- st_linestring(rbind(p1, p2))
    buffer_dist <- min((r_max - d) / 2, 5000)
    st_buffer(st_sfc(line, crs = 32736), dist = buffer_dist)
  })

  valid <- ellipses[!sapply(ellipses, is.null)]
  if (length(valid) > 0) {
    tgd_poly <- st_union(st_sf(geometry = do.call(c, valid)))
    area     <- as.numeric(st_area(tgd_poly)) / 1e6
    model_comparison <- dplyr::bind_rows(model_comparison,
      data.frame(ID = id, Method = "TGD", Area_km2 = area, AIC = NA))
    st_write(tgd_poly, file.path(TS, paste0("TGD_", id, ".shp")),
             delete_layer = TRUE, quiet = TRUE)
    save_png({
      plot(st_geometry(tgd_poly),
           main = paste("TGD Potential Path Area — Elephant", id),
           col = adjustcolor("steelblue", alpha.f = 0.4), border = "steelblue",
           xlab = "Easting (m)", ylab = "Northing (m)")
      lines(sub$x_utm, sub$y_utm, col = "grey30", lwd = 0.6)
      points(sub$x_utm, sub$y_utm, pch = 20, cex = 0.4, col = "black")
      legend("bottomright",
             legend = c(sprintf("Path area: %.1f km² | v_max=%d m/h", area, v_max),
                        "Track", "GPS fixes"),
             fill = c(adjustcolor("steelblue", 0.4), NA, NA),
             lty = c(NA,1,NA), pch = c(NA,NA,20), col = c(NA,"grey30","black"), bty = "n")
    }, file.path(TG, paste0("TGD_", id, ".png")))

    # ── Diel accessibility curve: PPA buffer area by hour of day ──────────────
    sub$hour <- as.integer(format(sub$timestamp, "%H", tz = "Africa/Harare"))
    diel_area <- tapply(seq_len(nrow(sub) - 1), sub$hour[-nrow(sub)], function(idx) {
      ells <- lapply(idx, function(i) {
        p1 <- c(sub$x_utm[i], sub$y_utm[i]); p2 <- c(sub$x_utm[i+1], sub$y_utm[i+1])
        dt <- sub$dt[i + 1]; if (is.na(dt) || dt <= 0 || dt > 24) return(NULL)
        d <- sqrt(sum((p2 - p1)^2)); r_max <- v_max * dt
        if (d > r_max) return(NULL)
        st_buffer(st_sfc(st_linestring(rbind(p1, p2)), crs = 32736),
                  dist = min((r_max - d) / 2, 5000))
      })
      ells <- ells[!sapply(ells, is.null)]
      if (length(ells) == 0) return(NA_real_)
      as.numeric(st_area(st_union(st_sf(geometry = do.call(c, ells))))) / 1e6
    })
    diel_df <- data.frame(hour = as.integer(names(diel_area)),
                          area_km2 = as.numeric(diel_area))
    diel_df  <- diel_df[!is.na(diel_df$area_km2), ]

    if (nrow(diel_df) >= 4) {
      save_png({
        plot(diel_df$hour, diel_df$area_km2, type = "b", pch = 19,
             col = "#1D9E75", lwd = 2,
             main  = paste("TGD Diel Accessibility — Elephant", id),
             xlab  = "Hour of day (local)", ylab = expression("Accessible area (km"^2*")"),
             xaxt  = "n")
        axis(1, at = seq(0, 23, 3), labels = paste0(seq(0, 23, 3), "h"))
        abline(v = c(6, 18), lty = 2, col = "orange", lwd = 1)
        legend("topright", c("PPA area per hour", "Dawn / Dusk"),
               lty = c(NA, 2), pch = c(19, NA),
               col = c("#1D9E75", "orange"), bty = "n", cex = 0.85)
      }, file.path(TG, paste0("TGD_", id, "_diel.png")))
    }

    # ── Speed threshold sensitivity: how v_max alters PPA area ────────────────
    v_seq   <- c(1000, 2000, 4000, 8000, 16000)   # m/h
    v_areas <- vapply(v_seq, function(vm) {
      ells2 <- lapply(1:(nrow(sub) - 1), function(i) {
        p1 <- c(sub$x_utm[i], sub$y_utm[i]); p2 <- c(sub$x_utm[i+1], sub$y_utm[i+1])
        dt <- sub$dt[i + 1]; if (is.na(dt) || dt <= 0 || dt > 24) return(NULL)
        d <- sqrt(sum((p2 - p1)^2)); r_max <- vm * dt
        if (d > r_max) return(NULL)
        st_buffer(st_sfc(st_linestring(rbind(p1, p2)), crs = 32736),
                  dist = min((r_max - d) / 2, 5000))
      })
      valid2 <- ells2[!sapply(ells2, is.null)]
      if (length(valid2) == 0) return(NA_real_)
      as.numeric(st_area(st_union(st_sf(geometry = do.call(c, valid2))))) / 1e6
    }, numeric(1))

    save_png({
      plot(v_seq / 1000, v_areas, type = "b", pch = 19, col = "#5A6BBA", lwd = 2,
           main  = paste("TGD Speed Threshold Sensitivity — Elephant", id),
           xlab  = "v_max (km/h)", ylab = expression("PPA Area (km"^2*")"),
           log   = "x")
      abline(v = v_max / 1000, lty = 2, col = "#D85A30", lwd = 1.5)
      legend("topleft", c("PPA area", sprintf("Chosen v_max = %d m/h", v_max)),
             lty = c(NA, 2), pch = c(19, NA),
             col = c("#5A6BBA", "#D85A30"), bty = "n", cex = 0.85)
    }, file.path(TG, paste0("TGD_", id, "_speed_sensitivity.png")))
  }
  cat("Finished TGD for Elephant", id, "\n")
}


# ============================================================
# SECTION G2 — ETD  (Elliptical Time Density)
# ============================================================

cat("\n=== Running ETD ===\n")

etd_results    <- list()
etd_comparison <- data.frame()

for (id in unique(df_clean$elephant_id)) {


  sub <- prep_sub(df_clean, id, max_pts = Inf)
  if (nrow(sub) < 10) { cat("Skipping ETD for", id, "(too few points)\n"); next }

  sens <- sigma_sensitivity(sub, sigmas = seq(250, 2000, 250))
  tryCatch(write.csv(sens, file.path(ET, paste0("ETD_sigma_sensitivity_", id, ".csv")),
                     row.names = FALSE), error = function(e) NULL)

  # sigma fixed for all animals (revision, 2026-10-28). With the corrected
  # sensitivity curve, area grows ~linearly in sigma with no plateau, so an
  # elbow rule only returned the grid's second value; the sensitivity curve is
  # kept as a diagnostic instead.
  sigma_use <- ETD_SIGMA
  cat(sprintf("  Elephant %s: sigma selected = %d m\n", id, sigma_use))

  ud <- tryCatch(compute_etd(sub, sigma = sigma_use, res = 250),
                 error = function(e) { cat("  compute_etd failed:", e$message, "\n"); NULL })
  if (is.null(ud)) next

  etd_results[[as.character(id)]] <- ud

  iso50 <- tryCatch(etd_isopleth(ud, 0.50), error = function(e) NULL)
  iso95 <- tryCatch(etd_isopleth(ud, 0.95), error = function(e) NULL)
  area50_e <- if (!is.null(iso50)) as.numeric(sf::st_area(iso50)) / 1e6 else NA_real_
  area95_e <- if (!is.null(iso95)) as.numeric(sf::st_area(iso95)) / 1e6 else NA_real_

  # ── Mean ellipse eccentricity ─────────────────────────────────────────────
  # For each movement segment: semi-major a = dist/2 + sigma, semi-minor b = sigma
  # Eccentricity e = sqrt(1 - (b/a)^2); e=0 → circle, e→1 → highly elongated
  ecc_vals <- tryCatch({
    sub_o <- sub[order(sub$timestamp), ]
    dx_e  <- diff(sub_o$x_utm); dy_e <- diff(sub_o$y_utm)
    dist_e <- sqrt(dx_e^2 + dy_e^2)
    a_e    <- dist_e / 2 + sigma_use
    b_e    <- rep(sigma_use, length(a_e))
    ecc    <- sqrt(pmax(0, 1 - (b_e / a_e)^2))
    ecc[is.finite(ecc)]
  }, error = function(e) numeric(0))
  mean_ecc <- if (length(ecc_vals) > 0) mean(ecc_vals, na.rm = TRUE) else NA_real_

  etd_comparison <- rbind(etd_comparison, data.frame(
    ID = id, Method = "ETD", sigma_m = sigma_use, Area_50 = area50_e, Area_95 = area95_e,
    Mean_Eccentricity = mean_ecc
  ))
  model_comparison <- dplyr::bind_rows(model_comparison,
    data.frame(ID = id, Method = "ETD", Area_km2 = area95_e, AIC = NA))

  sub$step_km <- c(NA_real_, sqrt(diff(sub$x_utm)^2 + diff(sub$y_utm)^2) / 1000)
  sub$dt_hr   <- c(NA_real_, as.numeric(diff(sub$timestamp), units = "secs") / 3600)
  sub$speed_c <- sub$step_km / sub$dt_hr
  sub_slow <- sub[!is.na(sub$speed_c) & sub$speed_c <  1, ]
  sub_fast <- sub[!is.na(sub$speed_c) & sub$speed_c >= 2, ]
  ud_slow  <- if (nrow(sub_slow) > 5)
    tryCatch(compute_etd(sub_slow, sigma_use), error = function(e) NULL) else NULL
  ud_fast  <- if (nrow(sub_fast) > 5)
    tryCatch(compute_etd(sub_fast, sigma_use), error = function(e) NULL) else NULL

  tryCatch(raster::writeRaster(ud,
    file.path(ES, paste0("ETD_", id, "_full.tif")), overwrite = TRUE), error = function(e) NULL)
  if (!is.null(ud_slow))
    tryCatch(raster::writeRaster(ud_slow,
      file.path(ES, paste0("ETD_", id, "_slow.tif")), overwrite = TRUE), error = function(e) NULL)
  if (!is.null(ud_fast))
    tryCatch(raster::writeRaster(ud_fast,
      file.path(ES, paste0("ETD_", id, "_fast.tif")), overwrite = TRUE), error = function(e) NULL)

  if (!is.null(iso50)) sf::st_write(iso50, file.path(ES, paste0("ETD_50_", id, ".shp")),
                                    delete_layer = TRUE, quiet = TRUE)
  if (!is.null(iso95)) sf::st_write(iso95, file.path(ES, paste0("ETD_95_", id, ".shp")),
                                    delete_layer = TRUE, quiet = TRUE)

  save_png({
    raster::plot(ud, col = rev(heat.colors(100)),
                 main = paste("ETD Utilization Distribution — Elephant", id),
                 xlab = "Easting (m)", ylab = "Northing (m)")
    if (!is.null(iso95)) plot(sf::st_geometry(iso95), add = TRUE, border = "white",  lwd = 2)
    if (!is.null(iso50)) plot(sf::st_geometry(iso50), add = TRUE, border = "yellow", lwd = 2)
    points(sub$x_utm, sub$y_utm, pch = 20, cex = 0.3, col = "black")
    legend("bottomright",
           legend = c("95% HR", "50% Core",
                      sprintf("50%%: %.1f km²", area50_e),
                      sprintf("95%%: %.1f km²", area95_e),
                      sprintf("sigma = %d m", sigma_use)),
           lty = c(1,1,NA,NA,NA), col = c("white","yellow",NA,NA,NA),
           lwd = c(2,2,NA,NA,NA), bty = "n", cex = 0.8)
  }, file.path(EG, paste0("ETD_", id, "_full.png")))

  if (!is.null(ud_slow) && !is.null(ud_fast)) {
    save_png({
      par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
      raster::plot(ud_slow, col = rev(heat.colors(100)),
                   main = "Slow steps (< 1 km/h)\nResting / foraging")
      raster::plot(ud_fast, col = rev(heat.colors(100)),
                   main = "Fast steps (>= 2 km/h)\nTransit / directional")
    }, file.path(EG, paste0("ETD_", id, "_speed_stratified.png")), width = 14, height = 7)
  }

  save_png({
    plot(sens$sigma_m, sens$area_km2, type = "b", pch = 19, col = "#1D9E75",
         xlab = "sigma (m)", ylab = "95% area (km²)",
         main = paste("ETD sigma sensitivity — Elephant", id))
    abline(v = sigma_use, col = "#D85A30", lty = 2)
    legend("bottomright", legend = sprintf("Selected sigma = %d m", sigma_use),
           lty = 2, col = "#D85A30", bty = "n")
  }, file.path(EG, paste0("ETD_sigma_sensitivity_", id, ".png")))

  cat("Finished ETD for Elephant", id, "\n")
}

tryCatch(write.csv(etd_comparison, file.path(ET, "ETD_Comparison.csv"), row.names = FALSE),
         error = function(e) NULL)
if (nrow(etd_comparison) > 0)
  cat("ETD summary saved → ETD/tables/\n")


# ============================================================
# SECTION H — SEASONAL COMPARISON
# ============================================================

cat("\n=== Seasonal comparison ===\n")

seasonal_df <- do.call(rbind, lapply(unique(df_clean$elephant_id), function(id) {
  ud_w <- akde_wet[[as.character(id)]]
  ud_d <- akde_dry[[as.character(id)]]
  if (is.null(ud_w) || is.null(ud_d)) return(NULL)
  ci_w <- summary(ud_w)$CI
  ci_d <- summary(ud_d)$CI
  sig  <- ci_w[1] > ci_d[3] | ci_d[1] > ci_w[3]
  sub     <- prep_sub(df_clean, id, max_pts = Inf)
  wet_pts <- sub[is_wet(sub$timestamp), ]
  dry_pts <- sub[!is_wet(sub$timestamp), ]
  cx_w <- mean(wet_pts$x_utm); cy_w <- mean(wet_pts$y_utm)
  cx_d <- mean(dry_pts$x_utm); cy_d <- mean(dry_pts$y_utm)
  shift_km <- sqrt((cx_w - cx_d)^2 + (cy_w - cy_d)^2) / 1000
  data.frame(
    ID = id,
    wet_low = ci_w[1], wet_est = ci_w[2], wet_high = ci_w[3],
    dry_low = ci_d[1], dry_est = ci_d[2], dry_high = ci_d[3],
    ratio_wet_dry     = round(ci_w[2] / ci_d[2], 2),
    CI_non_overlap    = sig,
    centroid_shift_km = round(shift_km, 2)
  )
}))

if (!is.null(seasonal_df) && nrow(seasonal_df) > 0) {
  write.csv(seasonal_df, file.path(CT, "Seasonal_Comparison.csv"), row.names = FALSE)

  if (nrow(seasonal_df) >= 5) {
    wt <- wilcox.test(seasonal_df$wet_est, seasonal_df$dry_est,
                      paired = TRUE, alternative = "two.sided")
    eff <- tryCatch(
      wilcox_effsize(
        data.frame(area = c(seasonal_df$wet_est, seasonal_df$dry_est),
                   season = rep(c("wet","dry"), each = nrow(seasonal_df)),
                   id = rep(seasonal_df$ID, 2)),
        area ~ season, paired = TRUE
      )$effsize, error = function(e) NA)
    sink(file.path(CT, "Seasonal_WilcoxTest.txt"))
    cat("=== Wilcoxon signed-rank test: wet vs dry AKDE area ===\n")
    cat(sprintf("W = %.0f,  p = %.4f\n", wt$statistic, wt$p.value))
    cat(sprintf("Effect size r = %.3f\n", eff))
    cat(sprintf("Mean wet area: %.1f km²\n", mean(seasonal_df$wet_est, na.rm = TRUE)))
    cat(sprintf("Mean dry area: %.1f km²\n", mean(seasonal_df$dry_est, na.rm = TRUE)))
    cat(sprintf("Elephants with non-overlapping CIs: %d of %d\n",
                sum(seasonal_df$CI_non_overlap, na.rm = TRUE), nrow(seasonal_df)))
    sink()
    cat("Wilcoxon test saved → comparison/tables/\n")
  }

  seas_long <- seasonal_df %>%
    select(ID, wet_low, wet_est, wet_high, dry_low, dry_est, dry_high) %>%
    pivot_longer(cols = -ID, names_to = c("season", ".value"), names_sep = "_")

  p_seas <- ggplot(seas_long, aes(x = factor(ID), y = est, colour = season)) +
    geom_point(position = position_dodge(0.5), size = 3) +
    geom_errorbar(aes(ymin = low, ymax = high),
                  position = position_dodge(0.5), width = 0.3, linewidth = 0.8) +
    scale_colour_manual(values = c(wet = "#1D9E75", dry = "#D85A30"),
                        labels = c(wet = "Wet (Nov-Apr)", dry = "Dry (May-Oct)")) +
    labs(title = "Seasonal AKDE home range areas with 95% confidence intervals",
         x = "Elephant ID", y = expression("95% AKDE area (km"^2*")"), colour = "Season") +
    theme_bw(base_size = 12) + theme(legend.position = "top")
  ggsave(file.path(CG, "Seasonal_Comparison_Plot.png"), p_seas,
         width = 10, height = 5, dpi = PNG_DPI)

  p_shift <- ggplot(seasonal_df, aes(x = reorder(factor(ID), -centroid_shift_km),
                                     y = centroid_shift_km)) +
    geom_col(fill = "#7F77DD", alpha = 0.8, width = 0.6) +
    geom_text(aes(label = sprintf("%.1f", centroid_shift_km)), vjust = -0.4, size = 3.5) +
    labs(title = "Seasonal centroid shift (wet to dry)",
         x = "Elephant ID", y = "Centroid shift (km)") +
    theme_bw(base_size = 12)
  ggsave(file.path(CG, "Seasonal_CentroidShift.png"), p_shift,
         width = 8, height = 5, dpi = PNG_DPI)
  cat("Seasonal plots saved → comparison/graphs/\n")
}


# ============================================================
# SECTION I — METHOD COMPARISON
# ============================================================

cat("\n=== Method comparison ===\n")

if (nrow(model_comparison) > 0) {

  mc_wide <- model_comparison %>%
    filter(Method %in% c("AKDE", "dBBMM", "aLOCOH", "TGD", "ETD")) %>%
    group_by(ID, Method) %>%
    summarise(Area_km2 = mean(Area_km2, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(id_cols = ID, names_from = Method, values_from = Area_km2)

  if ("aLOCOH" %in% names(mc_wide) && "AKDE" %in% names(mc_wide))
    mc_wide$ratio_AKDE_aLOCOH <- round(mc_wide$AKDE / mc_wide$aLOCOH, 2)
  if ("dBBMM"  %in% names(mc_wide) && "AKDE" %in% names(mc_wide))
    mc_wide$ratio_AKDE_dBBMM  <- round(mc_wide$AKDE / mc_wide$dBBMM,  2)
  if ("TGD"    %in% names(mc_wide) && "AKDE" %in% names(mc_wide))
    mc_wide$ratio_AKDE_TGD    <- round(mc_wide$AKDE / mc_wide$TGD,    2)
  if ("ETD"    %in% names(mc_wide) && "AKDE" %in% names(mc_wide))
    mc_wide$ratio_AKDE_ETD    <- round(mc_wide$AKDE / mc_wide$ETD,    2)

  tryCatch(write.csv(mc_wide, file.path(CT, "MethodComparison_Wide.csv"), row.names = FALSE),
           error = function(e) NULL)

  sink(file.path(CT, "MethodCorrelations.txt"))
  cat("=== Spearman correlation between home range estimators ===\n\n")
  # Full pairwise matrix (all 10 combinations of the 5 methods), matching the
  # complete 5x5 enumeration already used for Bhattacharyya's Affinity below —
  # previously this only tested 6 of the 10 possible pairs (every pair
  # involving TGD except AKDE-TGD, and every pair involving ETD except
  # AKDE-ETD and dBBMM-ETD, were silently never computed).
  method_pairs <- combn(c("AKDE", "dBBMM", "aLOCOH", "TGD", "ETD"), 2, simplify = FALSE)
  for (pr in method_pairs) {
    if (!all(pr %in% names(mc_wide))) next
    x <- mc_wide[[pr[1]]];  y <- mc_wide[[pr[2]]]
    complete <- complete.cases(x, y)
    if (sum(complete) < 4) {
      cat(sprintf("%s vs %s: insufficient complete pairs\n", pr[1], pr[2])); next
    }
    ct <- cor.test(x[complete], y[complete], method = "spearman")
    cat(sprintf("%s vs %s:  rho = %.3f,  p = %.4f,  n = %d\n",
                pr[1], pr[2], ct$estimate, ct$p.value, sum(complete)))
  }
  sink()
  cat("Method correlations saved → comparison/tables/\n")

  if ("AKDE" %in% names(mc_wide) && "dBBMM" %in% names(mc_wide)) {
    ba_df <- mc_wide %>%
      filter(!is.na(AKDE) & !is.na(dBBMM)) %>%
      mutate(mean_AD = (AKDE + dBBMM) / 2, diff_AD = AKDE - dBBMM)
    if (nrow(ba_df) >= 3) {
      bias   <- mean(ba_df$diff_AD)
      loa_up <- bias + 1.96 * sd(ba_df$diff_AD)
      loa_lo <- bias - 1.96 * sd(ba_df$diff_AD)
      p_ba <- ggplot(ba_df, aes(x = mean_AD, y = diff_AD, label = ID)) +
        geom_point(colour = "#1D9E75", size = 3) +
        ggrepel::geom_text_repel(size = 3, colour = "grey40") +
        geom_hline(yintercept = bias,   colour = "#D85A30", linewidth = 1) +
        geom_hline(yintercept = loa_up, colour = "grey50",  linewidth = 0.7, linetype = "dashed") +
        geom_hline(yintercept = loa_lo, colour = "grey50",  linewidth = 0.7, linetype = "dashed") +
        annotate("text", x = max(ba_df$mean_AD), y = bias,
                 label = sprintf("Bias = %.1f km²", bias),
                 hjust = 1, vjust = -0.5, size = 3.5, colour = "#D85A30") +
        annotate("text", x = max(ba_df$mean_AD), y = loa_up,
                 label = sprintf("+1.96 SD = %.1f km²", loa_up),
                 hjust = 1, vjust = -0.4, size = 3, colour = "grey40") +
        annotate("text", x = max(ba_df$mean_AD), y = loa_lo,
                 label = sprintf("-1.96 SD = %.1f km²", loa_lo),
                 hjust = 1, vjust = 1.2, size = 3, colour = "grey40") +
        labs(title = "Bland-Altman plot: AKDE vs dBBMM agreement",
             x = "Mean of AKDE and dBBMM (km²)", y = "AKDE minus dBBMM (km²)") +
        theme_bw(base_size = 12)
      ggsave(file.path(CG, "BlandAltman_AKDE_dBBMM.png"), p_ba,
             width = 8, height = 6, dpi = PNG_DPI)
      cat("Bland-Altman plot saved → comparison/graphs/\n")
    }
  }

  mc_plot <- model_comparison %>%
    filter(Method %in% c("AKDE", "dBBMM", "aLOCOH", "TGD", "ETD"))
  method_colors <- c(AKDE="#1D9E75", dBBMM="#D85A30", aLOCOH="#7F77DD",
                     TGD="#BA7517", ETD="#2196F3")
  p_mc <- ggplot(mc_plot, aes(x = factor(ID), y = Area_km2, fill = Method)) +
    geom_col(position = "dodge", alpha = 0.85, width = 0.7) +
    scale_fill_manual(values = method_colors[names(method_colors) %in% mc_plot$Method]) +
    labs(title = "Home range area by method and individual",
         x = "Elephant ID", y = expression("Area (km"^2*")"), fill = "Method") +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")
  ggsave(file.path(CG, "MethodComparison_Plot.png"), p_mc,
         width = 12, height = 6, dpi = PNG_DPI)
  cat("Method comparison plot saved → comparison/graphs/\n")

  # ── Bhattacharyya's Affinity (BA) — spatial concordance between methods ────
  # BA ∈ [0,1]: 0 = no overlap, 1 = perfect overlap.
  # Computed pairwise from the 95% area shapefiles of each method per elephant.
  ba_methods <- c("AKDE", "dBBMM", "aLOCOH", "TGD", "ETD")
  shp_roots  <- list(
    AKDE   = file.path(AS),
    dBBMM  = file.path(DS),
    aLOCOH = file.path(LS),
    TGD    = file.path(TS),
    ETD    = file.path(ES)
  )
  shp_stems  <- list(
    AKDE   = "AKDE_95_",
    dBBMM  = "dBBMM_95_",
    aLOCOH = "aLOCOH_95_",
    TGD    = "TGD_",
    ETD    = "ETD_95_"
  )

  # Load 95% isopleths as sf objects (one per method per elephant)
  load_95sf <- function(method, id) {
    p <- file.path(shp_roots[[method]], paste0(shp_stems[[method]], id, ".shp"))
    if (!file.exists(p)) return(NULL)
    tryCatch(sf::st_read(p, quiet = TRUE), error = function(e) NULL)
  }

  # Bhattacharyya's Affinity via intersection-over-union of polygon areas
  ba_pair <- function(sf1, sf2) {
    if (is.null(sf1) || is.null(sf2)) return(NA_real_)
    tryCatch({
      inter  <- suppressWarnings(sf::st_intersection(sf::st_union(sf1), sf::st_union(sf2)))
      a_i    <- as.numeric(sf::st_area(inter))
      a1     <- as.numeric(sf::st_area(sf::st_union(sf1)))
      a2     <- as.numeric(sf::st_area(sf::st_union(sf2)))
      if (a1 <= 0 || a2 <= 0) return(NA_real_)
      # Bhattacharyya's Affinity approximation: sqrt(intersection/A1 * intersection/A2)
      sqrt((sum(a_i) / a1) * (sum(a_i) / a2))
    }, error = function(e) NA_real_)
  }

  all_ids_ba  <- unique(df_clean$elephant_id)
  ba_rows <- do.call(rbind, lapply(all_ids_ba, function(bid) {
    do.call(rbind, lapply(seq_along(ba_methods), function(i) {
      do.call(rbind, lapply(seq_along(ba_methods), function(j) {
        if (i >= j) return(NULL)
        m1 <- ba_methods[i]; m2 <- ba_methods[j]
        sf1 <- load_95sf(m1, bid); sf2 <- load_95sf(m2, bid)
        ba  <- ba_pair(sf1, sf2)
        data.frame(ID = bid, Method1 = m1, Method2 = m2, BA = ba,
                   stringsAsFactors = FALSE)
      }))
    }))
  }))

  if (!is.null(ba_rows) && nrow(ba_rows) > 0) {
    write.csv(ba_rows, file.path(CT, "Method_BA_Overlap.csv"), row.names = FALSE)

    # Mean BA across individuals, one cell per method pair
    ba_mean <- aggregate(BA ~ Method1 + Method2, data = ba_rows, FUN = mean, na.rm = TRUE)
    ba_sym  <- rbind(ba_mean,
                     data.frame(Method1 = ba_mean$Method2,
                                Method2 = ba_mean$Method1, BA = ba_mean$BA))
    ba_sym  <- rbind(ba_sym,
                     data.frame(Method1 = ba_methods, Method2 = ba_methods, BA = 1))
    ba_sym$Method1 <- factor(ba_sym$Method1, levels = ba_methods)
    ba_sym$Method2 <- factor(ba_sym$Method2, levels = rev(ba_methods))

    p_ba_heat <- ggplot(ba_sym, aes(x = Method1, y = Method2, fill = BA)) +
      geom_tile(colour = "white", linewidth = 0.5) +
      geom_text(aes(label = ifelse(is.na(BA), "NA", sprintf("%.2f", BA))),
                size = 3.5) +
      scale_fill_gradient2(low = "#D85A30", mid = "#FFFACD", high = "#1D9E75",
                           midpoint = 0.5, na.value = "grey80",
                           name = "BA overlap") +
      labs(title = "Bhattacharyya's Affinity — Method Spatial Concordance",
           x = NULL, y = NULL) +
      theme_bw(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(file.path(CG, "Method_BA_Concordance.png"), p_ba_heat,
           width = 7, height = 6, dpi = PNG_DPI)
    cat("BA concordance heatmap saved → comparison/graphs/\n")
  }
}


# ============================================================
# SECTION J — INTER-INDIVIDUAL SPACE OVERLAP
# ============================================================

cat("\n=== Inter-individual space overlap ===\n")

ids_with_akde <- names(akde_results)
n_ids         <- length(ids_with_akde)

if (n_ids >= 2) {

  udoi_mat <- matrix(NA_real_, n_ids, n_ids,
                     dimnames = list(ids_with_akde, ids_with_akde))
  for (i in seq_along(ids_with_akde)) {
    udoi_mat[i, i] <- 1
    if (i == n_ids) next
    for (j in (i + 1):n_ids) {
      ov <- tryCatch(
        ctmm::overlap(list(akde_results[[ids_with_akde[i]]],
                           akde_results[[ids_with_akde[j]]]))$CI[, "est"],
        error = function(e) NA_real_
      )
      udoi_mat[i, j] <- ov;  udoi_mat[j, i] <- ov
    }
  }
  write.csv(as.data.frame(udoi_mat), file.path(CT, "Pairwise_UDOI_Matrix.csv"))

  udoi_melt <- melt(udoi_mat, varnames = c("Elephant_A", "Elephant_B"), value.name = "UDOI")
  p_udoi <- ggplot(udoi_melt, aes(x = factor(Elephant_A), y = factor(Elephant_B), fill = UDOI)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = ifelse(!is.na(UDOI), sprintf("%.2f", UDOI), "NA")), size = 3) +
    scale_fill_viridis_c(name = "UDOI", limits = c(0, 1), na.value = "grey80") +
    labs(title = "Pairwise home range overlap (UDOI) — all elephants",
         x = "Elephant ID", y = "Elephant ID") +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(CG, "Pairwise_UDOI_Heatmap.png"), p_udoi,
         width = 8, height = 7, dpi = PNG_DPI)
  cat("UDOI heatmap saved → comparison/graphs/\n")

  udoi_vals <- udoi_mat[upper.tri(udoi_mat)]
  udoi_vals <- udoi_vals[!is.na(udoi_vals)]
  sink(file.path(CT, "UDOI_Summary.txt"))
  cat("=== Pairwise UDOI summary ===\n")
  cat(sprintf("Mean UDOI:   %.3f\n", mean(udoi_vals)))
  cat(sprintf("Median UDOI: %.3f\n", median(udoi_vals)))
  cat(sprintf("Range:       %.3f - %.3f\n", min(udoi_vals), max(udoi_vals)))
  cat(sprintf("Pairs with UDOI > 0.5: %d of %d\n", sum(udoi_vals > 0.5), length(udoi_vals)))
  sink()

  cat("Computing population-level AKDE...\n")
  tryCatch({
    # Step 1: build telemetry for all individuals
    pop_tels_raw <- lapply(ids_with_akde, function(i)
      make_tel(prep_sub(df_clean, as.integer(i))))
    names(pop_tels_raw) <- ids_with_akde

    # Step 2: project ALL telemetry to a single common frame.
    # Use the first individual's projection as reference; ctmm::median on a list
    # is not universally supported, so we fall back to individual 1's CRS.
    pop_proj <- tryCatch(
      ctmm::projection(pop_tels_raw[[1]]),
      error = function(e) { cat("  projection extraction failed:", e$message, "\n"); NULL })
    cat(sprintf("  Common projection: %s\n",
                if (!is.null(pop_proj)) substr(pop_proj, 1, 80) else "NULL"))

    pop_tels <- lapply(seq_along(pop_tels_raw), function(.k) {
      tel <- pop_tels_raw[[.k]]
      if (!is.null(pop_proj))
        tryCatch({ ctmm::projection(tel) <- pop_proj; tel },
                 error = function(e) {
                   cat(sprintf("  reprojection failed for %s: %s\n",
                               ids_with_akde[.k], e$message))
                   tel
                 })
      else tel
    })

    # Verify all projections are identical after alignment
    projs_after <- sapply(pop_tels, function(tel)
      tryCatch(ctmm::projection(tel), error = function(e) "ERR"))
    cat(sprintf("  Unique projections after alignment: %d\n", length(unique(projs_after))))

    # Step 3: fit models on the projected telemetry (ctmm.fit — no model selection)
    cat("  Fitting individual models...\n")
    pop_fits <- lapply(seq_along(ids_with_akde), function(.k) {
      .i <- ids_with_akde[.k]
      .g <- tryCatch(
        ctmm.guess(pop_tels[[.k]], interactive = FALSE),
        error = function(e) { cat(sprintf("    %s guess failed: %s\n", .i, e$message)); NULL })
      if (is.null(.g)) return(NULL)
      .f <- tryCatch(
        ctmm.fit(pop_tels[[.k]], .g),
        error = function(e) { cat(sprintf("    %s fit failed: %s\n", .i, e$message)); NULL })
      if (!is.null(.f)) cat(sprintf("    %s fit OK\n", .i))
      .f
    })
    cat(sprintf("  Valid fits: %d / %d\n",
                sum(!sapply(pop_fits, is.null)), length(ids_with_akde)))

    # Step 4: individual AKDEs on projected telemetry
    cat("  Computing individual AKDEs...\n")
    pop_uds <- lapply(seq_along(ids_with_akde), function(.k) {
      .i  <- ids_with_akde[.k]
      fit <- pop_fits[[.k]]
      if (is.null(fit)) return(NULL)
      ud <- tryCatch(
        akde(pop_tels[[.k]], fit),
        error = function(e) { cat(sprintf("    %s akde failed: %s\n", .i, e$message)); NULL })
      if (!is.null(ud)) cat(sprintf("    %s akde OK\n", .i))
      ud
    })
    cat(sprintf("  Valid UDs: %d / %d\n",
                sum(!sapply(pop_uds, is.null)), length(ids_with_akde)))

    valid <- !sapply(pop_fits, is.null) & !sapply(pop_uds, is.null)
    if (sum(valid) < 2)
      stop("fewer than 2 valid models/UDs for population AKDE")

    # Ensure consistent names so pkde internal checks get length-1 results
    names(pop_tels) <- ids_with_akde
    names(pop_fits) <- ids_with_akde
    names(pop_uds)  <- ids_with_akde

    # Correct ctmm >= 0.7 signature: pkde(data, UD, ...) — no CTMM argument.
    # The individual UDs embed the CTMM model internally.
    pop_ud <- tryCatch(
      pkde(pop_tels[valid], pop_uds[valid]),
      error = function(e) { cat(sprintf("  pkde failed: %s\n", e$message)); NULL })

    if (is.null(pop_ud)) {
      # Fallback: stack individual AKDE rasters and take the mean density
      cat("  Falling back to mean-raster population UD...\n")
      pop_ras_list <- Filter(Negate(is.null), lapply(pop_uds[valid], function(ud)
        tryCatch(raster::raster(ud, DF = "PDF"), error = function(e) NULL)))
      if (length(pop_ras_list) >= 2) {
        pop_stk <- tryCatch(raster::stack(pop_ras_list), error = function(e) NULL)
        pop_ud  <- if (!is.null(pop_stk))
          tryCatch(raster::calc(pop_stk, fun = mean, na.rm = TRUE), error = function(e) NULL)
      }
    }

    if (is.null(pop_ud)) stop("all population AKDE methods failed")

    # Plot — handle both ctmm UD and plain raster objects
    pop_area <- tryCatch(summary(pop_ud)$CI[2], error = function(e) NA_real_)
    save_png({
      if (inherits(pop_ud, "UD")) {
        plot(pop_ud, main = "Population-level AKDE — all elephants",
             level.UD = c(0.50, 0.95), level = 0.95)
        if (!is.na(pop_area))
          title(sub = sprintf("Population 95%% area: %.0f km²", pop_area), cex.sub = 0.85)
      } else {
        raster::plot(pop_ud, col = rev(heat.colors(100)),
                     main = "Population-level AKDE (mean density) — all elephants",
                     xlab = "Easting (m)", ylab = "Northing (m)")
      }
    }, file.path(CG, "Population_AKDE.png"), width = 9, height = 8)
    cat("Population AKDE saved → comparison/graphs/\n")
  }, error = function(e) cat("Population AKDE failed:", e$message, "\n"))
}


# ============================================================
# SECTION K — MOVEMENT BEHAVIOUR
# ============================================================

cat("\n=== Movement behaviour analysis ===\n")

step_summary <- do.call(rbind, lapply(unique(df_clean$elephant_id), function(id) {
  sub <- df_clean[df_clean$elephant_id == id, ]
  sub <- sub[order(sub$timestamp), ]
  dx <- diff(sub$x_utm);  dy <- diff(sub$y_utm)
  dt <- as.numeric(diff(sub$timestamp), units = "secs") / 3600
  sl_km    <- sqrt(dx^2 + dy^2) / 1000
  spd      <- sl_km / pmax(dt, 0.01)
  headings <- atan2(dy, dx)
  ta_rad   <- diff(headings)
  ta_rad   <- ((ta_rad + pi) %% (2 * pi)) - pi

  save_png({
    par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))
    hist(sl_km, breaks = 50, col = "#1D9E75", border = "white",
         main = paste("Step lengths\nElephant", id), xlab = "Step length (km)", freq = FALSE)
    lines(density(sl_km, na.rm = TRUE), col = "#D85A30", lwd = 2)
    hist(ta_rad * 180 / pi, breaks = 36, col = "#7F77DD", border = "white",
         main = "Turning angles", xlab = "Turning angle (deg)", freq = FALSE)
    abline(v = 0, col = "#D85A30", lwd = 2, lty = 2)
    hist(spd[spd <= 10], breaks = 40, col = "#BA7517", border = "white",
         main = "Speed distribution\n(<= 10 km/h shown)", xlab = "Speed (km/h)", freq = FALSE)
    abline(v = 4, col = "#D85A30", lwd = 2, lty = 2)
    legend("topright", legend = "4 km/h ceiling", lty = 2, col = "#D85A30", bty = "n", cex = 0.85)
  }, file.path(MG, paste0("MovementBehaviour_", id, ".png")), width = 12, height = 5)

  data.frame(
    ID = id, n_steps = length(sl_km),
    median_sl_km     = round(median(sl_km, na.rm = TRUE), 3),
    mean_sl_km       = round(mean(sl_km, na.rm = TRUE), 3),
    sd_sl_km         = round(sd(sl_km, na.rm = TRUE), 3),
    max_sl_km        = round(max(sl_km, na.rm = TRUE), 2),
    median_speed_kmh = round(median(spd, na.rm = TRUE), 3),
    mean_abs_ta_deg  = round(mean(abs(ta_rad), na.rm = TRUE) * 180 / pi, 1),
    pct_backtrack    = round(mean(abs(ta_rad) > pi / 2, na.rm = TRUE) * 100, 1)
  )
}))

write.csv(step_summary, file.path(MT, "StepLength_Summary.csv"), row.names = FALSE)

# Local time (Africa/Harare, UTC+2) and step speed over steps <= 4 h, so that
# 6-hourly and gap-spanning steps do not inflate the profile (revision, 2026-10-30)
df_clean$hour <- as.integer(format(df_clean$timestamp, "%H", tz = "Africa/Harare"))
df_clean <- df_clean %>%
  arrange(elephant_id, timestamp) %>%
  group_by(elephant_id) %>%
  mutate(sl_km = c(NA_real_, sqrt(diff(x_utm)^2 + diff(y_utm)^2) / 1000),
         dt_h  = c(NA_real_, as.numeric(diff(timestamp), units = "secs") / 3600)) %>%
  ungroup()

diurnal <- df_clean %>%
  filter(!is.na(sl_km), dt_h > 0, dt_h <= 4) %>%
  mutate(spd = sl_km / dt_h) %>%
  group_by(elephant_id, hour) %>%
  summarise(mean_sl = mean(spd, na.rm = TRUE),
            se_sl   = sd(spd, na.rm = TRUE) / sqrt(n()), .groups = "drop")

p_diurnal <- ggplot(diurnal, aes(x = hour, y = mean_sl,
                                 colour = factor(elephant_id), group = factor(elephant_id))) +
  geom_ribbon(aes(ymin = mean_sl - se_sl, ymax = mean_sl + se_sl,
                  fill = factor(elephant_id)), alpha = 0.1, colour = NA) +
  geom_line(linewidth = 0.9) + geom_point(size = 1.5) +
  scale_x_continuous(breaks = seq(0, 23, 3)) +
  annotate("rect", xmin = 18, xmax = 23, ymin = -Inf, ymax = Inf, alpha = 0.06, fill = "navy") +
  annotate("rect", xmin = 0,  xmax = 6,  ymin = -Inf, ymax = Inf, alpha = 0.06, fill = "navy") +
  labs(title = "Diurnal movement pattern by elephant (shading = night)",
       x = "Hour of day (local time, UTC+2)", y = "Mean step speed (km/h +/- SE)",
       colour = "Elephant ID", fill = "Elephant ID") +
  theme_bw(base_size = 12) + theme(legend.position = "right")
ggsave(file.path(MG, "Diurnal_Pattern.png"), p_diurnal, width = 12, height = 5, dpi = PNG_DPI)

cat("Computing area accumulation curves...\n")
accum_list <- lapply(ids_with_akde, function(id) {
  sub <- prep_sub(df_clean, id)
  if (nrow(sub) < 40) return(NULL)
  tel <- make_tel(sub)
  fit <- tryCatch(ctmm.select(tel, ctmm.guess(tel, interactive = FALSE),
                              cores = N_CORES, verbose = TRUE)[[1]],
                  error = function(e) ctmm.guess(tel, interactive = FALSE))
  ns    <- unique(round(seq(30, length(tel), length.out = 12)))
  ns    <- ns[ns >= 20]
  areas <- sapply(ns, function(n) {
    tryCatch({ idx <- sort(sample(length(tel), n)); summary(akde(tel[idx], fit))$CI[2] },
             error = function(e) NA_real_)
  })
  data.frame(ID = id, n_fixes = ns, area_km2 = areas)
})
accum_df <- do.call(rbind, Filter(Negate(is.null), accum_list))
if (!is.null(accum_df) && nrow(accum_df) > 0) {
  p_accum <- ggplot(accum_df, aes(x = n_fixes, y = area_km2,
                                  colour = factor(ID), group = factor(ID))) +
    geom_line(linewidth = 0.9) + geom_point(size = 2) +
    labs(title = "Area accumulation curve — sampling adequacy",
         x = "Number of GPS fixes used", y = expression("AKDE 95% area (km"^2*")"),
         colour = "Elephant ID") +
    theme_bw(base_size = 12) + theme(legend.position = "right")
  ggsave(file.path(MG, "AreaAccumulation_Curves.png"), p_accum,
         width = 10, height = 6, dpi = PNG_DPI)
  cat("Area accumulation curves saved → movement/graphs/\n")
}

cat("Computing first-passage time (FPT)...\n")
for (id in ids_with_akde) {
  out_fpt <- file.path(MG, paste0("FPT_", id, ".png"))
  tryCatch({
    # Use the full speed-filtered dataset (no thinning) for maximum temporal density
    sub_fpt <- as.data.frame(df_clean[df_clean$elephant_id == as.integer(id), ])
    sub_fpt <- sub_fpt[order(sub_fpt$timestamp), ]
    sub_fpt <- sub_fpt[!duplicated(sub_fpt$timestamp), ]
    if (nrow(sub_fpt) < 30) {
      cat(sprintf("  FPT skipped for %s: too few points\n", id)); next
    }
    radii <- seq(500, 50000, by = 1000)   # 500 m to 50 km in 1-km steps
    vl1   <- .fpt_varlog(sub_fpt, radii)
    valid_idx <- which(!is.na(vl1))
    if (length(valid_idx) < 2) {
      cat(sprintf("  FPT skipped for %s: no valid variance values\n", id)); next
    }
    radii_p <- radii[valid_idx] / 1000   # km for plotting
    vl1     <- vl1[valid_idx]
    r_peak  <- radii_p[which.max(vl1)]
    save_png({
      plot(radii_p, vl1, type = "b", pch = 19, col = "#1D9E75",
           main = paste("First-passage time — Elephant", id),
           xlab = "Radius (km)", ylab = "Var(log FPT)")
      abline(v = r_peak, col = "#D85A30", lty = 2, lwd = 1.5)
      text(r_peak, max(vl1) * 0.9,
           labels = sprintf("Peak: %.1f km", r_peak),
           pos = 4, cex = 0.85, col = "#D85A30")
    }, out_fpt)
    cat(sprintf("  FPT done for Elephant %s — peak: %.1f km\n", id, r_peak))
  }, error = function(e) cat(sprintf("  FPT failed for %s: %s\n", id, e$message)))
}


# ============================================================
# SECTION L — FINAL SUMMARY OUTPUTS
# ============================================================

cat("\n=== Saving final summary outputs ===\n")

if (nrow(model_comparison) > 0) {
  tryCatch(
    write.csv(model_comparison, file.path(CT, "Model_Comparison_Stats.csv"), row.names = FALSE),
    error = function(e) {
      alt <- file.path(CT, paste0("Model_Comparison_Stats_", format(Sys.time(), "%H%M%S"), ".csv"))
      write.csv(model_comparison, alt, row.names = FALSE)
      cat("Saved as", basename(alt), "(original was locked)\n")
    }
  )
}

sink(file.path(CT, "SessionInfo.txt"))
cat("Analysis completed:", format(Sys.time()), "\n\n")
print(sessionInfo())
sink()

cat("\n=== All analyses complete ===\n")
cat(sprintf("Results saved to: %s\n", FINAL_DIR))
cat("  AKDE/       graphs/  tables/  spatial/\n")
cat("  dBBMM/      graphs/  spatial/\n")
cat("  aLOCOH/     graphs/  spatial/\n")
cat("  TGD/        graphs/  spatial/\n")
cat("  ETD/        graphs/  tables/  spatial/\n")
cat("  comparison/ graphs/  tables/\n")
cat("  data_quality/ graphs/ tables/\n")
cat("  movement/   graphs/  tables/\n")
