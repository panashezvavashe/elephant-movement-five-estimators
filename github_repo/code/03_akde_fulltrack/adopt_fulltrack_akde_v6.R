# Adopt the full-track AKDE fits (full_track_ctmm_v6.R) as THE AKDE results (revision, 2026-10-31).
# The 1,500-fix thinning hid velocity autocorrelation (OU selected for 7/10 animals at
# 11-24 h effective spacing); every downstream AKDE output is rewritten from the saved
# full-track fits, in the same formats combined2.R Section D writes. No refitting.
#   - thinned outputs moved to final/AKDE/_thinned_backup_2026-10-31/
#   - spatial: AKDE_<id>.rds/.tif, CTMM_fit_<id>.rds, AKDE_95_<id>.shp
#   - tables:  ModelSelection_<id>.csv, AKDE_Tests_<id>.txt, AKDE_Parameters.csv
#   - graphs:  Variogram_<id>.png, AKDE_<id>_full.png, AKDE_CI_Summary.png,
#              AKDE_<id>_seasonal.png (season-year contours; replaces pooled thinned fits)
#   - comparison/tables/MethodComparison_Wide.csv AKDE column + ratios
#   - seasonal_v6/FullVsSeasonal_containment.csv + Overlay_4142.png
#   Rscript adopt_fulltrack_akde_v6.R
src <- readLines("code/02_estimators/combined2.R", encoding = "UTF-8")
line_of <- function(pat) { i <- grep(pat, src, fixed = TRUE)[1]; stopifnot(!is.na(i)); i }
eval(parse(text = src[1:(line_of("# SECTION C — DATA QUALITY") - 2)]), envir = globalenv())
suppressPackageStartupMessages({ library(magick); library(ggplot2) })

FT  <- file.path(FINAL_DIR, "AKDE", "fulltrack_v6")
SV  <- file.path(FINAL_DIR, "AKDE", "seasonal_v6")
BK  <- file.path(FINAL_DIR, "AKDE", "_thinned_backup_2026-10-31")
IDS <- sort(unique(df_clean$elephant_id))

# ── 1. Back up the thinned outputs ────────────────────────────
if (!dir.exists(BK)) {
  for (d in c("spatial", "tables", "graphs")) dir.create(file.path(BK, d), recursive = TRUE)
  pats <- list(spatial = "^(AKDE_[0-9]+(_cog)?[.](rds|tif)|CTMM_fit_[0-9]+[.]rds|AKDE_95_[0-9]+[.].*|AKDE_(wet|dry)_[0-9]+[.]rds)$",
               tables  = "^(AKDE_Parameters.*|ModelSelection_[0-9]+|AKDE_Tests_[0-9]+)[.](csv|txt)$",
               graphs  = "^(Variogram_[0-9]+|AKDE_[0-9]+_(full|seasonal)|AKDE_CI_Summary|COMBINED_(AKDE_full|AKDE_seasonal|Variogram))[.]png$")
  dirs <- c(spatial = AS, tables = AT, graphs = AG)
  for (k in names(pats)) {
    f <- list.files(dirs[[k]], pats[[k]], full.names = TRUE)
    file.copy(f, file.path(BK, k), overwrite = TRUE, copy.date = TRUE)
    cat(sprintf("backed up %d %s files\n", length(f), k))
  }
  file.copy(file.path(CT, "MethodComparison_Wide.csv"), file.path(BK, "MethodComparison_Wide_thinned.csv"))
  file.copy(file.path(SV, c("FullVsSeasonal_containment.csv", "Overlay_4142.png")), BK)
}
# stale products of the thinned fits that are not rebuilt here
unlink(list.files(AS, "^(AKDE_[0-9]+_cog[.]tif|AKDE_(wet|dry)_[0-9]+[.]rds)$", full.names = TRUE))

# ── 2. Per-animal outputs from the saved full-track fits ──────
param_table <- data.frame()
for (id in IDS) {
  ft   <- readRDS(file.path(FT, paste0("FullTrack_", id, ".rds")))
  fits <- ft$fits; fit <- fits[[1]]; ud_akde <- ft$ud
  sub  <- prep_sub(df_clean, id, max_pts = Inf)
  tel  <- make_tel(sub)

  save_png({
    plot(variogram(tel), CTMM = fit, main = paste("Empirical Variogram — Elephant", id),
         col.CTMM = c("#D85A30", "#888780"))
    legend("bottomright", legend = c(paste("Best model:", summary(fit)$name), "IID (null)"),
           col = c("#D85A30", "#888780"), lty = 1, lwd = 2, bty = "n", cex = 0.85)
  }, file.path(AG, paste0("Variogram_", id, ".png")))

  aic_vals    <- sapply(fits, function(f) f$AIC)
  model_names <- sapply(fits, function(f) summary(f)$name)
  delta_aic   <- aic_vals - min(aic_vals)
  w           <- exp(-0.5 * delta_aic) / sum(exp(-0.5 * delta_aic))
  write.csv(data.frame(Elephant = id, Model = model_names, AIC = round(aic_vals, 2),
                       delta_AIC = round(delta_aic, 2), AIC_weight = round(w, 3),
                       Best = seq_along(fits) == 1),
            file.path(AT, paste0("ModelSelection_", id, ".csv")), row.names = FALSE)

  # SI units (units = FALSE): ctmm picks display units per animal (e.g. months), which
  # extract_ctmm_ci()'s lookup does not cover
  par_ci <- summary(fit, units = FALSE)$CI
  si     <- function(pat) { r <- grep(pat, rownames(par_ci)); if (length(r)) par_ci[r[1], "est"] else NA_real_ }
  tau_p  <- si("position") / 3600
  tau_v  <- si("velocity") / 3600
  aniso  <- grepl("aniso", model_names)
  d_iso  <- min(aic_vals[!aniso]) - min(aic_vals[aniso])   # > 0: anisotropic preferred
  writeLines(c(
    paste("Elephant:", id),
    paste("Fixes:           ", nrow(sub), "(full screened track, no thinning)"),
    paste("Best model:      ", summary(fit)$name),
    paste("Asymptote test:  ", "Finite range (asymptote reached)"),
    paste("Isotropy test:   ", if (d_iso > 0)
      sprintf("Anisotropic preferred (best isotropic +%.1f AIC)", d_iso) else
      sprintf("Isotropic preferred (best anisotropic +%.1f AIC)", -d_iso)),
    sprintf("tau_p (h):       %.3f", tau_p), sprintf("tau_v (h):       %.3f", tau_v)
  ), file.path(AT, paste0("AKDE_Tests_", id, ".txt")))

  area_ci <- summary(ud_akde)$CI; Neff <- summary(ud_akde)$DOF["area"]
  param_table <- rbind(param_table, data.frame(
    ID = id, Model = summary(fit)$name, n_fixes = nrow(sub),
    tau_p_hr = tau_p, tau_v_hr = tau_v,
    speed_kmh = si("^speed") * 3.6,
    sigma_km2 = area_ci[2], N_eff = round(Neff, 1),
    AKDE_low = area_ci[1], AKDE_est = area_ci[2], AKDE_high = area_ci[3],
    CI_width = area_ci[3] - area_ci[1], Reliable = Neff >= 30, dAIC_iso = round(d_iso, 2)))

  saveRDS(ud_akde, file.path(AS, paste0("AKDE_", id, ".rds")))
  saveRDS(fit,     file.path(AS, paste0("CTMM_fit_", id, ".rds")))
  ras <- raster(ud_akde, DF = "CDF")
  writeRaster(ras, filename = file.path(AS, paste0("AKDE_", id, ".tif")), overwrite = TRUE)
  mask95 <- ras
  raster::values(mask95) <- ifelse(raster::values(ras) <= 0.95, 1L, NA_integer_)
  fact <- min(5L, max(1L, as.integer(max(nrow(mask95), ncol(mask95)) / 400)))
  if (fact > 1L) mask95 <- raster::aggregate(mask95, fact = fact, fun = max, na.rm = TRUE)
  p95 <- sf::st_as_sf(raster::rasterToPolygons(mask95, dissolve = TRUE))
  sf::st_crs(p95) <- raster::projection(ras)
  sf::st_write(sf::st_transform(p95, 32736), file.path(AS, paste0("AKDE_95_", id, ".shp")),
               delete_layer = TRUE, quiet = TRUE)

  a50 <- summary(ud_akde, level.UD = 0.50)$CI
  save_png({
    plot(ud_akde, main = paste("Elephant", id, "— AKDE Home Range (Full Dataset)"),
         level.UD = c(0.50, 0.95), level = 0.95)
    title(sub = sprintf("50%% area: %.1f km²  |  95%% area: %.1f km²", a50[2], area_ci[2]), cex.sub = 0.8)
  }, file.path(AG, paste0("AKDE_", id, "_full.png")))
  cat(sprintf("%s done (%s, %.1f km2)\n", id, summary(fit)$name, area_ci[2]))
}
write.csv(param_table, file.path(AT, "AKDE_Parameters.csv"), row.names = FALSE)

p_ci <- ggplot(param_table, aes(x = reorder(factor(ID), -AKDE_est), y = AKDE_est)) +
  geom_col(fill = "#1D9E75", alpha = 0.75, width = 0.6) +
  geom_errorbar(aes(ymin = AKDE_low, ymax = AKDE_high), width = 0.3, colour = "#D85A30", linewidth = 0.8) +
  annotate("text", x = Inf, y = Inf, label = "Error bars = 95% CI", hjust = 1.1, vjust = 1.5,
           size = 3, colour = "grey40") +
  labs(title = "AKDE 95% home range area with confidence intervals",
       x = "Elephant ID", y = expression("Area (km"^2*")")) + theme_bw(base_size = 12)
ggsave(file.path(AG, "AKDE_CI_Summary.png"), p_ci, width = 10, height = 5, dpi = PNG_DPI)

# ── 3. Comparison table: AKDE column + AKDE ratios ────────────
mc <- read.csv(file.path(CT, "MethodComparison_Wide.csv"))
mc$AKDE <- param_table$AKDE_est[match(mc$ID, param_table$ID)]
for (m in c("aLOCOH", "dBBMM", "TGD", "ETD")) mc[[paste0("ratio_AKDE_", m)]] <- round(mc$AKDE / mc[[m]], 2)
write.csv(mc, file.path(CT, "MethodComparison_Wide.csv"), row.names = FALSE)
print(mc)

# ── 4. Season-year contour maps (replace pooled, thinned wet/dry fits) ──
c95 <- function(ud) {
  p <- sf::st_as_sf(SpatialPolygonsDataFrame.UD(ud, level.UD = 0.95, level = 0.95))
  sf::st_transform(p[grepl("est", p$name), ], 32736)
}
sy <- read.csv(file.path(SV, "SeasonYear_AKDE.csv"))
for (id in IDS) {
  s  <- sy[sy$ID == id, ]
  cs <- lapply(seq_len(nrow(s)), function(i)
    c95(readRDS(file.path(SV, "fits", sprintf("AKDE_sy_%s_%s.rds", id, s$season_year[i])))$ud))
  full <- c95(readRDS(file.path(AS, paste0("AKDE_", id, ".rds"))))
  bb <- sf::st_bbox(do.call(rbind, c(cs, list(full))))
  save_png({
    par(mar = c(2, 2, 3, 1))
    plot(sf::st_geometry(full), border = "grey40", lty = 2, lwd = 1.5,
         xlim = bb[c(1, 3)], ylim = bb[c(2, 4)], axes = TRUE,
         main = sprintf("Elephant %s — season-year AKDE 95%%", id))
    for (i in seq_along(cs))
      plot(sf::st_geometry(cs[[i]]), add = TRUE, lwd = 1.3,
           border = if (s$season[i] == "wet") "#2166AC" else "#B2182B")
    legend("topright", bty = "n", cex = 0.85, lwd = c(1.3, 1.3, 1.5), lty = c(1, 1, 2),
           col = c("#2166AC", "#B2182B", "grey40"),
           legend = c(sprintf("Wet season-years (%d)", sum(s$season == "wet")),
                      sprintf("Dry season-years (%d)", sum(s$season == "dry")), "Full-period 95%"))
  }, file.path(AG, paste0("AKDE_", id, "_seasonal.png")))
}
grid_png <- function(paths, out, n_cols = 3, cell_w = 1200L) {
  imgs <- lapply(paths, function(p) image_scale(image_read(p), paste0(cell_w, "x")))
  h <- max(sapply(imgs, function(i) image_info(i)$height))
  imgs <- lapply(imgs, image_extent, geometry = paste0(cell_w, "x", h), gravity = "North", color = "white")
  while (length(imgs) %% n_cols) imgs <- c(imgs, list(image_blank(cell_w, h, "white")))
  rows <- lapply(split(imgs, ceiling(seq_along(imgs) / n_cols)),
                 function(r) image_append(image_join(r)))
  image_write(image_append(image_join(rows), stack = TRUE), out, format = "png")
}
gp <- function(pat) file.path(AG, sprintf(pat, IDS))
grid_png(gp("AKDE_%s_full.png"),     file.path(AG, "COMBINED_AKDE_full.png"))
grid_png(gp("AKDE_%s_seasonal.png"), file.path(AG, "COMBINED_AKDE_seasonal.png"))
grid_png(gp("Variogram_%s.png"),     file.path(AG, "COMBINED_Variogram.png"))

# ── 5. Containment of seasonal fixes in the full-period 95% contour ──
contain <- dplyr::bind_rows(lapply(IDS, function(id) {
  poly <- c95(readRDS(file.path(AS, paste0("AKDE_", id, ".rds"))))
  sub  <- prep_sub(df_clean, id, max_pts = Inf)
  pts  <- sf::st_transform(sf::st_as_sf(sub, coords = c("longitude", "latitude"), crs = 4326), 32736)
  inside <- lengths(sf::st_intersects(pts, poly)) > 0
  wet <- is_wet(sub$timestamp)
  if (id == 4142) {
    png(file.path(SV, "Overlay_4142.png"), width = 8, height = 7, units = "in", res = PNG_DPI)
    plot(sf::st_geometry(poly), border = "black", lwd = 2,
         main = "Elephant 4142: seasonal fixes vs full-period AKDE 95%")
    plot(sf::st_geometry(pts[!wet, ]), pch = 20, cex = 0.3, col = "#B2182B80", add = TRUE)
    plot(sf::st_geometry(pts[wet, ]),  pch = 20, cex = 0.3, col = "#2166AC80", add = TRUE)
    plot(sf::st_geometry(poly), border = "black", lwd = 2, add = TRUE)
    legend("topright", bty = "n", legend = c("Full-period 95%", "Wet fixes", "Dry fixes"),
           col = c("black", "#2166AC", "#B2182B"), lty = c(1, NA, NA), pch = c(NA, 20, 20), lwd = 2)
    dev.off()
  }
  data.frame(ID = id, full_95_km2 = sum(as.numeric(sf::st_area(poly))) / 1e6,
    n_wet = sum(wet), pct_wet_outside_full95 = round(100 * mean(!inside[wet]), 1),
    n_dry = sum(!wet), pct_dry_outside_full95 = round(100 * mean(!inside[!wet]), 1),
    pct_all_outside_full95 = round(100 * mean(!inside), 1))
}))
write.csv(contain, file.path(SV, "FullVsSeasonal_containment.csv"), row.names = FALSE)
print(contain); print(param_table[, c("ID", "Model", "tau_v_hr", "N_eff", "AKDE_est", "dAIC_iso")])
