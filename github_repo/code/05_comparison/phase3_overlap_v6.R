# ============================================================
# Phase 3 (revision, 2026-10-29): between-estimator spatial overlap, v6 outputs
#   1. Bhattacharyya coefficient (BC) between the three probabilistic UDs
#      (AKDE, dBBMM, ETD): each UD -> density on a common 250 m UTM 36S grid,
#      renormalised; BC = sum(sqrt(p1 * p2)). Replaces the failed UDOI step.
#   2. 95% polygons, all five estimators: Jaccard index |A∩B| / |A∪B| and
#      centroid distance. TGD and aLoCoH have no UD, so this is the only overlap
#      measure that includes them.
#   3. Probability mass of each UD that falls inside the TGD / aLoCoH polygons.
# BC/Jaccard are per animal; across-animal means get a case-resampling bootstrap CI.
#
# Input:  final/<method>/spatial/   (v6 screening)
# Output: final/comparison/v6/Overlap_*.csv, Overlap_figure.png
#
#   Rscript phase3_overlap_v6.R
# ============================================================

Sys.unsetenv(c("PROJ_LIB", "PROJ_DATA", "GDAL_DATA"))
suppressPackageStartupMessages({
  library(ctmm); library(raster); library(sf); library(dplyr); library(tidyr)
  library(ggplot2); library(cowplot)
})

DATA_DIR <- "."   # run from the repository root
setwd(DATA_DIR)
F   <- file.path(DATA_DIR, "final")
OUT <- file.path(F, "comparison", "v6")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
set.seed(20261029)
N_BOOT  <- 10000L
PNG_DPI <- 300L
RES     <- 250
UTM     <- "+proj=utm +zone=36 +south +datum=WGS84 +units=m +no_defs"
IDS     <- c(3905, 3906, 3907, 4141, 4142, 4501, 4556, 4558, 4559, 4560)
UDS     <- c("AKDE", "dBBMM", "ETD")
METHODS <- c("AKDE", "dBBMM", "ETD", "TGD", "aLOCOH")
LABELS  <- c(AKDE = "AKDE", dBBMM = "dBBMM", ETD = "SAEK", TGD = "TGD", aLOCOH = "aLoCoH")
GROUP   <- c(AKDE = "range", dBBMM = "occurrence", ETD = "occurrence",
             TGD = "occurrence", aLOCOH = "fix geometry")

# ── UD loaders: each returns a density raster (per m²) in its native grid ──
ud_akde <- function(id) {
  ud <- readRDS(file.path(F, "AKDE", "spatial", paste0("AKDE_", id, ".rds")))
  raster(ud, DF = "PDF")                 # ctmm local projection
}
ud_dbbmm <- function(id) {
  # stored as cumulative volume (0-1, outside = 1); cell mass = increment of
  # the sorted volume values
  r <- raster(file.path(F, "dBBMM", "spatial", paste0("dBBMM_", id, ".tif")))
  v <- values(r); o <- order(v)
  p <- numeric(length(v)); p[o] <- diff(c(0, v[o]))
  values(r) <- p / prod(res(r)); r
}
ud_etd <- function(id) {
  r <- raster(file.path(F, "ETD", "spatial", paste0("ETD_", id, "_full.tif")))
  values(r) <- values(r) / 1e6; r          # per km² -> per m²
}

to_template <- function(r, tmpl) {
  x <- if (compareCRS(r, tmpl)) resample(r, tmpl, method = "bilinear")
       else projectRaster(r, to = tmpl, method = "bilinear")
  v <- values(x); v[is.na(v) | v < 0] <- 0
  m <- sum(v) * RES^2                        # mass captured by the template
  list(p = v / sum(v), mass = m)
}

read_poly <- function(m, id) {
  f <- switch(m,
    AKDE   = file.path(F, "AKDE",   "spatial", paste0("AKDE_95_",   id, ".shp")),
    dBBMM  = file.path(F, "dBBMM",  "spatial", paste0("dBBMM_95_",  id, ".shp")),
    ETD    = file.path(F, "ETD",    "spatial", paste0("ETD_95_",    id, ".shp")),
    TGD    = file.path(F, "TGD",    "spatial", paste0("TGD_",       id, ".shp")),
    aLOCOH = file.path(F, "aLOCOH", "spatial", paste0("aLOCOH_95_", id, ".shp")))
  p <- st_read(f, quiet = TRUE)
  if (is.na(st_crs(p))) st_crs(p) <- 32736
  st_union(st_make_valid(st_transform(p, 32736)))
}

bc_rows <- list(); poly_rows <- list(); mass_rows <- list(); qc_rows <- list()
for (id in IDS) {
  cat("Elephant", id, "\n")
  dens <- list(AKDE = ud_akde(id), dBBMM = ud_dbbmm(id), ETD = ud_etd(id))
  ext_utm <- lapply(dens, function(r)
    if (compareCRS(r, CRS(UTM))) extent(r) else extent(projectExtent(r, crs = CRS(UTM))))
  bb <- sapply(ext_utm, function(x) c(x@xmin, x@xmax, x@ymin, x@ymax))
  e <- extent(floor(min(bb[1, ]) / RES) * RES, ceiling(max(bb[2, ]) / RES) * RES,
              floor(min(bb[3, ]) / RES) * RES, ceiling(max(bb[4, ]) / RES) * RES)
  tmpl <- raster(e, res = RES, crs = CRS(UTM))
  P <- lapply(dens, to_template, tmpl = tmpl)
  qc_rows[[length(qc_rows) + 1]] <- data.frame(ID = id, method = UDS,
    mass_on_grid = round(sapply(P, `[[`, "mass"), 4), grid_cells = ncell(tmpl))

  for (pr in combn(UDS, 2, simplify = FALSE)) {
    p1 <- P[[pr[1]]]$p; p2 <- P[[pr[2]]]$p
    bc_rows[[length(bc_rows) + 1]] <- data.frame(ID = id,
      pair = paste(LABELS[pr[1]], "vs", LABELS[pr[2]]),
      type = paste(GROUP[pr[1]], GROUP[pr[2]], sep = "-"),
      BC = sum(sqrt(p1 * p2)),
      VI = sum(pmin(p1, p2)))               # volume of intersection, for reference
  }

  polys <- lapply(setNames(METHODS, METHODS), read_poly, id = id)
  for (pr in combn(METHODS, 2, simplify = FALSE)) {
    a <- polys[[pr[1]]]; b <- polys[[pr[2]]]
    ai <- as.numeric(st_area(st_intersection(a, b))); if (!length(ai)) ai <- 0
    au <- as.numeric(st_area(st_union(a, b)))
    poly_rows[[length(poly_rows) + 1]] <- data.frame(ID = id,
      pair = paste(LABELS[pr[1]], "vs", LABELS[pr[2]]),
      type = paste(GROUP[pr[1]], GROUP[pr[2]], sep = "-"),
      Jaccard = ai / au,
      share_A_in_B = ai / as.numeric(st_area(a)),
      share_B_in_A = ai / as.numeric(st_area(b)),
      centroid_km = as.numeric(st_distance(st_centroid(a), st_centroid(b))) / 1000)
  }

  # UD mass inside the polygon-only estimators
  cells <- st_as_sf(as.data.frame(coordinates(tmpl)), coords = c("x", "y"), crs = 32736)
  for (pm in c("TGD", "aLOCOH")) {
    inside <- lengths(st_intersects(cells, polys[[pm]])) > 0
    for (u in UDS) mass_rows[[length(mass_rows) + 1]] <- data.frame(ID = id,
      polygon = LABELS[pm], UD = LABELS[u], UD_mass_inside = sum(P[[u]]$p[inside]))
  }
}
bc   <- do.call(rbind, bc_rows)
poly <- do.call(rbind, poly_rows)
mass <- do.call(rbind, mass_rows)
qc   <- do.call(rbind, qc_rows)

# ── Across-animal summaries (all / no 4501) ───────────────────
boot_mean <- function(x) quantile(replicate(N_BOOT, mean(sample(x, replace = TRUE))), c(.025, .975))
summ <- function(d, val, by) {
  subsets <- list(all = d, no4501 = d[d$ID != 4501, ])
  do.call(rbind, lapply(names(subsets), function(nm) {
    s <- subsets[[nm]] %>% group_by(across(all_of(by))) %>%
      summarise(n = n(), mean = mean(.data[[val]]),
                ci = list(boot_mean(.data[[val]])),
                median = median(.data[[val]]), min = min(.data[[val]]), max = max(.data[[val]]),
                .groups = "drop")
    s$boot_low  <- sapply(s$ci, `[`, 1)
    s$boot_high <- sapply(s$ci, `[`, 2)
    cbind(subset = nm, as.data.frame(dplyr::select(s, -ci)))
  }))
}
bc_sum   <- summ(bc, "BC", c("pair", "type"))
jac_sum  <- summ(poly, "Jaccard", c("pair", "type"))
cen_sum  <- summ(poly, "centroid_km", c("pair", "type"))
mass_sum <- summ(mass, "UD_mass_inside", c("polygon", "UD"))
rnd <- function(d) { n <- sapply(d, is.numeric) & names(d) != "n"; d[n] <- lapply(d[n], round, 3); d }

write.csv(rnd(bc),   file.path(OUT, "Overlap_BC_per_animal.csv"), row.names = FALSE)
write.csv(rnd(poly), file.path(OUT, "Overlap_polygon_per_animal.csv"), row.names = FALSE)
write.csv(rnd(mass), file.path(OUT, "Overlap_UDmass_in_polygons_per_animal.csv"), row.names = FALSE)
write.csv(qc,        file.path(OUT, "Overlap_grid_QC.csv"), row.names = FALSE)
write.csv(rnd(rbind(cbind(measure = "BC", bc_sum),
                    cbind(measure = "Jaccard", jac_sum),
                    cbind(measure = "centroid_km", cen_sum))),
          file.path(OUT, "Overlap_summary.csv"), row.names = FALSE)
write.csv(rnd(mass_sum), file.path(OUT, "Overlap_UDmass_summary.csv"), row.names = FALSE)

# ── Figure: A Jaccard matrix (lower) + BC (upper); B per-animal BC; C example map ──
lev <- LABELS[METHODS]
js  <- jac_sum[jac_sum$subset == "all", ]
bs  <- bc_sum[bc_sum$subset == "all", ]
split_pair <- function(p) do.call(rbind, strsplit(p, " vs "))
jm <- data.frame(split_pair(js$pair), val = js$mean, what = "Jaccard (95% polygons)")
bm <- data.frame(split_pair(bs$pair), val = bs$mean, what = "BC (UDs)")
names(jm)[1:2] <- names(bm)[1:2] <- c("a", "b")
mat <- rbind(transform(jm, x = a, y = b), transform(bm, x = b, y = a))
mat$x <- factor(mat$x, lev); mat$y <- factor(mat$y, rev(lev))
pA <- ggplot(mat, aes(x, y, fill = val)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", val)), size = 3.4) +
  scale_fill_gradient("Mean overlap", low = "#F7FBFF", high = "#2171B5", limits = c(0, 1)) +
  labs(title = "A  Mean overlap across 10 elephants",
       subtitle = "Lower: Jaccard of 95% areas  |  Upper: Bhattacharyya coefficient of UDs",
       x = NULL, y = NULL) +
  theme_bw(base_size = 10) + theme(panel.grid = element_blank())

pB <- ggplot(bc, aes(pair, BC)) +
  geom_boxplot(outlier.shape = NA, fill = "grey92", width = 0.5) +
  geom_point(aes(colour = ID == 4501), position = position_jitter(width = 0.08, seed = 1), size = 2) +
  scale_colour_manual(NULL, values = c(`FALSE` = "grey25", `TRUE` = "#B2182B"),
                      labels = c("Other elephants", "4501")) +
  scale_y_continuous("Bhattacharyya coefficient", limits = c(0, 1)) +
  labs(title = "B  UD overlap per elephant", x = NULL) +
  theme_bw(base_size = 10) + theme(legend.position = "bottom")

EX <- 4142
ex_polys <- do.call(rbind, lapply(METHODS, function(m)
  st_sf(method = factor(LABELS[m], lev), geometry = read_poly(m, EX))))
pC <- ggplot(ex_polys) +
  geom_sf(aes(colour = method), fill = NA, linewidth = 0.6) +
  scale_colour_manual(NULL, values = c(AKDE = "#2166AC", dBBMM = "#1B9E77", SAEK = "#D95F02",
                                       TGD = "#7570B3", aLoCoH = "#E7298A")) +
  labs(title = paste0("C  95% areas, elephant ", EX)) +
  theme_bw(base_size = 10) + theme(axis.text = element_text(size = 7))

ggsave(file.path(OUT, "Overlap_figure.png"),
       plot_grid(plot_grid(pA, pB, ncol = 1, rel_heights = c(1.1, 1)), pC, ncol = 2, rel_widths = c(1, 1.1)),
       width = 13, height = 9, dpi = PNG_DPI)

sink(file.path(OUT, "Overlap_summary.txt"), split = TRUE)
cat("Between-estimator overlap, v6 —", format(Sys.time()), "\n\n")
cat("== Grid QC: UD mass captured by the common 250 m grid ==\n"); print(qc, row.names = FALSE)
cat("\n== Bhattacharyya coefficient ==\n"); print(rnd(bc_sum), row.names = FALSE)
cat("\n== Jaccard (95% polygons) ==\n"); print(rnd(jac_sum), row.names = FALSE)
cat("\n== Centroid distance (km) ==\n"); print(rnd(cen_sum), row.names = FALSE)
cat("\n== UD mass inside TGD / aLoCoH polygons ==\n"); print(rnd(mass_sum), row.names = FALSE)
sink()
