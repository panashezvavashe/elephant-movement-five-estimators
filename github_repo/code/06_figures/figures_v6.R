# ============================================================
# Manuscript figures for the v6 revision (2026-10-30)
#   part A (no seasonal input): Fig 5, Appendix 5, Appendix 10
#   part B (after rerun_seasonal_v6.R): Fig 4, Appendix 11, Appendix 12
#
#   Rscript figures_v6.R A      |      Rscript figures_v6.R B
# ============================================================
Sys.unsetenv(c("PROJ_LIB", "PROJ_DATA", "GDAL_DATA"))
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2); library(ggrepel) })
part <- commandArgs(TRUE)[1]
F   <- "final"
OUT <- file.path(F, "comparison", "v6", "figures")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
DPI <- 300
IDS <- c(3905, 3906, 3907, 4141, 4142, 4501, 4556, 4558, 4559, 4560)
heat <- function(d, xlab, title) {
  d <- d %>% group_by(ID) %>% mutate(rel = (area - min(area)) / max(1e-9, max(area) - min(area))) %>% ungroup()
  ggplot(d, aes(factor(x), factor(ID), fill = rel)) +
    geom_tile(colour = "white") + geom_text(aes(label = round(area)), size = 3) +
    scale_fill_gradient("Relative area\n(within elephant)", low = "#FFF7EC", high = "#D7301F") +
    labs(x = xlab, y = "Elephant ID", title = title) + theme_bw(base_size = 11) +
    theme(panel.grid = element_blank())
}

if (part == "A") {
  mc <- read.csv(file.path(F, "comparison", "tables", "MethodComparison_Wide.csv"))
  long <- mc %>% select(ID, AKDE, dBBMM, ETD, TGD, aLOCOH) %>%
    pivot_longer(-ID, names_to = "method", values_to = "area") %>%
    mutate(method = recode(method, ETD = "SAEK", aLOCOH = "aLoCoH"))
  ord <- long %>% group_by(method) %>% summarise(m = median(area)) %>% arrange(m) %>% pull(method)
  long$method <- factor(long$method, ord)
  p5 <- ggplot(long, aes(method, area, group = factor(ID), colour = factor(ID))) +
    geom_line(linewidth = 0.7) + geom_point(size = 2) +
    scale_y_log10("95% area (km², log scale)") + labs(x = NULL, colour = "Elephant") +
    theme_bw(base_size = 12)
  ggsave(file.path(OUT, "Fig5_area_by_estimator.png"), p5, width = 9, height = 5.5, dpi = DPI)

  ds <- do.call(rbind, lapply(IDS, function(i) cbind(ID = i,
          read.csv(file.path(F, "dBBMM", "spatial", paste0("dBBMM_Sensitivity_", i, ".csv"))))))
  ds$x <- paste0("m=", ds$margin, " / ws=", ds$window_size); ds$area <- ds$area95_km2
  ggsave(file.path(OUT, "App5_dBBMM_sensitivity.png"),
         heat(ds, "Margin / window size", "dBBMM 95% area (km²) across smoothing settings"),
         width = 9, height = 5.5, dpi = DPI)

  es <- read.csv(file.path(F, "ETD", "tables", "ETD_sigma_sensitivity_all_v6.csv"))
  es$x <- es$sigma_m; es$area <- es$area_km2
  ggsave(file.path(OUT, "App10_SAEK_sigma.png"),
         heat(es, "σ (m)", "SAEK 95% area (km²) across σ (σ = 500 m used)"),
         width = 10, height = 5.5, dpi = DPI)
}

if (part == "B") {
  S  <- file.path(F, "AKDE", "seasonal_v6")
  pa <- read.csv(file.path(S, "PerAnimal_meta.csv"))
  pop <- read.csv(file.path(S, "Meta_population.csv"))
  t1 <- read.csv(file.path(F, "data_quality", "tables", "Table1_centroid_v6.csv"))
  w <- pa %>% filter(!is.na(wet_mean)) %>%
    transmute(ID, season = "Wet", est = wet_mean, lo = wet_low, hi = wet_high) %>%
    bind_rows(pa %>% filter(!is.na(dry_mean)) %>%
                transmute(ID, season = "Dry", est = dry_mean, lo = dry_low, hi = dry_high))
  p4 <- ggplot(w, aes(factor(ID), est, colour = season)) +
    geom_pointrange(aes(ymin = lo, ymax = hi), position = position_dodge(width = 0.5)) +
    scale_colour_manual(NULL, values = c(Wet = "#2166AC", Dry = "#B2182B")) +
    labs(x = "Elephant ID", y = "Mean season-year AKDE 95% area (km²)",
         title = "Season-year AKDE combined per elephant (ctmm::meta)") + theme_bw(base_size = 12)
  ggsave(file.path(OUT, "Fig4_seasonal_meta.png"), p4, width = 10, height = 5, dpi = DPI)

  r <- pa %>% filter(!is.na(ratio_wet_dry)) %>%
    transmute(lab = as.character(ID), ratio = ratio_wet_dry, lo = ratio_low, hi = ratio_high, grp = "Elephant") %>%
    bind_rows(pop %>% filter(tag %in% c("all_season_years", "one_per_animal")) %>%
      transmute(lab = recode(tag, all_season_years = "Population (all season-years)",
                             one_per_animal = "Population (one per elephant)"),
                ratio = ratio_wet_dry, lo = ratio_low, hi = ratio_high, grp = "Population"))
  r$lab <- factor(r$lab, rev(r$lab))
  p11 <- ggplot(r, aes(ratio, lab, colour = grp)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_pointrange(aes(xmin = lo, xmax = hi)) + scale_x_log10("Wet:dry area ratio (log scale, 95% CI)") +
    scale_colour_manual(NULL, values = c(Elephant = "grey20", Population = "#B2182B")) +
    labs(y = NULL) + theme_bw(base_size = 12) + theme(legend.position = "none")
  ggsave(file.path(OUT, "App11_ratio_forest.png"), p11, width = 8, height = 5.5, dpi = DPI)

  sc <- pa %>% filter(!is.na(ratio_wet_dry)) %>% left_join(t1, by = "ID")
  p12 <- ggplot(sc, aes(ratio_wet_dry, centroid_shift_sy_mean_km)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_point(aes(shape = CI_excludes_1), size = 3) +
    geom_text_repel(aes(label = ID), size = 3.2, seed = 1) +
    scale_x_log10("Wet:dry area ratio (log scale)") +
    scale_shape_manual("Ratio CI excludes 1", values = c(`FALSE` = 1, `TRUE` = 16)) +
    labs(y = "Mean wet–dry centroid shift per season-year (km)") + theme_bw(base_size = 12)
  ggsave(file.path(OUT, "App12_centroid_vs_ratio.png"), p12, width = 8, height = 5.5, dpi = DPI)
}
cat("figures part", part, "done\n")
