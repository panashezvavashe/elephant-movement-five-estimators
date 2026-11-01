# Appendix 9 redrawn from the saved SAEK (formerly ETD) rasters with the new name and
# visible 95%/50% isopleths (revision, 2026-10-31).
Sys.unsetenv(c("PROJ_LIB", "PROJ_DATA", "GDAL_DATA"))
suppressPackageStartupMessages({ library(terra); library(ggplot2); library(cowplot) })
IDS <- c(3905, 3906, 3907, 4141, 4142, 4501, 4556, 4558, 4559, 4560)
OUT <- "final/comparison/v6/figures"
mc  <- read.csv("final/comparison/tables/MethodComparison_Wide.csv")
fx  <- read.csv("data/processed/all_elephants_screened_v6.csv")
vol_threshold <- function(v, p) { s <- sort(v, decreasing = TRUE); s[which(cumsum(s) / sum(s) >= p)[1]] }
panel <- function(id) {
  r <- rast(sprintf("final/ETD/spatial/ETD_%d_full.tif", id))
  d <- as.data.frame(r, xy = TRUE, na.rm = TRUE); names(d)[3] <- "z"
  t95 <- vol_threshold(d$z, 0.95); t50 <- vol_threshold(d$z, 0.50)
  d$rel <- d$z / max(d$z)
  p <- fx[fx$elephant_id == id & fx$speed_km_h <= 5, ]
  e <- ext(r)   # fixes outside the raster are the ones the 45 km filter removed
  p <- p[p$x_utm >= e$xmin & p$x_utm <= e$xmax & p$y_utm >= e$ymin & p$y_utm <= e$ymax, ]
  ggplot(d, aes(x, y)) +
    geom_raster(aes(fill = rel)) +
    scale_fill_gradient(low = "white", high = "#B2182B", guide = "none") +
    geom_point(data = p, aes(x_utm, y_utm), size = 0.05, colour = "grey35", alpha = 0.4) +
    geom_contour(aes(z = z), breaks = t95, colour = "#2166AC", linewidth = 0.5) +
    geom_contour(aes(z = z), breaks = t50, colour = "black", linewidth = 0.4) +
    coord_equal(expand = FALSE) +
    labs(title = sprintf("SAEK — Elephant %d (95%%: %.0f km²)", id, mc$ETD[mc$ID == id]), x = NULL, y = NULL) +
    theme_bw(base_size = 8) + theme(axis.text = element_text(size = 5), plot.title = element_text(size = 8))
}
ps <- lapply(IDS, panel)
g <- plot_grid(plotlist = ps, ncol = 4)
ggsave(file.path(OUT, "App9_SAEK_panels.png"), g, width = 12, height = 9, dpi = 300, bg = "white")
cat("done\n")
