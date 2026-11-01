# How much of each dBBMM isopleth lies far (> 2 km) from any fix? Flags the
# circular gap artefacts seen in Appendix 4 (revision, 2026-10-31).
Sys.unsetenv(c("PROJ_LIB", "PROJ_DATA", "GDAL_DATA"))
suppressPackageStartupMessages({ library(sf); library(dplyr) })
df <- read.csv("data/processed/all_elephants_screened_v6.csv")
out <- list()
for (id in c(3905, 3906, 3907, 4141, 4142, 4501, 4556, 4558, 4559, 4560)) {
  s <- df[df$elephant_id == id & df$speed_km_h <= 5, ]
  pts <- st_as_sf(s, coords = c("x_utm", "y_utm"), crs = 32736)
  near <- st_union(st_buffer(pts, 2000))
  for (lev in c(95, 99)) {
    p <- st_union(st_read(sprintf("final/dBBMM/spatial/dBBMM_%d_%d.shp", lev, id), quiet = TRUE))
    p <- st_transform(p, 32736)
    a <- as.numeric(st_area(p)) / 1e6
    far <- sum(as.numeric(st_area(st_difference(p, near)))) / 1e6
    out[[length(out) + 1]] <- data.frame(ID = id, level = lev, area_km2 = round(a, 1),
      far2km_km2 = round(far, 1), pct_far = round(100 * far / a, 1))
  }
}
r <- do.call(rbind, out); print(r)
write.csv(r, "final/dBBMM/tables/dBBMM_artefact_check_v6.csv", row.names = FALSE)
