# Sensitivity (Reviewer 1, comment 3): repeat the wet:dry meta-analysis without
# season-year fits with N_eff < 10 (revision, 2026-10-31).
Sys.unsetenv(c("PROJ_LIB", "PROJ_DATA", "GDAL_DATA"))
suppressPackageStartupMessages({ library(ctmm); library(dplyr) })
# run from the repository root
OUT <- "final/AKDE/seasonal_v6"
sy_tab <- read.csv(file.path(OUT, "SeasonYear_AKDE.csv"))
key <- paste(sy_tab$ID, sy_tab$season_year, sep = "_")
uds <- setNames(lapply(key, function(k) readRDS(file.path(OUT, "fits", paste0("AKDE_sy_", k, ".rds")))$ud), key)
src <- readLines("code/04_seasonal/rerun_seasonal_v6.R", encoding = "UTF-8")
i0 <- grep("^meta_summary <- function", src); i1 <- grep("^run_meta <- function", src) - 1
eval(parse(text = src[i0:i1]))
ok <- sy_tab$N_eff >= 10; w <- sy_tab$season == "wet"
cat("kept", sum(ok), "of", nrow(sy_tab), "\n")
pop <- rbind(meta_summary(key[w & ok], key[!w & ok], "all_season_years_neff10"),
             meta_summary(key[w & ok & sy_tab$ID != 4501], key[!w & ok & sy_tab$ID != 4501], "no4501_neff10"))
pa <- bind_rows(lapply(unique(sy_tab$ID), function(id) {
  kw <- key[w & ok & sy_tab$ID == id]; kd <- key[!w & ok & sy_tab$ID == id]
  if (!length(kw) || !length(kd)) return(data.frame(ID = id, n_wet = length(kw), n_dry = length(kd)))
  s <- meta_summary(kw, kd, id); s$ID <- id
  s$CI_excludes_1 <- !is.na(s$ratio_low) & (s$ratio_low > 1 | s$ratio_high < 1)
  s[, c("ID", "n_wet", "n_dry", "ratio_wet_dry", "ratio_low", "ratio_high", "CI_excludes_1")]
}))
print(pop[, names(pop) != "log"]); print(pa)
write.csv(pop[, names(pop) != "log"], file.path(OUT, "Meta_population_neff10.csv"), row.names = FALSE)
write.csv(pa, file.path(OUT, "PerAnimal_meta_neff10.csv"), row.names = FALSE)
write.csv(sy_tab[!ok, c("ID", "season_year", "n_fixes", "N_eff", "model")], file.path(OUT, "SeasonYear_neff_lt10.csv"), row.names = FALSE)
