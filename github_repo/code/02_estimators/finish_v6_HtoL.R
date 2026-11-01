# ============================================================
# Finish combined2.R from saved per-estimator outputs (revision, 2026-10-30)
# Sections D (AKDE), E (dBBMM) and F (aLoCoH) ran in resume_v6_from_dBBMM.R; TGD,
# ETD and the 4558 aLoCoH ran per animal in run_section_one_id.R. This script runs
# the set-up and Sections A-C, rebuilds model_comparison / etd_results from the
# saved files, then runs Section H to the end of combined2.R unchanged.
#
#   Rscript finish_v6_HtoL.R
# ============================================================

src <- readLines("code/02_estimators/combined2.R", encoding = "UTF-8")
line_of <- function(pat) { i <- grep(pat, src, fixed = TRUE)[1]; stopifnot(!is.na(i)); i }
eval(parse(text = src[1:(line_of("# SECTION D — AKDE") - 2)]), envir = globalenv())

poly_area <- function(f) {
  if (!file.exists(f)) return(NA_real_)
  sum(as.numeric(sf::st_area(sf::st_read(f, quiet = TRUE)))) / 1e6
}
param_table <- read.csv(file.path(AT, "AKDE_Parameters.csv"))
etd_rows <- list()
for (id in unique(df_clean$elephant_id)) {
  k  <- as.character(id)
  ud <- readRDS(file.path(AS, paste0("AKDE_", id, ".rds")))
  ft <- readRDS(file.path(AS, paste0("CTMM_fit_", id, ".rds")))
  akde_results[[k]] <- ud; ctmm_fits[[k]] <- ft
  ci <- summary(ud)$CI
  ci_df <- rbind(ci_df, data.frame(ID = id, low = ci[1], est = ci[2], high = ci[3], CI_width = ci[3] - ci[1]))
  rows <- data.frame(ID = id, Method = "AKDE", Area_km2 = ci[2], AIC = ft$AIC)
  for (s in c("wet", "dry")) {
    f <- file.path(AS, paste0("AKDE_", s, "_", id, ".rds"))
    if (!file.exists(f)) next
    u <- readRDS(f)
    if (s == "wet") akde_wet[[k]] <- u else akde_dry[[k]] <- u
    rows <- rbind(rows, data.frame(ID = id, Method = if (s == "wet") "AKDE_Wet" else "AKDE_Dry",
                                   Area_km2 = summary(u)$CI[2], AIC = NA))
  }
  rows <- rbind(rows,
    data.frame(ID = id, Method = "dBBMM",  Area_km2 = poly_area(file.path(DS, paste0("dBBMM_95_", id, ".shp"))),  AIC = NA),
    data.frame(ID = id, Method = "aLOCOH", Area_km2 = poly_area(file.path(LS, paste0("aLOCOH_95_", id, ".shp"))), AIC = NA),
    data.frame(ID = id, Method = "TGD",    Area_km2 = poly_area(file.path(TS, paste0("TGD_", id, ".shp"))),       AIC = NA))
  ef <- file.path(ET, paste0("ETD_Comparison_", id, ".csv"))
  if (file.exists(ef)) {
    e <- read.csv(ef); etd_rows[[k]] <- e
    rows <- rbind(rows, data.frame(ID = id, Method = "ETD", Area_km2 = e$Area_95, AIC = NA))
    etd_results[[k]] <- raster::raster(file.path(ES, paste0("ETD_", id, "_full.tif")))
  }
  model_comparison <- dplyr::bind_rows(model_comparison, rows)
}
etd_comparison <- do.call(rbind, etd_rows)
write.csv(etd_comparison, file.path(ET, "ETD_Comparison.csv"), row.names = FALSE)
cat("\n=== Rebuilt model_comparison ===\n"); print(model_comparison)

# Optional first argument "K": resume at Section K (H-J outputs already written)
start <- if (length(commandArgs(TRUE)) && commandArgs(TRUE)[1] == "K") "# SECTION K — MOVEMENT" else "# SECTION H — SEASONAL"
ids_with_akde <- names(akde_results)
eval(parse(text = src[(line_of(start) - 1):length(src)]), envir = globalenv())
