# ============================================================
# Resume combined2.R after Section D (revision, 2026-10-30)
# Section D (AKDE) completed on the v6 data at 14:57; the run was stopped in
# Section E to switch dBBMM/TGD/ETD to full-resolution tracks. This script runs
# combined2.R's set-up, Sections A-C, reloads the saved Section D outputs, then
# runs Section E to the end of the file unchanged.
#
#   Rscript resume_v6_from_dBBMM.R
# ============================================================

src <- readLines("code/02_estimators/combined2.R", encoding = "UTF-8")
line_of <- function(pat) { i <- grep(pat, src, fixed = TRUE)[1]; stopifnot(!is.na(i)); i }
d_start <- line_of("# SECTION D — AKDE")
e_start <- line_of("# SECTION E — dBBMM")

eval(parse(text = src[1:(d_start - 2)]), envir = globalenv())   # set-up, A, B, C

# ── Reload Section D outputs ──────────────────────────────────
cat("\n=== Reloading Section D (AKDE) outputs ===\n")
param_table <- read.csv(file.path(AT, "AKDE_Parameters.csv"))
for (id in unique(df_clean$elephant_id)) {
  k  <- as.character(id)
  ud <- readRDS(file.path(AS, paste0("AKDE_", id, ".rds")))
  ft <- readRDS(file.path(AS, paste0("CTMM_fit_", id, ".rds")))
  akde_results[[k]] <- ud; ctmm_fits[[k]] <- ft
  ci <- summary(ud)$CI
  model_comparison <- dplyr::bind_rows(model_comparison,
    data.frame(ID = id, Method = "AKDE", Area_km2 = ci[2], AIC = ft$AIC))
  ci_df <- rbind(ci_df, data.frame(ID = id, low = ci[1], est = ci[2], high = ci[3], CI_width = ci[3] - ci[1]))
  for (s in c("wet", "dry")) {
    f <- file.path(AS, paste0("AKDE_", s, "_", id, ".rds"))
    if (!file.exists(f)) next
    u <- readRDS(f)
    if (s == "wet") akde_wet[[k]] <- u else akde_dry[[k]] <- u
    model_comparison <- dplyr::bind_rows(model_comparison,
      data.frame(ID = id, Method = if (s == "wet") "AKDE_Wet" else "AKDE_Dry",
                 Area_km2 = summary(u)$CI[2], AIC = NA))
  }
}
print(model_comparison)

eval(parse(text = src[(e_start - 1):length(src)]), envir = globalenv())   # E to end
