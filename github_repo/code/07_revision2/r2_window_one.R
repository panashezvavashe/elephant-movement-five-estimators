# ============================================================
# Second-round review, Reviewer 1 comment 3 (2026-11-09):
# rerun one estimator for one animal on the COMMON tracking window only
# (latest track start to earliest track end across the ten animals), so that
# differences in tracking length (34-62 months) cannot drive the ranking.
#
#   Rscript --max-ppsize=500000 code/07_revision2/r2_window_one.R <AKDE|E|F|G|G2> <elephant_id>
#     AKDE = full-track ctmm model selection + AKDE (as full_track_ctmm_v6.R)
#     E = dBBMM, F = aLoCoH, G = TGD, G2 = SAEK (section text verbatim from combined2.R)
#
# All outputs go to final/r2_window/ (nothing under final/ is overwritten);
# the 95% areas are written to final/r2_window/areas/Area_<sec>_<id>.csv
# ============================================================

args <- commandArgs(trailingOnly = TRUE)
sec <- args[1]; one_id <- as.integer(args[2])

src <- readLines("code/02_estimators/combined2.R", encoding = "UTF-8")
line_of <- function(pat) { i <- grep(pat, src, fixed = TRUE)[1]; stopifnot(!is.na(i)); i }

# Set-up, helpers and data load (sections 0, A, B) with two changes:
# outputs redirected to final/r2_window, and the legacy "migrate" block removed
setup <- src[1:(line_of("# SECTION C — DATA QUALITY") - 2)]
i_fd  <- grep('FINAL_DIR <- file.path(DATA_DIR, "final")', setup, fixed = TRUE)
stopifnot(length(i_fd) == 1)
setup[i_fd] <- 'FINAL_DIR <- file.path(DATA_DIR, "final", "r2_window")'
m0 <- grep("# ── Migrate pre-existing results", setup, fixed = TRUE)
m1 <- grep('cat("Migration of existing outputs', setup, fixed = TRUE)
stopifnot(length(m0) == 1, length(m1) == 1)
setup <- setup[-(m0:m1)]
eval(parse(text = setup), envir = globalenv())
N_CORES <- 1L   # jobs run four at a time in parallel

# Common window across all ten animals
rng <- tapply(df_clean$timestamp, df_clean$elephant_id, range)
win_start <- max(do.call(c, lapply(rng, `[`, 1)))
win_end   <- min(do.call(c, lapply(rng, `[`, 2)))
df_clean <- df_clean[df_clean$elephant_id == one_id &
                     df_clean$timestamp >= win_start & df_clean$timestamp <= win_end, ]
cat(sprintf("== %s, elephant %d: window %s to %s, %d fixes ==\n", sec, one_id,
            format(win_start, "%Y-%m-%d %H:%M"), format(win_end, "%Y-%m-%d %H:%M"), nrow(df_clean)))

AR <- file.path(FINAL_DIR, "areas"); dir.create(AR, recursive = TRUE, showWarnings = FALSE)
t0 <- Sys.time()

if (sec == "AKDE") {
  sub  <- prep_sub(df_clean, one_id, max_pts = Inf)
  tel  <- make_tel(sub)
  fits <- ctmm.select(tel, ctmm.guess(tel, interactive = FALSE), cores = N_CORES, verbose = TRUE)
  fit  <- fits[[1]]
  ud   <- akde(tel, fit); a <- summary(ud, units = FALSE)$CI[1, ] / 1e6
  saveRDS(list(fits = fits, ud = ud), file.path(AS, paste0("AKDE_window_", one_id, ".rds")))
  model_comparison <- data.frame(ID = one_id, Method = "AKDE", Area_km2 = a[2], AIC = fit$AIC,
                                 area_low = a[1], area_high = a[3],
                                 model = summary(fit)$name, N_eff = summary(ud)$DOF["area"])
} else {
  bounds <- list(
    E  = c("# SECTION E — dBBMM",  "# SECTION F — aLOCOH"),
    F  = c("# SECTION F — aLOCOH", "# SECTION G — TGD"),
    G  = c("# SECTION G — TGD",    "# SECTION G2 — ETD"),
    G2 = c("# SECTION G2 — ETD",   "# SECTION H — SEASONAL"))[[sec]]
  eval(parse(text = src[(line_of(bounds[1]) - 1):(line_of(bounds[2]) - 2)]), envir = globalenv())
}

out <- model_comparison[model_comparison$Method %in% c("AKDE", "dBBMM", "aLOCOH", "TGD", "ETD"), ]
out$n_fixes   <- nrow(df_clean)
out$win_start <- format(win_start, "%Y-%m-%d")
out$win_end   <- format(win_end, "%Y-%m-%d")
out$mins      <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
write.csv(out, file.path(AR, sprintf("Area_%s_%d.csv", sec, one_id)), row.names = FALSE)
cat("DONE", sec, one_id, "\n")
