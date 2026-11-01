# Sensitivity: ctmm model selection + AKDE on the FULL screened tracks (no thinning),
# to test whether thinning to 1,500 fixes hid OUF velocity autocorrelation (revision, 2026-10-31).
#   Rscript full_track_ctmm_v6.R [id ...]
src <- readLines("code/02_estimators/combined2.R", encoding = "UTF-8")
line_of <- function(pat) { i <- grep(pat, src, fixed = TRUE)[1]; stopifnot(!is.na(i)); i }
eval(parse(text = src[1:(line_of("# SECTION C — DATA QUALITY") - 2)]), envir = globalenv())
OUTD <- file.path("final", "AKDE", "fulltrack_v6"); dir.create(OUTD, showWarnings = FALSE, recursive = TRUE)
ids <- as.integer(commandArgs(TRUE)); if (!length(ids)) ids <- sort(unique(df_clean$elephant_id))
tau_row <- function(fit, nm) {
  s <- summary(fit, units = FALSE)$CI; r <- grep(nm, rownames(s))
  if (length(r)) s[r[1], ] / 3600 else c(NA, NA, NA)
}
for (id in ids) {
  f_out <- file.path(OUTD, paste0("FullTrack_", id, ".csv"))
  if (file.exists(f_out)) next
  t0  <- Sys.time()
  sub <- prep_sub(df_clean, id, max_pts = Inf)
  th  <- prep_sub(df_clean, id)
  tel <- make_tel(sub)
  fits <- ctmm.select(tel, ctmm.guess(tel, interactive = FALSE), cores = N_CORES, verbose = TRUE)
  fit <- fits[[1]]; nm <- sapply(fits, function(f) summary(f)$name)
  aic <- sapply(fits, function(f) f$AIC)
  ud  <- akde(tel, fit); a <- summary(ud, units = FALSE)$CI[1, ] / 1e6
  tp  <- tau_row(fit, "position"); tv <- tau_row(fit, "velocity")
  saveRDS(list(fits = fits, ud = ud), file.path(OUTD, paste0("FullTrack_", id, ".rds")))
  res <- data.frame(ID = id, n_full = nrow(sub), n_thin = nrow(th),
    med_dt_full_h = median(diff(as.numeric(sub$timestamp))) / 3600,
    med_dt_thin_h = median(diff(as.numeric(th$timestamp))) / 3600,
    best = nm[1], second = nm[2], dAIC_second = round(aic[2] - aic[1], 2),
    tau_p_h = tp[2], tau_p_lo = tp[1], tau_p_hi = tp[3],
    tau_v_h = tv[2], tau_v_lo = tv[1], tau_v_hi = tv[3],
    area95 = a[2], area95_lo = a[1], area95_hi = a[3],
    N_eff_area = summary(ud)$DOF["area"],
    mins = as.numeric(difftime(Sys.time(), t0, units = "mins")))
  write.csv(res, f_out, row.names = FALSE); print(res)
}
