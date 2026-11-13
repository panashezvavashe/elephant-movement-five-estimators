# ============================================================
# Second-round review (Reviewer 1, ms 2029790), quick checks on saved outputs
#   C1. Permutation p-values for Kendall's W (raw and after log(spread))
#   C2. Season-year areas larger than the full-period area: role of N_eff,
#       CI overlap, and variogram asymptote (from the cached season-year fits)
#   C3. Pre-check: do areas grow with tracking duration?
#
# Inputs:  final/comparison/tables/MethodComparison_Wide.csv
#          final/AKDE/seasonal_v6/{SeasonYear_AKDE.csv, FullVsSeasonal_containment.csv, fits/}
#          data/processed/all_elephants_screened_v6.csv
# Output:  final/comparison/r2/
#
#   Rscript code/07_revision2/r2_quick_checks.R
# ============================================================

Sys.unsetenv(c("PROJ_LIB", "PROJ_DATA", "GDAL_DATA"))
suppressPackageStartupMessages({
  library(ctmm); library(sp); library(sf); library(dplyr)
})

DATA_DIR <- "."   # run from the repository root
setwd(DATA_DIR)
CT  <- file.path(DATA_DIR, "final", "comparison", "tables")
SY  <- file.path(DATA_DIR, "final", "AKDE", "seasonal_v6")
OUT <- file.path(DATA_DIR, "final", "comparison", "r2")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
set.seed(20261109)
N_PERM <- 10000L
FIT_PTS <- 1500L

METHODS <- c("AKDE", "dBBMM", "ETD", "TGD", "aLOCOH")
LABELS  <- c(AKDE = "AKDE", dBBMM = "dBBMM", ETD = "SAEK", TGD = "TGD", aLOCOH = "aLoCoH")

# ── Helper functions from combined2.R (as in rerun_seasonal_v6.R) ──
src <- readLines("combined2.R", encoding = "UTF-8")
for (e in parse(text = src, keep.source = FALSE)) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) &&
      is.call(e[[3]]) && identical(e[[3]][[1]], as.name("function")))
    eval(e, envir = globalenv())
}
MAX_PTS <- .Machine$integer.max
season_year <- function(ts) {
  m <- as.integer(format(ts, "%m")); y <- as.integer(format(ts, "%Y"))
  ifelse(m >= 11, paste0("W", y), ifelse(m <= 4, paste0("W", y - 1), paste0("D", y)))
}

df <- read.csv("data/processed/all_elephants_screened_v6.csv")
df$timestamp <- as.POSIXct(df$timestamp, format = "%m/%d/%Y %H:%M", tz = "UTC")
df_clean <- df %>%
  filter(speed_km_h <= 5, !elephant_id %in% c(4557L)) %>%
  arrange(elephant_id, timestamp) %>%
  filter(!(elephant_id == 3905 & timestamp >= as.POSIXct("2023-08-08", tz = "UTC")))

mc <- read.csv(file.path(CT, "MethodComparison_Wide.csv"))

# ── C1. Permutation p for Kendall's W ─────────────────────────
kendall_w <- function(m) {
  r <- apply(m, 2, rank); n <- nrow(r); k <- ncol(r)
  12 * sum((rowSums(r) - k * (n + 1) / 2)^2) / (k^2 * (n^3 - n))
}
# Null: ranks independent across estimators -> permute each column separately
perm_w <- function(m) {
  W  <- kendall_w(m)
  nl <- replicate(N_PERM, kendall_w(apply(m, 2, sample)))
  c(W = W, null_mean = mean(nl), null_q95 = unname(quantile(nl, 0.95)),
    p = (sum(nl >= W) + 1) / (N_PERM + 1))
}
spread <- df_clean %>% group_by(ID = elephant_id) %>%
  summarise(spread_km2 = sqrt(det(cov(cbind(x_utm, y_utm)))) / 1e6, .groups = "drop")
c1 <- do.call(rbind, lapply(c("all", "no4501"), function(s) {
  d <- left_join(if (s == "all") mc else mc[mc$ID != 4501, ], spread, by = "ID")
  raw <- as.matrix(d[, METHODS])
  res <- sapply(METHODS, function(m) resid(lm(log(d[[m]]) ~ log(d$spread_km2))))
  rbind(data.frame(subset = s, measure = "W (areas)",              t(perm_w(raw))),
        data.frame(subset = s, measure = "W after log(spread)",    t(perm_w(res))))
}))
c1[, 3:6] <- round(c1[, 3:6], 4)
write.csv(c1, file.path(OUT, "C1_KendallW_permutation.csv"), row.names = FALSE)

# ── C3 pre-check. Area vs tracking duration ───────────────────
dur <- df_clean %>% group_by(ID = elephant_id) %>%
  summarise(start = min(timestamp), end = max(timestamp),
            months = as.numeric(difftime(end, start, units = "days")) / 30.44, .groups = "drop")
d3 <- left_join(left_join(mc, dur, by = "ID"), spread, by = "ID")
c3 <- do.call(rbind, lapply(METHODS, function(m) {
  r1 <- suppressWarnings(cor.test(d3[[m]], d3$months, method = "spearman"))
  rs <- resid(lm(log(d3[[m]]) ~ log(d3$spread_km2)))
  r2 <- suppressWarnings(cor.test(rs, d3$months, method = "spearman"))
  data.frame(method = LABELS[m], rho_area_vs_months = round(unname(r1$estimate), 3), p = signif(r1$p.value, 3),
             rho_resid_vs_months = round(unname(r2$estimate), 3), p_resid = signif(r2$p.value, 3))
}))
c3 <- rbind(c3, data.frame(method = "spread of fixes",
  rho_area_vs_months = round(cor(d3$spread_km2, d3$months, method = "spearman"), 3), p = NA,
  rho_resid_vs_months = NA, p_resid = NA))
write.csv(c3, file.path(OUT, "C3_area_vs_duration.csv"), row.names = FALSE)
write.csv(dur %>% mutate(across(c(start, end), ~format(.x, "%Y-%m-%d")), months = round(months, 1)),
          file.path(OUT, "C3_track_duration.csv"), row.names = FALSE)

# ── C2. Season-year vs full-period area ───────────────────────
sy   <- read.csv(file.path(SY, "SeasonYear_AKDE.csv"))
full <- read.csv(file.path(SY, "FullVsSeasonal_containment.csv"))[, c("ID", "full_95_km2")]
sy <- left_join(sy, full, by = "ID") %>%
  mutate(larger = area_est > full_95_km2,
         sig_larger = area_low > full_95_km2,          # CI excludes the full-period area
         ratio = area_est / full_95_km2,
         neff_class = ifelse(N_eff < 10, "N_eff < 10", ifelse(N_eff < 30, "10-30", ">= 30")))

# Asymptote: tau_p and its CI from each cached fit; variogram on the same track
sy$tau_p_d <- sy$tau_p_hi_d <- sy$span_over_tau <- NA_real_
pdf(file.path(OUT, "C2_seasonyear_variograms.pdf"), width = 11, height = 8.5)
par(mfrow = c(3, 4), mar = c(4, 4, 2.5, 1))
for (i in seq_len(nrow(sy))) {
  id <- sy$ID[i]; syl <- sy$season_year[i]
  r  <- readRDS(file.path(SY, "fits", sprintf("AKDE_sy_%s_%s.rds", id, syl)))
  ci <- summary(r$fit, units = FALSE)$CI
  k  <- grep("position", rownames(ci))
  if (length(k)) {
    sy$tau_p_d[i]    <- ci[k[1], "est"] / 86400
    sy$tau_p_hi_d[i] <- ci[k[1], "high"] / 86400
    sy$span_over_tau[i] <- sy$span_days[i] / sy$tau_p_d[i]
  }
  sub <- prep_sub(df_clean, id); sub$sy <- season_year(sub$timestamp)
  s   <- thin_to(sub[sub$sy == syl, ], FIT_PTS)
  svf <- variogram(make_tel(s))
  plot(svf, CTMM = r$fit, fraction = 0.65, units = TRUE,
       main = sprintf("%s %s  N_eff %.1f  ratio %.2f", id, syl, sy$N_eff[i], sy$ratio[i]),
       cex.main = 0.85)
}
dev.off()
sy$asymptote <- with(sy, ifelse(is.na(tau_p_d), "no range model",
                       ifelse(is.finite(tau_p_hi_d) & span_over_tau >= 10, "reached (span >= 10 tau_p)",
                       ifelse(is.finite(tau_p_hi_d), "weak (span < 10 tau_p)", "not reached (tau_p CI open)"))))
write.csv(sy, file.path(OUT, "C2_seasonyear_vs_full.csv"), row.names = FALSE)

c2_tab <- sy %>% group_by(neff_class) %>%
  summarise(n = n(), pct_larger = round(100 * mean(larger), 1), larger = sum(larger),
            sig_larger = sum(sig_larger), median_ratio = round(median(ratio), 2), .groups = "drop")
c2_seas <- sy %>% group_by(season) %>%
  summarise(n = n(), pct_larger = round(100 * mean(larger), 1), larger = sum(larger),
            sig_larger = sum(sig_larger), median_ratio = round(median(ratio), 2), .groups = "drop")
c2_asym <- sy %>% group_by(asymptote) %>%
  summarise(n = n(), larger = sum(larger), median_ratio = round(median(ratio), 2), .groups = "drop")
rho_neff <- suppressWarnings(cor.test(log(sy$ratio), log(sy$N_eff), method = "spearman"))
fis <- fisher.test(table(sy$N_eff < 10, sy$larger))
glm1 <- glm(larger ~ log(N_eff) + season, family = binomial, data = sy)

# ── Summary ───────────────────────────────────────────────────
sink(file.path(OUT, "R2_quick_checks_summary.txt"), split = TRUE)
cat("Second-round quick checks —", format(Sys.time()), "\n\n")
cat("== C1. Kendall's W with permutation p (", N_PERM, "column-wise permutations) ==\n"); print(c1, row.names = FALSE)
cat("\n== C3 pre-check. Spearman: area vs tracking duration (months) ==\n"); print(c3, row.names = FALSE)
cat("\nTrack durations (months):\n"); print(dur %>% transmute(ID, months = round(months, 1)), row.names = FALSE)
cat("\n== C2. Season-year area vs full-period AKDE 95% ==\n")
cat(sprintf("Larger than full-period: %d of %d; CI lower bound above full-period: %d\n",
            sum(sy$larger), nrow(sy), sum(sy$sig_larger)))
cat("\nBy N_eff class:\n"); print(as.data.frame(c2_tab), row.names = FALSE)
cat("\nBy season:\n"); print(as.data.frame(c2_seas), row.names = FALSE)
cat("\nBy asymptote class:\n"); print(as.data.frame(c2_asym), row.names = FALSE)
cat(sprintf("\nSpearman log(ratio) vs log(N_eff): rho = %.3f, p = %.3g\n", rho_neff$estimate, rho_neff$p.value))
cat(sprintf("Fisher (N_eff<10 vs larger): OR = %.2f, p = %.3g\n", fis$estimate, fis$p.value))
cat("\nLogistic: larger ~ log(N_eff) + season\n"); print(round(summary(glm1)$coefficients, 3))
cat("\nspan / tau_p: median", round(median(sy$span_over_tau, na.rm = TRUE), 1),
    "range", paste(round(range(sy$span_over_tau, na.rm = TRUE), 1), collapse = "-"), "\n")
sink()
