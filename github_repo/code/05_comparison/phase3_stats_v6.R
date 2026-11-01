# ============================================================
# Phase 3 (revision, 2026-10-29): comparison statistics on the v6 table (rebuilt screening)
#   1. Rank stability: Kendall's W (all five estimators) + bootstrap CI
#   2. Pairwise Spearman rho, bootstrap CI, Holm-adjusted p
#   3. Log-scale Bland-Altman: AKDE vs each other estimator (AKDE vs dBBMM = figure)
#   4. "Mechanical" rank agreement: each estimator vs the spread of the fixes
# Every analysis is run with all ten animals and without 4501.
#
# Input:  final/comparison/tables/MethodComparison_Wide.csv  (frozen)
# Output: final/comparison/v6/
#
#   Rscript phase3_stats_v6.R
# ============================================================

Sys.unsetenv(c("PROJ_LIB", "PROJ_DATA", "GDAL_DATA"))
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(cowplot)
})

DATA_DIR <- "."   # run from the repository root
setwd(DATA_DIR)
CT  <- file.path(DATA_DIR, "final", "comparison", "tables")
OUT <- file.path(DATA_DIR, "final", "comparison", "v6")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
set.seed(20261028)
N_BOOT  <- 10000L
PNG_DPI <- 300L

METHODS <- c("AKDE", "dBBMM", "ETD", "TGD", "aLOCOH")
LABELS  <- c(AKDE = "AKDE", dBBMM = "dBBMM", ETD = "SAEK", TGD = "TGD", aLOCOH = "aLoCoH")
# Range vs occurrence grouping used in the revised paper
GROUP   <- c(AKDE = "range", dBBMM = "occurrence", ETD = "occurrence",
             TGD = "occurrence", aLOCOH = "fix geometry")

mc <- read.csv(file.path(CT, "MethodComparison_Wide.csv"))
subsets <- list(all = mc, no4501 = mc[mc$ID != 4501, ])

# ── 1. Kendall's W ────────────────────────────────────────────
kendall_w <- function(m) {          # m: animals x methods; ranks within each method
  r <- apply(m, 2, rank)
  n <- nrow(r); k <- ncol(r)
  S <- sum((rowSums(r) - k * (n + 1) / 2)^2)
  T <- sum(apply(r, 2, function(x) { t <- table(x); sum(t^3 - t) }))  # tie correction
  12 * S / (k^2 * (n^3 - n) - k * T)
}
w_rows <- lapply(names(subsets), function(s) {
  m <- as.matrix(subsets[[s]][, METHODS]); n <- nrow(m); k <- ncol(m)
  W <- kendall_w(m)
  chi2 <- k * (n - 1) * W
  boot <- replicate(N_BOOT, kendall_w(m[sample(n, n, replace = TRUE), , drop = FALSE]))
  ci <- quantile(boot, c(0.025, 0.975), na.rm = TRUE)
  data.frame(subset = s, n = n, k = k, W = round(W, 3),
             boot_low = round(ci[1], 3), boot_high = round(ci[2], 3),
             chi2 = round(chi2, 2), df = n - 1,
             p = signif(pchisq(chi2, n - 1, lower.tail = FALSE), 3))
})
w_tab <- do.call(rbind, w_rows); rownames(w_tab) <- NULL
write.csv(w_tab, file.path(OUT, "KendallW.csv"), row.names = FALSE)

# ── 2. Pairwise Spearman ──────────────────────────────────────
pairs <- combn(METHODS, 2, simplify = FALSE)
sp_rows <- lapply(names(subsets), function(s) {
  d <- subsets[[s]]; n <- nrow(d)
  out <- do.call(rbind, lapply(pairs, function(p) {
    x <- d[[p[1]]]; y <- d[[p[2]]]
    ct <- suppressWarnings(cor.test(x, y, method = "spearman"))
    b  <- replicate(N_BOOT, { i <- sample(n, n, replace = TRUE)
                              suppressWarnings(cor(x[i], y[i], method = "spearman")) })
    ci <- quantile(b, c(0.025, 0.975), na.rm = TRUE)
    data.frame(subset = s, pair = paste(LABELS[p[1]], "vs", LABELS[p[2]]),
               type = paste(GROUP[p[1]], GROUP[p[2]], sep = "-"), n = n,
               rho = round(unname(ct$estimate), 3),
               boot_low = round(ci[1], 3), boot_high = round(ci[2], 3),
               p = signif(ct$p.value, 3))
  }))
  out$p_holm <- signif(p.adjust(out$p, "holm"), 3)
  out[order(-out$rho), ]
})
sp_tab <- do.call(rbind, sp_rows); rownames(sp_tab) <- NULL
write.csv(sp_tab, file.path(OUT, "Spearman_Holm.csv"), row.names = FALSE)

# ── 3. Log-scale Bland–Altman: AKDE vs each other method ──────
ba_stats <- function(a, b) {
  d <- log(a) - log(b); m <- mean(d); s <- sd(d); n <- length(d)
  se_loa <- s * sqrt(3 / n)                       # Bland & Altman (1999) approx.
  tq <- qt(0.975, n - 1)
  data.frame(n = n,
    ratio_mean = exp(m), ratio_low_CI = exp(m - tq * s / sqrt(n)), ratio_high_CI = exp(m + tq * s / sqrt(n)),
    LoA_low = exp(m - 1.96 * s), LoA_high = exp(m + 1.96 * s),
    LoA_low_CI = paste(round(exp(m - 1.96 * s - tq * se_loa), 2), round(exp(m - 1.96 * s + tq * se_loa), 2), sep = "-"),
    LoA_high_CI = paste(round(exp(m + 1.96 * s - tq * se_loa), 2), round(exp(m + 1.96 * s + tq * se_loa), 2), sep = "-"),
    # proportional bias: does the log ratio change with size?
    trend_p = signif(summary(lm(d ~ I((log(a) + log(b)) / 2)))$coefficients[2, 4], 3))
}
ba_tab <- do.call(rbind, lapply(names(subsets), function(s) {
  d <- subsets[[s]]
  do.call(rbind, lapply(setdiff(METHODS, "AKDE"), function(mth)
    cbind(subset = s, comparison = paste0("AKDE / ", LABELS[mth]), ba_stats(d$AKDE, d[[mth]]))))
}))
num <- sapply(ba_tab, is.numeric) & !names(ba_tab) %in% c("n", "trend_p")
ba_tab[num] <- lapply(ba_tab[num], round, 2)
write.csv(ba_tab, file.path(OUT, "BlandAltman_log.csv"), row.names = FALSE)

ba_plot <- function(d, mth, title) {
  df <- data.frame(ID = d$ID, mean_log = (log(d$AKDE) + log(d[[mth]])) / 2,
                   diff = log(d$AKDE) - log(d[[mth]]))
  m <- mean(df$diff); s <- sd(df$diff)
  brks <- log(c(0.25, 0.5, 1, 2, 4, 8))
  ggplot(df, aes(mean_log, diff)) +
    geom_hline(yintercept = 0, colour = "grey70") +
    geom_hline(yintercept = m, colour = "#2166AC", linewidth = 0.8) +
    geom_hline(yintercept = m + c(-1.96, 1.96) * s, colour = "#B2182B", linetype = "dashed") +
    geom_point(size = 2.4) +
    ggrepel::geom_text_repel(aes(label = ID), size = 3, seed = 1) +
    scale_x_continuous("Mean area, geometric (km²)", breaks = log(c(50, 100, 200, 400, 800, 1600)),
                       labels = function(x) round(exp(x))) +
    scale_y_continuous(paste0("Ratio AKDE / ", LABELS[mth]), breaks = brks,
                       labels = function(x) format(exp(x), digits = 2)) +
    labs(title = title,
         subtitle = sprintf("Mean ratio %.2f×; 95%% limits of agreement %.2f–%.2f×",
                            exp(m), exp(m - 1.96 * s), exp(m + 1.96 * s))) +
    theme_bw(base_size = 11)
}
p_all <- ba_plot(subsets$all, "dBBMM", "A  All ten elephants")
p_no  <- ba_plot(subsets$no4501, "dBBMM", "B  Excluding 4501")
ggsave(file.path(OUT, "BlandAltman_log_AKDE_dBBMM.png"),
       plot_grid(p_all, p_no, ncol = 2), width = 12, height = 5, dpi = PNG_DPI)
ggsave(file.path(OUT, "BlandAltman_log_AKDE_vs_all.png"),
       plot_grid(plotlist = lapply(setdiff(METHODS, "AKDE"), function(m)
         ba_plot(subsets$all, m, LABELS[m])), ncol = 2),
       width = 12, height = 10, dpi = PNG_DPI)

# ── 4. Mechanical agreement: area vs spread of the fixes ──────
# Spread = generalised SD of the screened fixes, sqrt(det(cov(x, y))), in km².
df <- read.csv("data/processed/all_elephants_screened_v6.csv")
df$timestamp <- as.POSIXct(df$timestamp, format = "%m/%d/%Y %H:%M", tz = "UTC")
spread <- df %>% filter(speed_km_h <= 5, elephant_id != 4557,
                        !(elephant_id == 3905 & timestamp >= as.POSIXct("2023-08-08", tz = "UTC"))) %>%
  group_by(ID = elephant_id) %>%
  summarise(spread_km2 = sqrt(det(cov(cbind(x_utm, y_utm)))) / 1e6, .groups = "drop")
mech <- do.call(rbind, lapply(names(subsets), function(s) {
  d <- left_join(subsets[[s]], spread, by = "ID")
  do.call(rbind, lapply(METHODS, function(m) {
    ct <- suppressWarnings(cor.test(d[[m]], d$spread_km2, method = "spearman"))
    data.frame(subset = s, method = LABELS[m], rho_vs_spread = round(unname(ct$estimate), 3),
               p = signif(ct$p.value, 3))
  }))
}))
write.csv(mech, file.path(OUT, "Mechanical_rank_vs_spread.csv"), row.names = FALSE)

# Kendall's W after removing the shared spread component (ranks of residuals
# from log(area) ~ log(spread), per method): what agreement remains beyond
# "bigger cloud of fixes -> bigger area" for every estimator.
d  <- left_join(mc, spread, by = "ID")
res <- sapply(METHODS, function(m) resid(lm(log(d[[m]]) ~ log(d$spread_km2))))
W_res <- kendall_w(res)
b_res <- replicate(N_BOOT, { i <- sample(nrow(res), replace = TRUE); kendall_w(res[i, ]) })
mech_w <- data.frame(W_residual = round(W_res, 3),
                     boot_low = round(quantile(b_res, 0.025), 3),
                     boot_high = round(quantile(b_res, 0.975), 3))
write.csv(mech_w, file.path(OUT, "KendallW_after_spread.csv"), row.names = FALSE)

# ── Summary ───────────────────────────────────────────────────
sink(file.path(OUT, "Phase3_summary.txt"), split = TRUE)
cat("Phase 3 statistics on MethodComparison_Wide.csv —", format(Sys.time()), "\n\n")
cat("== Kendall's W (5 estimators) ==\n"); print(w_tab, row.names = FALSE)
cat("\n== Kendall's W on residuals after log(spread) ==\n"); print(mech_w, row.names = FALSE)
cat("\n== Spearman (Holm over 10 pairs) ==\n"); print(sp_tab, row.names = FALSE)
cat("\n== Log Bland–Altman (AKDE / other) ==\n"); print(ba_tab, row.names = FALSE)
cat("\n== Area rank vs spread of fixes ==\n"); print(mech, row.names = FALSE)
sink()
