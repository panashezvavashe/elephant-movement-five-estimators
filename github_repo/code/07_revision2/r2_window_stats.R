# ============================================================
# Second-round review, Reviewer 1 comment 3 (2026-11-09):
# rank agreement on the common tracking window (outputs of r2_window_one.R)
#   - Kendall's W (bootstrap CI, permutation p), with and without 4501
#   - pairwise Spearman with bootstrap CI and Holm correction
#   - W after removing log(spread) of the window fixes, permutation p
#   - window vs full-period rank for each estimator
#
# Input:  final/r2_window/areas/Area_<sec>_<id>.csv
#         final/comparison/tables/MethodComparison_Wide.csv (full period)
# Output: final/r2_window/R2_window_summary.txt and CSVs
#
#   Rscript code/07_revision2/r2_window_stats.R
# ============================================================

suppressPackageStartupMessages(library(dplyr))
DATA_DIR <- "."   # run from the repository root
setwd(DATA_DIR)
WD <- file.path("final", "r2_window")
set.seed(20261109)
N_BOOT <- 10000L; N_PERM <- 10000L
METHODS <- c("AKDE", "dBBMM", "ETD", "TGD", "aLOCOH")
LABELS  <- c(AKDE = "AKDE", dBBMM = "dBBMM", ETD = "SAEK", TGD = "TGD", aLOCOH = "aLoCoH")

a <- do.call(rbind, lapply(list.files(file.path(WD, "areas"), full.names = TRUE), function(f)
  read.csv(f)[, c("ID", "Method", "Area_km2", "n_fixes", "win_start", "win_end")]))
win <- tidyr::pivot_wider(a[, c("ID", "Method", "Area_km2")], names_from = Method, values_from = Area_km2)
win <- as.data.frame(win[order(win$ID), c("ID", METHODS)])
stopifnot(nrow(win) == 10, !anyNA(win))
write.csv(win, file.path(WD, "MethodComparison_Window.csv"), row.names = FALSE)
full <- read.csv(file.path("final", "comparison", "tables", "MethodComparison_Wide.csv"))[, c("ID", METHODS)]

kendall_w <- function(m) {
  r <- apply(m, 2, rank); n <- nrow(r); k <- ncol(r)
  S <- sum((rowSums(r) - k * (n + 1) / 2)^2)
  T <- sum(apply(r, 2, function(x) { t <- table(x); sum(t^3 - t) }))
  12 * S / (k^2 * (n^3 - n) - k * T)
}
w_stats <- function(m, label) {
  W <- kendall_w(m)
  b <- replicate(N_BOOT, kendall_w(m[sample(nrow(m), replace = TRUE), , drop = FALSE]))
  nl <- replicate(N_PERM, kendall_w(apply(m, 2, sample)))
  data.frame(measure = label, n = nrow(m), W = round(W, 3),
             boot_low = round(quantile(b, 0.025, na.rm = TRUE), 3),
             boot_high = round(quantile(b, 0.975, na.rm = TRUE), 3),
             p_perm = signif((sum(nl >= W) + 1) / (N_PERM + 1), 3))
}

# Spread of the window fixes (same filters as the estimators' input)
df <- read.csv("data/processed/all_elephants_screened_v6.csv")
df$timestamp <- as.POSIXct(df$timestamp, format = "%m/%d/%Y %H:%M", tz = "UTC")
ws <- as.POSIXct(a$win_start[1], tz = "UTC"); we <- as.POSIXct(a$win_end[1], tz = "UTC") + 86400
spread <- df %>% filter(speed_km_h <= 5, elephant_id != 4557, timestamp >= ws, timestamp < we) %>%
  group_by(ID = elephant_id) %>%
  summarise(spread_km2 = sqrt(det(cov(cbind(x_utm, y_utm)))) / 1e6, .groups = "drop")

w_tab <- list(); sp_tab <- list()
for (s in c("all", "no4501")) {
  d <- if (s == "all") win else win[win$ID != 4501, ]
  d <- left_join(d, spread, by = "ID")
  m <- as.matrix(d[, METHODS])
  res <- sapply(METHODS, function(k) resid(lm(log(d[[k]]) ~ log(d$spread_km2))))
  w_tab[[s]] <- cbind(subset = s, rbind(w_stats(m, "W (window areas)"),
                                        w_stats(res, "W after log(spread)")))
  pairs <- combn(METHODS, 2, simplify = FALSE)
  sp <- do.call(rbind, lapply(pairs, function(p) {
    ct <- suppressWarnings(cor.test(d[[p[1]]], d[[p[2]]], method = "spearman", exact = TRUE))
    b  <- replicate(N_BOOT, { i <- sample(nrow(d), replace = TRUE)
      suppressWarnings(cor(d[[p[1]]][i], d[[p[2]]][i], method = "spearman")) })
    data.frame(subset = s, pair = paste(LABELS[p[1]], "vs", LABELS[p[2]]),
               rho = round(unname(ct$estimate), 3),
               boot_low = round(quantile(b, 0.025, na.rm = TRUE), 3),
               boot_high = round(quantile(b, 0.975, na.rm = TRUE), 3), p = ct$p.value)
  }))
  sp$p_holm <- signif(p.adjust(sp$p, "holm"), 3); sp$p <- signif(sp$p, 3)
  sp_tab[[s]] <- sp
}
w_tab <- do.call(rbind, w_tab); sp_tab <- do.call(rbind, sp_tab)
rownames(w_tab) <- rownames(sp_tab) <- NULL
write.csv(w_tab, file.path(WD, "KendallW_window.csv"), row.names = FALSE)
write.csv(sp_tab, file.path(WD, "Spearman_Holm_window.csv"), row.names = FALSE)

# Window vs full-period, per estimator
fw <- inner_join(full, win, by = "ID", suffix = c("_full", "_win"))
vs <- do.call(rbind, lapply(METHODS, function(k) {
  ct <- suppressWarnings(cor.test(fw[[paste0(k, "_full")]], fw[[paste0(k, "_win")]], method = "spearman"))
  data.frame(method = LABELS[k], rho_window_vs_full = round(unname(ct$estimate), 3),
             p = signif(ct$p.value, 3),
             median_ratio_window_to_full = round(median(fw[[paste0(k, "_win")]] / fw[[paste0(k, "_full")]]), 2))
}))
write.csv(vs, file.path(WD, "Window_vs_full.csv"), row.names = FALSE)
ranks <- data.frame(ID = win$ID, sapply(METHODS, function(k) rank(win[[k]])))
names(ranks)[-1] <- LABELS[METHODS]
write.csv(ranks, file.path(WD, "Ranks_window.csv"), row.names = FALSE)

sink(file.path(WD, "R2_window_summary.txt"), split = TRUE)
cat("Common-window rank agreement —", format(Sys.time()), "\n")
cat("Window:", a$win_start[1], "to", a$win_end[1], "\n\n")
win_r <- win; win_r[-1] <- round(win_r[-1], 1)
cat("== 95% areas in the window (km²) ==\n"); print(win_r, row.names = FALSE)
nf <- unique(a[, c("ID", "n_fixes")])
cat("\n== Fixes per animal in the window ==\n"); print(nf[order(nf$ID), ], row.names = FALSE)
cat("\n== Kendall's W ==\n"); print(w_tab, row.names = FALSE)
cat("\n== Spearman (Holm over 10 pairs) ==\n"); print(sp_tab, row.names = FALSE)
cat("\n== Window vs full period, per estimator ==\n"); print(vs, row.names = FALSE)
cat("\n== Ranks in the window (1 = smallest) ==\n"); print(ranks, row.names = FALSE)
sink()
