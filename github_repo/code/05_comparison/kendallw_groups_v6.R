# Kendall's W within estimator groups (Reviewer 1, comment 1; revision 2026-10-31)
set.seed(20261028); N_BOOT <- 10000L
mc <- read.csv("final/comparison/tables/MethodComparison_Wide.csv")
kendall_w <- function(m) {
  r <- apply(m, 2, rank); n <- nrow(r); k <- ncol(r)
  S <- sum((rowSums(r) - k * (n + 1) / 2)^2)
  T <- sum(apply(r, 2, function(x) { t <- table(x); sum(t^3 - t) }))
  12 * S / (k^2 * (n^3 - n) - k * T)
}
G <- list(occurrence = c("dBBMM", "ETD", "TGD"), UD = c("AKDE", "dBBMM", "ETD"),
          all5 = c("AKDE", "dBBMM", "ETD", "TGD", "aLOCOH"))
out <- do.call(rbind, lapply(names(G), function(g) do.call(rbind, lapply(c("all", "no4501"), function(s) {
  d <- if (s == "all") mc else mc[mc$ID != 4501, ]
  m <- as.matrix(d[, G[[g]]]); n <- nrow(m); k <- ncol(m); W <- kendall_w(m)
  b <- replicate(N_BOOT, kendall_w(m[sample(n, n, TRUE), , drop = FALSE]))
  chi2 <- k * (n - 1) * W
  data.frame(group = g, subset = s, methods = paste(G[[g]], collapse = "+"), n, k, W = round(W, 3),
             boot_low = round(quantile(b, .025, na.rm = TRUE), 3), boot_high = round(quantile(b, .975, na.rm = TRUE), 3),
             p = signif(pchisq(chi2, n - 1, lower.tail = FALSE), 3))
}))))
rownames(out) <- NULL; print(out)
write.csv(out, "final/comparison/v6/KendallW_groups.csv", row.names = FALSE)
