# ============================================================
# Run one estimator section of combined2.R for one animal (revision, 2026-10-30)
# Lets the slow per-animal sections (TGD, ETD) run as parallel processes.
# The section text is taken verbatim from combined2.R.
#
#   Rscript --max-ppsize=500000 run_section_one_id.R <F|G|G2> <elephant_id>
# Writes the section's usual per-animal outputs; for G2 also
# final/ETD/tables/ETD_Comparison_<id>.csv
# ============================================================

args <- commandArgs(trailingOnly = TRUE)
sec <- args[1]; one_id <- as.integer(args[2])

src <- readLines("code/02_estimators/combined2.R", encoding = "UTF-8")
line_of <- function(pat) { i <- grep(pat, src, fixed = TRUE)[1]; stopifnot(!is.na(i)); i }
c_start <- line_of("# SECTION C — DATA QUALITY")
bounds <- list(
  F  = c("# SECTION F — aLOCOH", "# SECTION G — TGD"),
  G  = c("# SECTION G — TGD",    "# SECTION G2 — ETD"),
  G2 = c("# SECTION G2 — ETD",   "# SECTION H — SEASONAL"))[[sec]]

eval(parse(text = src[1:(c_start - 2)]), envir = globalenv())    # set-up, A, B
df_clean <- df_clean[df_clean$elephant_id == one_id, ]
cat(sprintf("== Section %s, elephant %d (%d fixes) ==\n", sec, one_id, nrow(df_clean)))
eval(parse(text = src[(line_of(bounds[1]) - 1):(line_of(bounds[2]) - 2)]), envir = globalenv())

if (sec == "G2" && nrow(etd_comparison) > 0)
  write.csv(etd_comparison, file.path(ET, paste0("ETD_Comparison_", one_id, ".csv")), row.names = FALSE)
cat("DONE", sec, one_id, "\n")
