# Regenerate Appendix 1 (diurnal speed profile, local time) from combined2.R's
# Section K code, and summarise the population profile (revision, 2026-10-30).
#   Rscript diurnal_v6.R
src <- readLines("code/02_estimators/combined2.R", encoding = "UTF-8")
line_of <- function(pat) { i <- grep(pat, src, fixed = TRUE)[1]; stopifnot(!is.na(i)); i }
eval(parse(text = src[1:(line_of("# SECTION C — DATA QUALITY") - 2)]), envir = globalenv())
eval(parse(text = src[line_of("# Local time (Africa/Harare, UTC+2)"):line_of("ggsave(file.path(MG, \"Diurnal_Pattern.png\")")]),
     envir = globalenv())
pop <- diurnal %>% filter(elephant_id != 4501) %>% group_by(hour) %>%
  summarise(mean_kmh = round(mean(mean_sl), 3), .groups = "drop")
write.csv(pop, file.path(MT, "Diurnal_population_v6.csv"), row.names = FALSE)
print(as.data.frame(pop))
