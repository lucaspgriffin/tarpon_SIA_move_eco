# Run the full analysis in order. The fin and plasma assignment scripts write the probability rasters and
# result tables that the later scripts read.
library(here)
scripts <- c("script_1_tar_sia_f.R", "script_6_assignment_models_fin.R", "script_7_assignment_models_blood.R",
             "script_6_7_assignment_models_output_summary.R", "script_6_7_sensitivity.R", "manuscript_figures.R")
for (s in scripts) {
  message("Running ", s)
  source(here("scripts", s), echo = FALSE)
}
writeLines(capture.output(sessionInfo()), here("output/session_info.txt"))  # R and package versions used
