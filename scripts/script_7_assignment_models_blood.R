# Blood plasma assignment models
# Three-month window, 75% latitude prediction interval, and TDFs of 1.1 (d13C) and 1.4 (d15N).

# Load libraries
pacman::p_load(dplyr, readr, here)
source(here("scripts/source_script_assignment_functions.R"))

# Needed objects ####
e <- load_assignment_inputs()

# Plasma samples ####
blood <- read_csv(here("data/data_cleaning/tarpon_iso_corrected.csv"), show_col_types = FALSE) %>%
  filter(sample_material == "plasma") %>%
  arrange(ID)

# Run assignments ####
raster_dir <- here("output/blood/rasters")
dir.create(raster_dir, recursive = TRUE, showWarnings = FALSE)

blood_assign <- assign_samples_parallel(e, blood, Troph.fract.C = -1.1, Troph.fract.N = -1.4,
                                        lat_ci = "75", month_count = 3, raster_dir = raster_dir, seed = 1)

blood_final <- blood %>%
  left_join(blood_assign, by = "ID") %>%
  mutate(Region_Original_fine = categorize_location(Latitude, Longitude),
         Region_New_fine = categorize_location(Latitude_new, Longitude_new),
         dist = round(dist, 2))

# Save ####
write_csv(blood_final, here("output/blood/blood_final.csv"))
