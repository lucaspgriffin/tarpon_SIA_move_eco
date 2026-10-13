# Fin tissue assignment models
# Five-month window, 95% latitude prediction interval, and TDFs of 3.0 (d13C) and 3.5 (d15N). Uses the
# ethanol-adjusted fin values from script_1_tar_sia_f.R.

# Load libraries
pacman::p_load(dplyr, readr, here)
source(here("scripts/source_script_assignment_functions.R"))

# Needed objects ####
e <- load_assignment_inputs()

# Fin samples ####
fin <- read_csv(here("data/data_cleaning/tarpon_iso_corrected.csv"), show_col_types = FALSE) %>%
  filter(sample_material == "fin") %>%
  arrange(ID)
fin %>% summarise(mean(d13C), sd(d13C), min(d13C), max(d13C), mean(d15N), sd(d15N), min(d15N), max(d15N))

# Variances ####
# Tarpon, lizardfish, calibration and pooled variances (Supplementary Materials, Sum of variances)
variance_fin <- summarise_variance(variance_terms(e, fin, Troph.fract.C = -3, Troph.fract.N = -3.5, seed = 1))
variance_fin

# Run assignments ####
raster_dir <- here("output/fin/rasters")
dir.create(raster_dir, recursive = TRUE, showWarnings = FALSE)

fin_assign <- assign_samples_parallel(e, fin, Troph.fract.C = -3, Troph.fract.N = -3.5,
                                      lat_ci = "95", month_count = 5, raster_dir = raster_dir, seed = 1)

finclip_final <- fin %>%
  left_join(fin_assign, by = "ID") %>%
  mutate(month = month(ymd(Date), label = TRUE, abbr = FALSE),
         Region_Original_fine = categorize_location(Latitude, Longitude),
         Region_New_fine = categorize_location(Latitude_new, Longitude_new),
         dist = round(dist, 2))
nrow(finclip_final)

# Save ####
write_csv(finclip_final, here("output/fin/finclip_final.csv"))
write_csv(variance_fin, here("output/fin/variance_terms.csv"))
