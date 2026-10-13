# Sensitivity of assignments to trophic discrimination factors
# Paired fin and plasma samples (n = 30), 121 TDF combinations per tissue (+/- 0.5 per mil in 0.1 steps
# around 3.0/3.5 for fin and 1.1/1.4 for plasma), with the 75% latitude prediction interval for both tissues
# and tarpon SD terms from the 30 paired samples. The fin grid is also run with the 95% interval used for the
# main fin assignments.

# Load libraries
pacman::p_load(dplyr, readr, tidyr, here)
source(here("scripts/source_script_assignment_functions.R"))

# Settings ####
CONSISTENT_KM <- 100                       # an assignment is consistent if its distance varies by less than this across the grid
THRESHOLDS_KM <- c(50, 75, 100, 150, 200)  # to show how the count depends on that choice

# Needed objects ####
e <- load_assignment_inputs()
iso <- read_csv(here("data/data_cleaning/tarpon_iso_corrected.csv"), show_col_types = FALSE)
dir.create(here("output/sensitivity"), recursive = TRUE, showWarnings = FALSE)

# TDF grids ####
# Fish with both a fin clip and a plasma sample
blood_paired <- iso %>% filter(sample_material == "plasma") %>% arrange(ID)
fin_paired <- iso %>% filter(sample_material == "fin", ID_raw %in% blood_paired$ID_raw) %>% arrange(ID)
stopifnot(nrow(fin_paired) == 30, nrow(blood_paired) == 30)

sens_fin <- assign_grid(e, fin_paired, -3, -3.5, lat_ci = "75", month_count = 5)
sens_blood <- assign_grid(e, blood_paired, -1.1, -1.4, lat_ci = "75", month_count = 3)
sens_fin_95 <- assign_grid(e, fin_paired, -3, -3.5, lat_ci = "95", month_count = 5)

write_csv(sens_fin, here("output/sensitivity/sensitivity_fin.csv"))
write_csv(sens_blood, here("output/sensitivity/sensitivity_blood.csv"))
write_csv(sens_fin_95, here("output/sensitivity/sensitivity_fin_95.csv"))

# Consistent assignments ####
# An assignment is consistent if its distance to the highest-probability cell varies by less than
# CONSISTENT_KM across all 121 combinations. The rule measures distance, so consistent_peak_moved counts
# consistent samples whose highest-probability cell moves more than CONSISTENT_KM from the cell at the base
# TDFs (to another cell at a similar distance).
grids <- bind_rows(sens_fin %>% mutate(tissue = "fin", lat_interval = "75%"),
                   sens_blood %>% mutate(tissue = "plasma", lat_interval = "75%"),
                   sens_fin_95 %>% mutate(tissue = "fin", lat_interval = "95%")) %>%
  mutate(base = ifelse(tissue == "fin", Troph.fract.C == -3 & Troph.fract.N == -3.5, Troph.fract.C == -1.1 & Troph.fract.N == -1.4))
stopifnot(sum(grids$base) == 90)
dist_ranges <- grids %>%
  group_by(tissue, lat_interval, ID) %>%
  summarise(dist_range = max(dist) - min(dist),
            peak_shift = max(gc_dist_km(Longitude_new[base], Latitude_new[base], Longitude_new, Latitude_new)), .groups = "drop")
stability <- dist_ranges %>%
  group_by(tissue, lat_interval) %>%
  summarise(n = n(), consistent = sum(dist_range < CONSISTENT_KM),
            consistent_peak_moved = sum(dist_range < CONSISTENT_KM & peak_shift > CONSISTENT_KM), .groups = "drop")
stability
stability_by_threshold <- expand_grid(distinct(dist_ranges, tissue, lat_interval), threshold_km = THRESHOLDS_KM) %>%
  rowwise() %>%
  mutate(consistent = sum(dist_ranges$dist_range[dist_ranges$tissue == tissue & dist_ranges$lat_interval == lat_interval] < threshold_km)) %>%
  ungroup()
stability_by_threshold
write_csv(stability, here("output/sensitivity/stability_summary.csv"))
write_csv(stability_by_threshold, here("output/sensitivity/stability_by_threshold.csv"))
