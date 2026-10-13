# Assignment model output summary
# Values reported in the Results, the regional averages in Figure 2, and the fit of the fin values to the
# isoscape (Supplementary Materials, Isotopic discrimination and sensitivity).

# Load libraries
pacman::p_load(dplyr, tidyr, readr, raster, here)
source(here("scripts/source_script_assignment_functions.R"))

# Load assignment outputs ####
finclip_final <- read_csv(here("output/fin/finclip_final.csv"), show_col_types = FALSE)
blood_final <- read_csv(here("output/blood/blood_final.csv"), show_col_types = FALSE)
dir.create(here("output/summary"), recursive = TRUE, showWarnings = FALSE)

msd <- function(x, d = 1) sprintf(paste0("%.", d, "f +/- %.", d, "f"), mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
rng <- function(x, d = 1) sprintf(paste0("%.", d, "f to %.", d, "f"), min(x, na.rm = TRUE), max(x, na.rm = TRUE))

# Results values ####
# Isotope values and fork length (cm)
iso_summary <- bind_rows(
  finclip_final %>% summarise(tissue = "fin", n = n(), d13C_mean_sd = msd(d13C), d13C_range = rng(d13C), d15N_mean_sd = msd(d15N), d15N_range = rng(d15N),
                              FL_mean_sd = msd(FL_cm, 2), FL_range = rng(FL_cm)),
  blood_final %>% summarise(tissue = "plasma", n = n(), d13C_mean_sd = msd(d13C), d13C_range = rng(d13C), d15N_mean_sd = msd(d15N), d15N_range = rng(d15N),
                            FL_mean_sd = msd(FL_cm, 2), FL_range = rng(FL_cm)))
iso_summary

# Distance to the highest-probability cell
dist_summary <- bind_rows(
  finclip_final %>% summarise(tissue = "fin", n = n(), dist_mean_sd = msd(dist), dist_range = rng(dist)),
  blood_final %>% summarise(tissue = "plasma", n = n(), dist_mean_sd = msd(dist), dist_range = rng(dist)))
dist_summary

# Paired fin and plasma samples
paired <- inner_join(blood_final %>% transmute(fish = sub("B$", "", ID), dist_b = dist),
                     finclip_final %>% transmute(fish = sub("[CE]$", "", ID), dist_f = dist), by = "fish")
paired_test <- t.test(paired$dist_b, paired$dist_f, paired = TRUE)
paired_summary <- data.frame(n = nrow(paired), plasma = msd(paired$dist_b), fin = msd(paired$dist_f),
                             t = round(unname(paired_test$statistic), 2), df = unname(paired_test$parameter),
                             p = signif(paired_test$p.value, 2))
paired_summary

# Contingent assignments and region of the highest-probability cell, by capture region
table(finclip_final$assigned_polygon)
region_table <- finclip_final %>% count(Region_Original_fine, Region_New_fine) %>%
  pivot_wider(names_from = Region_New_fine, values_from = n, values_fill = 0)
region_table

# Latitude of the highest-probability cells by capture region and peak region
peak_latitudes <- finclip_final %>%
  group_by(Region_Original_fine, Region_New_fine) %>%
  summarise(n = n(), lat_min = min(Latitude_new), lat_median = median(Latitude_new), lat_max = max(Latitude_new), .groups = "drop")
peak_latitudes

# Regional averages (Figure 2) ####
# Regional averages of the resource probability surfaces (Figure 2). Each cell is averaged over all fish
# sampled in the region, with a cell outside a fish's search area counting as 0 for that fish. Cells outside
# every fish's search area are left empty.
regions <- c("N. Gulf of Mexico", "CW FL", "E FL", "S FL", "Mid Atlantic")
region_stack <- function(region) {
  ids <- finclip_final$ID[finclip_final$Region_Original_fine == region]
  stack(here("output/fin/rasters", sprintf("probability_raster_%s_C_-3.00_N_-3.50.tif", ids)))
}
stacks <- lapply(regions, region_stack)
names(stacks) <- regions
search_area_list <- lapply(stacks, function(s) calc(s, mean, na.rm = TRUE))  # empty outside every search area
avg_list <- lapply(regions, function(r) mask(calc(stacks[[r]], function(x) sum(x, na.rm = TRUE) / length(x)), search_area_list[[r]]))
names(avg_list) <- regions
n_cover_list <- lapply(stacks, function(s) calc(s, function(x) sum(!is.na(x))))

dir.create(here("output/fin/regional_averages"), showWarnings = FALSE)
for (r in regions) writeRaster(avg_list[[r]], here("output/fin/regional_averages", paste0(gsub("[ .]", "", r), ".tif")), overwrite = TRUE)

xy <- xyFromCell(avg_list[[1]], 1:ncell(avg_list[[1]]))
cell_region <- categorize_location(xy[, 2], xy[, 1])

# Mean and maximum of each surface within each region
regional_probability <- bind_rows(lapply(regions, function(r) {
  data.frame(capture_region = r, cell_region = cell_region, p = getValues(avg_list[[r]])) %>%
    filter(!is.na(p), cell_region != "Unlabeled") %>%
    group_by(capture_region, cell_region) %>%
    summarise(mean_probability = round(mean(p), 2), max_probability = round(max(p), 2), .groups = "drop")
}))
regional_probability

# Where each surface peaks, the latitude span of its highest cells, and how many fish reach each cell
surface_summary <- bind_rows(lapply(regions, function(r) {
  p <- getValues(avg_list[[r]]); n_cover <- getValues(n_cover_list[[r]])
  top <- which.max(p)
  top_decile <- which(p >= quantile(p, 0.9, na.rm = TRUE))
  data.frame(capture_region = r, n_fish = nlayers(stacks[[r]]),
             p_min = min(p, na.rm = TRUE), p_max = max(p, na.rm = TRUE),
             peak_lat = xy[top, 2], peak_lon = xy[top, 1], peak_region = cell_region[top],
             top_decile_lat_min = min(xy[top_decile, 2]), top_decile_lat_max = max(xy[top_decile, 2]),
             median_fish_per_cell = median(n_cover[!is.na(p)]))
}))
surface_summary

# Fit to the isoscape ####
# TDF-adjusted fin d13C (TDF 3.0 per mil) against the lowest d13C of the 0-20 m isoscape
e <- load_assignment_inputs()
depth_mask <- e$depth_raster <= 0 & e$depth_raster >= -20
depth_mask[depth_mask == 0] <- NA
iso_c_min <- min(getValues(mask(e$lmodC, depth_mask)), na.rm = TRUE)
iso_fit <- data.frame(n = nrow(finclip_final), isoscape_d13C_min = round(iso_c_min, 2),
                      below_isoscape = sum(finclip_final$d13C - 3 < iso_c_min))
iso_fit

# Save ####
write_csv(iso_summary, here("output/summary/isotope_summary.csv"))
write_csv(dist_summary, here("output/summary/distance_summary.csv"))
write_csv(paired_summary, here("output/summary/fin_vs_plasma.csv"))
write_csv(region_table, here("output/summary/peak_region_by_capture_region.csv"))
write_csv(peak_latitudes, here("output/summary/peak_latitudes.csv"))
write_csv(regional_probability, here("output/summary/regional_probability.csv"))
write_csv(surface_summary, here("output/summary/regional_surfaces.csv"))
write_csv(iso_fit, here("output/summary/isoscape_fit.csv"))
