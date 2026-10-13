# Ethanol adjustment of fin isotope values
# Most fin clips were stored in ethanol. Paired ethanol-stored and frozen clips (45 fish) are compared with a
# paired t-test, and ethanol-stored values are shifted by the mean paired difference when it is significant
# (d13C; d15N did not differ and is used as measured). Reads data/data_cleaning/tarpon_iso.csv and writes the
# adjusted values to data/data_cleaning/tarpon_iso_corrected.csv and the tests to output/ethanol.

# Load libraries
pacman::p_load(dplyr, tidyr, readr, here)

ALPHA <- 0.05  # offset applied only when the paired difference is significant

# Load our isotope data ####
tarpon_iso <- read_csv(here("data/data_cleaning/tarpon_iso.csv"), show_col_types = FALSE)
dir.create(here("output/ethanol"), recursive = TRUE, showWarnings = FALSE)

# Paired ethanol and frozen clips ####
# Pivot the data to one row per fish with the ethanol and frozen values for one isotope
pair_values <- function(data, iso) {
  data %>%
    dplyr::select(ID_raw, Type, value = all_of(iso)) %>%
    filter(Type %in% c("ethanol", "no ethanol")) %>%
    pivot_wider(names_from = Type, values_from = value) %>%
    rename(frozen = `no ethanol`) %>%
    filter(!is.na(ethanol), !is.na(frozen))
}

# Paired t-test of ethanol minus frozen (Wilcoxon signed-rank test as a check against outliers);
# the offset is what gets added to ethanol-stored values
paired_test <- function(pairs, iso) {
  d <- pairs$ethanol - pairs$frozen
  tt <- t.test(d)
  data.frame(isotope = iso, n = length(d), mean_diff_ethanol_minus_frozen = mean(d), sd_diff = sd(d),
             t = unname(tt$statistic), df = unname(tt$parameter), p = tt$p.value,
             wilcoxon_p = wilcox.test(d, exact = FALSE)$p.value,
             offset_applied = ifelse(tt$p.value < ALPHA, -mean(d), 0))
}

offsets <- bind_rows(lapply(c("d13C", "d15N"), function(iso) paired_test(pair_values(tarpon_iso, iso), iso)))
offsets

# T-2078 ####
# The ethanol clip of T-2078 is left out of the paired comparison (d13C 12.3 per mil below its frozen clip).
# Same tests with that pair included, for reference.
t2078 <- read_csv(here("data/data_cleaning/ethanol_outlier_T-2078.csv"), show_col_types = FALSE) %>%
  dplyr::select(ID, ID_raw, Type, d13C, d15N, CN_ratio_original)
offsets_with_t2078 <- bind_rows(lapply(c("d13C", "d15N"), function(iso) paired_test(pair_values(bind_rows(tarpon_iso, t2078), iso), iso)))
offsets_with_t2078

write_csv(offsets, here("output/ethanol/paired_tests.csv"))
write_csv(offsets_with_t2078, here("output/ethanol/paired_tests_with_T-2078.csv"))

# Adjust ethanol-stored values ####
# Shift ethanol-stored samples by the offsets (rounded once, to 0.1 per mil as reported)
for (iso in c("d13C", "d15N")) {
  off <- offsets$offset_applied[offsets$isotope == iso]
  raw <- tarpon_iso[[iso]]
  tarpon_iso[[paste0(iso, "_uncorrected")]] <- round(raw, 1)
  tarpon_iso[[paste0(iso, "_corrected")]] <- ifelse(tarpon_iso$Type == "ethanol", round(raw + off, 1), NA)
  tarpon_iso[[iso]] <- ifelse(tarpon_iso$Type == "ethanol", tarpon_iso[[paste0(iso, "_corrected")]], round(raw, 1))
}

# Remove ethanol rows if a frozen clip from the same fish is available
tarpon_iso <- tarpon_iso %>%
  group_by(ID_raw) %>%
  filter(!(Type == "ethanol" & any(Type == "no ethanol") & !all(Type == "no ethanol")) | Type == "plasma") %>%
  ungroup() %>%
  mutate(across(c(C, N, CN_ratio_original), ~ round(.x, 1)),
         sample_material = ifelse(Type != "plasma", "fin", Type))

table(tarpon_iso$sample_material)
stopifnot(sum(tarpon_iso$sample_material == "fin") == 417, sum(tarpon_iso$sample_material == "plasma") == 30)

# Save ####
write_csv(tarpon_iso, here("data/data_cleaning/tarpon_iso_corrected.csv"))
