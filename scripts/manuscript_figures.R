# Manuscript figures
# Figures 2-4 and S11-S14, written to figures/. S11 and S14 use a Google satellite basemap (ggmap) with the
# water labels turned off. Set GGMAP_GOOGLE_API_KEY in ~/.Renviron; without a key they are drawn on the
# coastline layer.

# Load libraries
pacman::p_load(dplyr, tidyr, readr, raster, sf, ggplot2, cowplot, ggmap, here)
source(here("scripts/source_script_assignment_functions.R"))

# Needed objects ####
load(here("data/isoscape_data/env_land_data/land_liz.RData"))
land_liz <- st_set_crs(land_liz, 4326)

finclip_final <- read_csv(here("output/fin/finclip_final.csv"), show_col_types = FALSE)
blood_final <- read_csv(here("output/blood/blood_final.csv"), show_col_types = FALSE)
sens_fin <- read_csv(here("output/sensitivity/sensitivity_fin.csv"), show_col_types = FALSE)
sens_blood <- read_csv(here("output/sensitivity/sensitivity_blood.csv"), show_col_types = FALSE)
dir.create(here("figures"), showWarnings = FALSE)

fin_raster <- function(id) raster(here("output/fin/rasters", sprintf("probability_raster_%s_C_-3.00_N_-3.50.tif", id)))
blood_raster <- function(id) raster(here("output/blood/rasters", sprintf("probability_raster_%s_C_-1.10_N_-1.40.tif", id)))
raster_df <- function(r) { df <- as.data.frame(r, xy = TRUE); names(df)[3] <- "p"; df[!is.na(df$p), ] }
plot_title <- function(s) paste0(s$ID, " (", round(s$FL_cm, digits = 0), " cm): ", as.character(as.Date(s$Date)))

legend_prob <- get_legend(
  ggplot(data.frame(x = 1:2, y = 1:2, z = c(0, 1)), aes(x, y, fill = z)) + geom_tile() +
    scale_fill_viridis_c(limits = c(0, 1), name = "Resource\nprobability") + theme_classic(base_size = 16) +
    theme(legend.key.height = unit(1.1, "cm"), legend.key.width = unit(0.7, "cm")))

# Resource probability map for one sample (red circle = capture location, arrow = highest-probability cell)
individual_map <- function(r, s, peak = NULL) {
  p <- ggplot() +
    geom_sf(data = land_liz, color = "black", fill = "white") +
    geom_tile(data = raster_df(r) %>% rename(Longitude = x, Latitude = y), aes(x = Longitude, y = Latitude, fill = p)) +
    scale_fill_viridis_c(limits = c(0, 1), option = "viridis", direction = 1) +
    coord_sf(xlim = c(-91, -75), ylim = c(24.1, 39), expand = FALSE) +
    labs(title = plot_title(s)) +
    theme_classic(base_size = 14) +
    geom_point(data = s, aes(x = Longitude, y = Latitude), color = "red", stroke = 1.5, size = 3, shape = 1) +
    guides(fill = "none")
  if (!is.null(peak)) {
    p <- p + geom_segment(data = peak, aes(x = Longitude, y = Latitude, xend = Longitude_new, yend = Latitude_new),
                          arrow = arrow(type = "closed", length = unit(0.1, "inches")), color = "black", linewidth = 0.5)
  }
  p
}

# Figure 2 ####
# Figure 2: sample locations and average resource probabilities by capture region (averaged over all fish
# sampled in the region, see script_6_7_assignment_models_output_summary.R). All panels use the 0-1 scale of the legend.
regional_map <- function(region, title) {
  r <- raster(here("output/fin/regional_averages", paste0(gsub("[ .]", "", region), ".tif")))
  ggplot() +
    geom_sf(data = land_liz, color = "black", fill = "white") +
    geom_raster(data = raster_df(r), aes(x = x, y = y, fill = p)) +
    scale_fill_viridis_c(limits = c(0, 1), option = "viridis", direction = 1) +
    coord_sf(xlim = c(-91, -75), ylim = c(24.1, 39), expand = FALSE) +
    labs(y = "Latitude", x = "Longitude", title = title) +
    theme_classic(base_size = 14) +
    guides(fill = "none")
}
set.seed(42)
plot_points <- ggplot() +
  geom_sf(data = land_liz, color = "black", fill = "white") +
  geom_jitter(data = finclip_final, aes(Longitude, Latitude), col = "red") +
  coord_sf(xlim = c(-91, -75), ylim = c(24.1, 39), expand = FALSE) +
  labs(x = "Longitude", y = "Latitude") +
  theme_classic(base_size = 14)
fig2 <- plot_grid(
  plot_grid(plot_points,
            regional_map("N. Gulf of Mexico", "Northern Gulf of Mexico"), regional_map("CW FL", "West Florida"),
            regional_map("E FL", "East Florida"), regional_map("S FL", "South Florida"), regional_map("Mid Atlantic", "Mid. Atlantic"),
            ncol = 3, labels = "auto", label_size = 24),
  legend_prob, rel_widths = c(1, 0.1))
ggsave(here("figures/Figure 2.png"), fig2, width = 16, height = 16 / 1.5601, dpi = 200, bg = "white")

# Figure 3 ####
# Figure 3: individual examples
examples <- c("T-1070", "T-673", "T-5475", "T-585", "T-349E", "T-5240")
fig3_maps <- lapply(examples, function(id) individual_map(fin_raster(id), filter(finclip_final, ID == id)))
fig3 <- plot_grid(plot_grid(plotlist = fig3_maps, ncol = 3, labels = "auto", label_size = 24), legend_prob, rel_widths = c(1, 0.12))
ggsave(here("figures/Figure 3.png"), fig3, width = 15, height = 15 / 1.3468, dpi = 200, bg = "white")

# Figure 4 ####
# Figure 4: distance to the highest-probability cell for paired fin and plasma samples, with examples
paired <- inner_join(blood_final %>% transmute(fish = sub("B$", "", ID), dist_b = dist),
                     finclip_final %>% transmute(fish = sub("[CE]$", "", ID), dist_f = dist), by = "fish") %>%
  pivot_longer(c(dist_b, dist_f), names_to = "SampleType", values_to = "Distance") %>%
  mutate(SampleType = factor(ifelse(SampleType == "dist_b", "Blood plasma", "Fin"), levels = c("Blood plasma", "Fin")))
set.seed(7)
fig4a <- ggplot(paired, aes(x = SampleType, y = Distance, fill = SampleType)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +  # every fish is drawn by geom_jitter, so outliers are not drawn twice
  geom_jitter(aes(color = SampleType), width = 0.2, size = 2, alpha = 0.8) +
  scale_fill_manual(values = c("Blood plasma" = "red", "Fin" = "blue")) +
  scale_color_manual(values = c("Blood plasma" = "darkred", "Fin" = "darkblue")) +
  theme_classic(base_size = 16) +
  labs(x = "", y = "Distance to highest assignment probability (km)") +
  theme(legend.position = "none")
fin_example <- function(id) { s <- filter(finclip_final, ID == id); individual_map(fin_raster(id), s, s) }
blood_example <- function(id) { s <- filter(blood_final, ID == id); individual_map(blood_raster(id), s, s) }
fig4_maps <- plot_grid(fin_example("T-2078E"), blood_example("T-2078B"), fin_example("T-2137E"), blood_example("T-2137B"),
                       ncol = 2, labels = c("b", "c", "d", "e"), label_size = 24)
fig4 <- plot_grid(fig4a, fig4_maps, legend_prob, nrow = 1, rel_widths = c(0.75, 1.25, 0.14), labels = c("a", "", ""), label_size = 24)
ggsave(here("figures/Figure 4.png"), fig4, width = 16, height = 16 / 1.5972, dpi = 200, bg = "white")

# Figure S11 ####
# Figure S11: capture location to highest-probability cell, by capture region
use_google <- nzchar(Sys.getenv("GGMAP_GOOGLE_API_KEY"))
if (use_google) {
  register_google(key = Sys.getenv("GGMAP_GOOGLE_API_KEY"))
  SEE <- get_googlemap(center = c(lon = -80, lat = 35), zoom = 5, size = c(640, 640), scale = 2, maptype = c("hybrid"),
                       style = "feature:water|element:labels|visibility:off")
}
to_from <- finclip_final %>%
  mutate(Region_Original_fine = factor(Region_Original_fine,
                                       levels = c("N. Gulf of Mexico", "CW FL", "S FL", "E FL", "Mid Atlantic"),
                                       labels = c("North Gulf of Mexico", "West Florida", "South Florida", "East Florida", "Mid. Atlantic")))
if (use_google) {
  figS11 <- ggmap(SEE, extent = "normal") +
    scale_x_continuous(limits = c(-92, -75)) +
    scale_y_continuous(limits = c(23, 39)) +
    ylab("Latitude") +
    xlab("Longitude") +
    geom_segment(data = to_from, aes(x = Longitude, y = Latitude, xend = Longitude_new, yend = Latitude_new),
                 arrow = arrow(length = unit(0.05, "inches")), color = "yellow", alpha = 0.25, linewidth = 0.25) +
    theme_minimal() +
    facet_wrap(~Region_Original_fine) +
    theme(strip.text = element_text(size = 10))
} else {
  figS11 <- ggplot() +
    geom_sf(data = land_liz, color = "grey40", fill = "grey88", linewidth = 0.2) +
    geom_segment(data = to_from, aes(x = Longitude, y = Latitude, xend = Longitude_new, yend = Latitude_new),
                 arrow = arrow(length = unit(0.05, "inches")), color = "darkorange3", alpha = 0.5, linewidth = 0.3) +
    coord_sf(xlim = c(-92, -75), ylim = c(23, 39), expand = FALSE) +
    facet_wrap(~Region_Original_fine) +
    labs(x = NULL, y = NULL) +
    theme_minimal() +
    theme(strip.text = element_text(size = 10), panel.background = element_rect(fill = "aliceblue", color = NA))
}
ggsave(here("figures/Figure S11.png"), figS11, width = 8, height = 8)

# Figure S12 ####
# Figure S12: original and TDF-adjusted values against the isoscape (fin and plasma)
e <- load_assignment_inputs()
overlap_plot <- function(dTC, dTN, offset_C, offset_N) {
  data_points <- data.frame(d13C = c(round(values(e$lmodC), 1), dTC, dTC + offset_C),
                            d15N = c(round(values(e$lmodN), 1), dTN, dTN + offset_N),
                            color = c(rep("grey", ncell(e$lmodC)), rep("red", length(dTC)), rep("purple", length(dTC))),
                            size = 1.5) %>%
    filter(!is.na(d13C))
  adjust_sizes <- function(data) {
    overlap_counts <- data %>% group_by(d13C, d15N) %>% summarize(count = n(), .groups = "drop")
    data <- merge(data, overlap_counts, by = c("d13C", "d15N"))
    data$size <- 1.25 * data$count
    data
  }
  combined_data <- bind_rows(filter(data_points, color == "grey"),
                             adjust_sizes(filter(data_points, color == "red")),
                             adjust_sizes(filter(data_points, color == "purple")))
  combined_data$color <- factor(combined_data$color, levels = c("grey", "red", "purple"),
                                labels = c("Isoscape values", "Original values", "Adjusted values"))
  # Legend labels give the number of overlapping samples (point size is 1.25 x count)
  k <- unique(round(pretty(range(combined_data$count, na.rm = TRUE), n = 4)))
  k <- k[k >= 1 & k <= max(combined_data$count, na.rm = TRUE)]
  ggplot() +
    geom_point(data = filter(combined_data, color == "Isoscape values"), aes(x = d13C, y = d15N, fill = color, size = size),
               shape = 21, color = "grey", stroke = 0, alpha = 0.7) +
    geom_point(data = filter(combined_data, color != "Isoscape values"), aes(x = d13C, y = d15N, fill = color, size = size),
               shape = 21, color = "black", stroke = 0.5, alpha = 0.7) +
    scale_fill_manual(values = c("Isoscape values" = "grey", "Original values" = "red", "Adjusted values" = "purple")) +
    scale_size_continuous(range = c(1.5, 5), breaks = 1.25 * k, labels = k) +
    theme_classic() +
    labs(x = expression("" * delta^13 * "C"), y = expression("" * delta^15 * "N"), fill = "Classification", size = "Overlap count") +
    theme(axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), axis.text = element_text(size = 12)) +
    guides(fill = guide_legend(order = 2), size = guide_legend(order = 1))
}
figS12 <- plot_grid(overlap_plot(finclip_final$d13C, finclip_final$d15N, -3, -3.5),
                    overlap_plot(blood_final$d13C, blood_final$d15N, -1.1, -1.4),
                    ncol = 2, labels = "auto", label_size = 28)
ggsave(here("figures/Figure S12.png"), figS12, width = 16, height = 16 / 1.9945, dpi = 200, bg = "white")

# Figure S13 ####
# Figure S13: distance across TDF combinations (fin a, plasma b).
# The 11 x 11 grid is smoothed for display by bilinear interpolation, which stays within the computed
# distances.
interp_group <- function(data, n = 100) {
  xs <- sort(unique(data$Troph.fract.C))
  ys <- sort(unique(data$Troph.fract.N))
  z <- matrix(NA_real_, length(xs), length(ys))
  z[cbind(match(data$Troph.fract.C, xs), match(data$Troph.fract.N, ys))] <- data$dist
  xo <- seq(min(xs), max(xs), length.out = n)
  yo <- seq(min(ys), max(ys), length.out = n)
  zx <- apply(z, 2, function(col) approx(xs, col, xo)$y)
  zxy <- t(apply(zx, 1, function(row) approx(ys, row, yo)$y))
  data.frame(expand.grid(X = xo, Y = yo), Z = as.vector(zxy), ID = data$ID[1])
}
sensitivity_surface <- function(sens, offset_C, offset_N) {
  df <- bind_rows(lapply(split(sens, sens$ID), interp_group)) %>%
    mutate(X = -X, Y = -Y)  # TDFs shown as positive values, as in the text
  stopifnot(min(df$Z) >= min(sens$dist), max(df$Z) <= max(sens$dist))
  df$ID <- sub("[BCE]$", "", df$ID)
  ggplot(df, aes(x = X, y = Y, fill = Z)) +
    geom_tile() +
    scale_fill_viridis_c() +
    labs(x = expression("Fractionation value (" * Delta^{13} * "C)"), y = expression("Fractionation value (" * Delta^{15} * "N)"),
         fill = "Distance \n(km)") +
    theme_classic() +
    scale_x_continuous(breaks = -offset_C + c(-0.3, 0, 0.3)) +  # ticks inside the panel, so labels of adjacent panels do not run together
    geom_vline(xintercept = -offset_C, linetype = "dashed", color = "red") +
    geom_hline(yintercept = -offset_N, linetype = "dashed", color = "red") +
    facet_wrap(~ID)
}
figS13 <- plot_grid(sensitivity_surface(sens_fin, -3, -3.5), sensitivity_surface(sens_blood, -1.1, -1.4),
                    ncol = 1, labels = "auto", label_size = 24)
ggsave(here("figures/Figure S13.png"), figS13, width = 12, height = 12 / 0.9170, dpi = 200, bg = "white")

# Figure S14 ####
# Figure S14: highest-probability cells for fin across TDF combinations (purple = 3.0 and 3.5 per mil)
sens_to_from <- sens_fin %>% left_join(dplyr::select(finclip_final, ID, Latitude, Longitude), by = "ID")
base_tdf <- filter(sens_to_from, Troph.fract.C == -3, Troph.fract.N == -3.5)
if (use_google) {
  figS14 <- ggmap(SEE, extent = "normal") +
    scale_x_continuous(limits = c(min(sens_to_from$Longitude_new) - 1.5, max(sens_to_from$Longitude_new) + 1.5)) +
    scale_y_continuous(limits = c(min(sens_to_from$Latitude_new) - 1, max(sens_to_from$Latitude_new) + 1)) +
    geom_segment(data = sens_to_from, aes(x = Longitude, y = Latitude, xend = Longitude_new, yend = Latitude_new),
                 arrow = arrow(length = unit(0.05, "inches")), color = "yellow", alpha = 0.5) +
    geom_segment(data = base_tdf, aes(x = Longitude, y = Latitude, xend = Longitude_new, yend = Latitude_new),
                 arrow = arrow(length = unit(0.1, "inches")), color = "purple", alpha = 1) +
    theme_classic() +
    theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank()) +
    facet_wrap(~ID)
} else {
  figS14 <- ggplot() +
    geom_sf(data = land_liz, color = "grey40", fill = "grey88", linewidth = 0.2) +
    geom_segment(data = sens_to_from, aes(x = Longitude, y = Latitude, xend = Longitude_new, yend = Latitude_new),
                 arrow = arrow(length = unit(0.05, "inches")), color = "darkorange2", alpha = 0.5) +
    geom_segment(data = base_tdf, aes(x = Longitude, y = Latitude, xend = Longitude_new, yend = Latitude_new),
                 arrow = arrow(length = unit(0.1, "inches")), color = "purple") +
    coord_sf(xlim = c(min(sens_to_from$Longitude_new) - 1.5, max(sens_to_from$Longitude_new) + 1.5),
             ylim = c(min(sens_to_from$Latitude_new) - 1, max(sens_to_from$Latitude_new) + 1), expand = FALSE) +
    theme_classic() +
    theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(),
          panel.background = element_rect(fill = "aliceblue", color = NA)) +
    facet_wrap(~ID)
}
ggsave(here("figures/Figure S14.png"), figS14, width = 10, height = 9)
