
# Load required packages using pacman
pacman::p_load(raster, INLA, ggplot2, cowplot, SpatialEpi, splancs, PBSmapping, maps, mapdata, 
               dplyr, sp, sf, inlabru, lubridate, ggpolypath, units, INLAutils, beepr, viridis, 
               automap, cowplot)

# 1. Data Input
# Load environmental data for the models
load("data/isoscape_data/env_land_data/land_liz.RData")
load("data/isoscape_data/env_land_data/depth_fin.RData")
load("data/isoscape_data/env_land_data/SST_fin.RData")
load("data/isoscape_data/env_land_data/CHL_fin.RData")
land_liz <- st_set_crs(land_liz, 4326) 

# Depth Data Preparation and Plotting
depth_df <- rasterToPoints(Depth)
depth_raster_data <- data.frame(depth_df)
names(depth_raster_data) <- c("Longitude", "Latitude", "Value")
depth_plot <- ggplot() +
  geom_raster(data = depth_raster_data, aes(x = Longitude, y = Latitude, fill = Value)) +
  geom_sf(data = land_liz, fill = "black", color = "black") +
  scale_fill_viridis_c(option = "mako") +
  labs(fill = "Depth (m)") +
  theme_classic()
depth_plot

# SST Data Preparation and Plotting
SST_df <- rasterToPoints(SST)
SST_raster_data <- data.frame(SST_df)
names(SST_raster_data) <- c("Longitude", "Latitude", "Value")
SST_plot <- ggplot() +
  geom_raster(data = SST_raster_data, aes(x = Longitude, y = Latitude, fill = Value)) +
  geom_sf(data = land_liz, fill = "black", color = "black") +
  scale_fill_viridis_c(option = "inferno") +
  labs(fill = expression("Sea surface \ntemperature (" * degree * "C)")) +
  theme_classic()
SST_plot

# CHL Data Preparation and Plotting
CHL_df <- rasterToPoints(CHL)
CHL_raster_data <- data.frame(CHL_df)
names(CHL_raster_data) <- c("Longitude", "Latitude", "Value")
CHL_plot <- ggplot() +
  geom_raster(data = CHL_raster_data, aes(x = Longitude, y = Latitude, fill = Value)) +
  geom_sf(data = land_liz, fill = "black", color = "black") +
  scale_fill_viridis_c(option = "turbo") +
  labs(fill = expression("Chlorophyll " * italic(a) * " (mg/" * m^3 * ")")) +  # Correctly formatted label with only "a" italicized
  theme_classic()
CHL_plot

# Combine and Save Environmental Data Plots
combined_env_plot <- plot_grid(SST_plot, CHL_plot, depth_plot, labels = "auto", ncol = 1, align = "hv")
ggsave("figures/isoscape/envdata.png", combined_env_plot)

# Distance Calculation
Loc <- cbind(iso_liz_sf_utm$Easting, iso_liz_sf_utm$Northing)
D <- dist(Loc)
par(mfrow = c(1, 2), mar = c(5, 5, 2, 2), cex.lab = 1.5)
hist(D, freq = TRUE, main = "", xlab = "Distance between sites (km)", ylab = "Frequency")
plot(x = sort(D), y = (1:length(D)) / length(D), type = "l", xlab = "Distance between sites (km)", ylab = "Cumulative proportion")

# 2. Generate Mesh
max.edge_liz <- diff(range(st_coordinates(iso_liz_sf_utm)[, 1])) / (1 * 5)
bound.outer_liz <- diff(range(st_coordinates(iso_liz_sf_utm)[, 1])) / 10
load("data/isoscape_data/env_land_data/sea.RData") # warning code will break if R is not updated
sea_utm <- st_transform(sea, crs = "+proj=utm +zone=17 +datum=WGS84 +units=km")
sea_utm.bdry <- as(sea_utm, "Spatial") %>% inla.sp2segment()
mesh_liz <- inla.mesh.2d(boundary = sea_utm.bdry,
                         loc = st_coordinates(iso_liz_sf_utm),
                         max.edge = c(0.075, 0.5) * max.edge_liz,
                         offset = c(max.edge_liz, bound.outer_liz),
                         cutoff = 0.0025,
                         crs = st_crs(sea_utm))
mesh_liz$n
mesh_plot <- ggplot() +
  gg(mesh_liz) +
  geom_sf(data = iso_liz_sf_utm, col = 'purple', size = 1.7, alpha = 0.5) +
  ylab("Latitude") + xlab("Longitude") + theme_bw()
ggsave("figures/isoscape/mesh.png", mesh_plot)

# 3. Generate Prediction Data
load("data/isoscape_data/env_land_data/env_data.RData")
sea_sp <- as(sea, "Spatial")
env_data_pred_temp1_liz <- mask(env_data, sea_sp)
env_data_pred_liz <- na.omit(as.data.frame(env_data_pred_temp1_liz, xy = TRUE))
names(env_data_pred_liz)[names(env_data_pred_liz) == "x"] <- "Longitude"
names(env_data_pred_liz)[names(env_data_pred_liz) == "y"] <- "Latitude"
env_data_pred_liz_sf <- st_as_sf(env_data_pred_liz, coords = c("Longitude", "Latitude"), crs = 4326)
env_data_pred_liz_sf_utm <- st_transform(env_data_pred_liz_sf, crs = "+proj=utm +zone=17 +datum=WGS84 +units=km")
utm_coords <- st_coordinates(env_data_pred_liz_sf_utm)
env_data_pred_liz$X <- utm_coords[, 'X']
env_data_pred_liz$Y <- utm_coords[, 'Y']
load("data/isoscape_data/iso_liz_fin.RData")
env_data_pred_liz <- env_data_pred_liz %>%
  filter(sst_org <= max(iso_liz$sst_org) & sst_org >= min(iso_liz$sst_org)) %>%
  filter(depth_org <= max(iso_liz$depth_org) & depth_org >= min(iso_liz$depth_org)) %>%
  filter(chl_org <= max(iso_liz$chl_org) & chl_org >= min(iso_liz$chl_org))
pred_locs_liz <- cbind(env_data_pred_liz$X, env_data_pred_liz$Y)
pred_grid_liz <- inla.mesh.projector(mesh_liz, loc = pred_locs_liz)
colnames(pred_locs_liz) <- c("X", "Y")
N15_pred_pts_liz <- SpatialPoints(pred_locs_liz)

# 4. Associate Observation and Prediction Data Locations with Mesh Vertices
load("data/isoscape_data/iso_liz_sf.RData")
loc.data <- cbind(iso_liz_sf_utm$Easting, iso_liz_sf_utm$Northing)
A_b_liz <- inla.spde.make.A(mesh_liz, loc = loc.data)
spde_b_liz <- inla.spde2.matern(mesh_liz, alpha = 2)
wb.index_liz <- inla.spde.make.index('w', n.spde = spde_b_liz$n.spde)
A_pred_liz <- pred_grid_liz$proj$A
N15_obs_liz <- nrow(iso_liz)
N15_pred_liz <- nrow(pred_locs_liz)

# 4. Models
# Create the stack
N15_stack_fit_liz <- inla.stack(
  tag = "Fit",
  data = list(y = iso_liz$d15N),
  A = list(1, 1, 1, 1, A_b_liz),
  effects = list(
    Intercept = rep(1, N15_obs_liz),
    sst = iso_liz$sst,
    depth = iso_liz$depth,
    chl = iso_liz$chl,
    w = wb.index_liz))

N15_stack_pred_liz <- inla.stack(data = list(y = NA),
                                 A = list(1, 1, 1, 1, A_pred_liz), tag = "Predicted",
                                 effects = list(
                                   Intercept = rep(1, N15_pred_liz),
                                   sst = env_data_pred_liz$sst,
                                   depth = env_data_pred_liz$depth,
                                   chl = env_data_pred_liz$chl,
                                   w = wb.index_liz))

N15_stack_all_liz <- inla.stack(N15_stack_fit_liz, N15_stack_pred_liz)

# Model selection and run best model
N15_iso_liz_model_selection <- iso_liz %>%
  dplyr::select(d15N, sst, depth, chl)
names(N15_iso_liz_model_selection)
set.seed(1)
N15_spatial_f <- y ~ -1 + Intercept + sst + depth + chl + f(w, model = spde_b_liz)
set.seed(1)
N15_spatial_m_liz <- inla(N15_spatial_f, # N15_result$best_formula
                          family = "gaussian",
                          data = inla.stack.data(N15_stack_fit_liz),
                          control.compute = list(dic = TRUE),
                          control.predictor = list(A = inla.stack.A(N15_stack_fit_liz)))
# Model Outputs
summary(N15_spatial_m_liz)
autoplot(N15_spatial_m_liz, which = c(1), CI = TRUE)
N15_Fit1_liz <- N15_spatial_m_liz$summary.fitted.values$mean[1:dim(iso_liz)[1]]
N15_E1_liz <- iso_liz$d15N - N15_Fit1_liz
plot(N15_Fit1_liz, iso_liz$d15N)
cor.test(N15_Fit1_liz, iso_liz$d15N)

# Homogeneity and Normality
par(mfrow = c(1, 1), mar = c(5, 5, 2, 2), cex.lab = 1.5)
plot(x = N15_Fit1_liz, y = N15_E1_liz)
abline(h = 0, v = 0)
hist(N15_E1_liz, breaks = 25)

# 4. Predict
N15_spatial_pred_liz <- inla(N15_spatial_f,
                             family = "gaussian",
                             data = inla.stack.data(N15_stack_all_liz),
                             control.compute = list(dic = TRUE),
                             control.predictor = list(A = inla.stack.A(N15_stack_all_liz)),
                             quantiles = NULL,
                             control.inla = list(int.strategy = "eb"))

# 5. Extract Prediction Nodes and Response
N15_id_pred_liz <- inla.stack.index(N15_stack_all_liz, "Predicted")$data
N15_mean_liz <- N15_spatial_pred_liz$summary.fitted.values$mean[N15_id_pred_liz]
N15_sd_liz <- N15_spatial_pred_liz$summary.fitted.values$sd[N15_id_pred_liz]
N15_res_liz <- data.frame(pred_locs_liz, N15_mean_liz, N15_var_liz = N15_sd_liz^2)
N15_res_liz_sf <- st_as_sf(N15_res_liz, coords = c("X", "Y"), crs = "+proj=utm +zone=17 +datum=WGS84 +units=km")
N15_res_liz_sf_latlon <- st_transform(N15_res_liz_sf, crs = 4326)
N15_lat_lon <- st_coordinates(N15_res_liz_sf_latlon)
N15_res_liz$Latitude <- N15_lat_lon[, 'Y']
N15_res_liz$Longitude <- N15_lat_lon[, 'X']

# Examine the Matern Correlation Function
N15_SpFi.w <- inla.spde2.result(inla = N15_spatial_m_liz, name = "w", spde = spde_b_liz, do.transfer = TRUE)
N15_Kappa <- inla.emarginal(function(x) x, N15_SpFi.w$marginals.kappa[[1]])
N15_sigmau <- inla.emarginal(function(x) sqrt(x), N15_SpFi.w$marginals.variance.nominal[[1]])
N15_r <- inla.emarginal(function(x) x, N15_SpFi.w$marginals.range.nominal[[1]])
N15_D <- as.matrix(dist(mesh_liz$loc[, 1:2]))
N15_d.vec <- seq(0, max(N15_D), length = 100)
N15_Cor.M <- (N15_Kappa * N15_d.vec) * besselK(N15_Kappa * N15_d.vec, 1)
N15_Cor.M[1] <- 1
par(mfrow = c(1, 1), mar = c(5, 5, 2, 2))
plot(x = N15_d.vec, y = N15_Cor.M, pch = 16, type = "l", cex.lab = 1.5, xlab = "Distance (km)", ylab = "Correlation")

# Plot the Spatial Random Field
N15_w.pm <- N15_spatial_m_liz$summary.random$w$mean
N15_w.sd <- N15_spatial_m_liz$summary.random$w$sd
N15_wproj <- inla.mesh.projector(mesh_liz)
N15_w.pm100_100 <- inla.mesh.project(N15_wproj, N15_w.pm)
N15_w.sd100_100 <- inla.mesh.project(N15_wproj, N15_w.sd)
N15_Grid <- expand.grid(Xkm = N15_wproj$x, Ykm = N15_wproj$y)
N15_Grid$N15_w.pm <- as.vector(N15_w.pm100_100)
N15_Grid$N15_w.sd <- as.vector(N15_w.sd100_100)
col.l <- colorRampPalette(c('red', 'green'))(50)
col.r <- rev(rainbow(30, alpha = 0.35))
N15_plot.wpm <- lattice::levelplot(N15_w.pm ~ Xkm * Ykm, data = N15_Grid, aspect = "iso", col.regions = col.r, scales = list(draw = TRUE), xlab = list("X-coordinates (km)", cex = 1.5), ylab = list("Y-coordinates (km)", cex = 1.5), main = list("Posterior mean spatial random field", cex = 1.5))
N15_plot.wpm
lattice::levelplot(N15_w.sd ~ Xkm * Ykm, col.regions = col.l, data = N15_Grid, scales = list(draw = TRUE), xlab = list("X-coordinates (km)", cex = 1.5), ylab = list("Y-coordinates (km)", cex = 1.5), main = list("Posterior sd spatial random field", cex = 1.5))

# Plotting the GMF
N15_Grid$Xm <- N15_Grid$Xkm * 1000
N15_Grid$Ym <- N15_Grid$Ykm * 1000
N15_Grid_sf <- st_as_sf(N15_Grid, coords = c("Xm", "Ym"), crs = 32617)
st_crs(N15_Grid_sf) <- 32617
N15_Grid_geo <- st_transform(N15_Grid_sf, crs = 4326)
coords <- st_coordinates(N15_Grid_geo)
N15_Grid_geo$Longitude <- coords[, 1]
N15_Grid_geo$Latitude <- coords[, 2]
GMRF_N15 <- ggplot(data = N15_Grid_geo %>% filter(!is.na(N15_w.pm))) +
  geom_tile(aes(x = Longitude, y = Latitude, fill = N15_w.pm), width = diff(range(N15_Grid_geo$Longitude)) / sqrt(nrow(N15_Grid_geo)), height = diff(range(N15_Grid_geo$Latitude)) / sqrt(nrow(N15_Grid_geo))) +
  coord_fixed() +
  labs(fill = "Posterior mean values \nof spatial random field", x = "Longitude", y = "Latitude") +
  scale_fill_viridis_c(option = "turbo") +
  theme_minimal() +
  geom_sf(data =  st_as_sf(land_liz), fill = "darkgrey", lwd = .05) +
  theme_classic()
GMRF_N15

# Plot the results as a raster
N15_mean_liz_r <- subset(N15_res_liz, select = c("Longitude", "Latitude", "N15_mean_liz"))
coordinates(N15_mean_liz_r) <- ~Longitude + Latitude
gridded(N15_mean_liz_r) <- TRUE
N15_mean_liz_raster <- raster(N15_mean_liz_r)
plot(N15_mean_liz_raster)
N15_mean_liz_raster_df <- as.data.frame(N15_mean_liz_raster, xy = TRUE)
N15_observed_range <- range(iso_liz$d15N, na.rm = TRUE)
hist(iso_liz$d15N)
N15_predicted_range <- range(N15_mean_liz_raster_df$N15_mean_liz, na.rm = TRUE)
hist(N15_mean_liz_raster_df$N15_mean_liz)

# Plot raster using ggplot
N15_plot_INLA_liz <- ggplot() +
  geom_sf(data =  st_as_sf(land_liz), fill = "darkgrey", lwd = .05) +
  geom_raster(data = N15_mean_liz_raster_df, aes(x = x, y = y, fill = N15_mean_liz), interpolate = TRUE) +
  scale_fill_viridis_c(name = "Lizardfish INLA \nd15N", na.value = NA, option = "viridis", direction = -1) +
  geom_point(data = iso_liz, aes(x = Longitude, y = Latitude), size = 0.75) +
  theme_classic() +
  labs(x = "Longitude", y = "Latitude")

# Explore the rasters more
colfunc <- colorRampPalette(c("#a50026", "#d73027", "#f46d43", "#fdae61", "#fee090", "#ffffbf", "#e0f3f8", "#abd9e9", "#74add1", "#4575b4", "#313695"))
plot(N15_mean_liz_raster, col = colfunc(100), useRaster = FALSE)
range(iso_liz$d15N)
N15_mean_liz_raster_2 <- N15_mean_liz_raster
N15_mean_liz_raster_2[N15_mean_liz_raster_2 < min(iso_liz$d15N)] <- NA
N15_mean_liz_raster_2[N15_mean_liz_raster_2 > max(iso_liz$d15N)] <- NA
plot(N15_mean_liz_raster_2, col = colfunc(100), useRaster = FALSE)
plot(land_liz, add = TRUE, col = "white")
x_liz <- iso_liz$Longitude
y_liz <- iso_liz$Latitude
liz_xy <- cbind(x_liz,y_liz)
points(liz_xy, pch = 16, cex = 0.5, col = "black")
N15_mean_liz_raster_2[N15_mean_liz_raster_2 < min(iso_liz$depth_org)] <- NA
plot(N15_mean_liz_raster_2, col = colfunc(100), useRaster = FALSE)
plot(land_liz, add = TRUE, col = "white")
points(liz_xy, pch = 16, cex = 0.5, col = "black")

# Variance
N15_var_liz_r <- subset(N15_res_liz, select = c("Longitude", "Latitude", "N15_var_liz"))
coordinates(N15_var_liz_r) <- ~Longitude + Latitude
gridded(N15_var_liz_r) <- TRUE
N15_var_raster_liz <- raster(N15_var_liz_r)
plot(N15_var_raster_liz, col = colfunc(100), useRaster = FALSE)
points(liz_xy, pch = 16, cex = 0.5, col = "black")
plot(land_liz, add = TRUE, col = "white")

# Kriging
iso_liz_krig <- iso_liz
coordinates(iso_liz_krig) <- ~Longitude + Latitude
proj4string(iso_liz_krig) <- CRS("+proj=longlat +datum=WGS84 +no_defs")
N15_iso_liz_krig_var_auto <- autofitVariogram(d15N ~ 1, iso_liz_krig)
summary(N15_iso_liz_krig_var_auto)
plot(N15_iso_liz_krig_var_auto)
crs_projected <- "+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs"
iso_liz_krig_projected <- spTransform(iso_liz_krig, CRS(crs_projected))
env_data_pred_liz_krige <- na.omit(as.data.frame(env_data_pred_temp1_liz, xy = TRUE))
pred_locs_liz_krige <- cbind(env_data_pred_liz_krige$x, env_data_pred_liz_krige$y)
colnames(pred_locs_liz_krige) <- c("x", "y")
N15_pred_pts_liz_krige <- SpatialPoints(pred_locs_liz_krige)
proj4string(N15_pred_pts_liz_krige) <- CRS("+proj=longlat +datum=WGS84 +no_defs")
N15_pred_pts_liz_projected_krige <- spTransform(N15_pred_pts_liz_krige, CRS(crs_projected))
N15_iso_liz__krig_fin <- autoKrige(d15N ~ 1, iso_liz_krig_projected, N15_pred_pts_liz_projected_krige)
N15_iso_liz__krig_fin_sf <- st_as_sf(N15_iso_liz__krig_fin$krige_output)
range(N15_iso_liz__krig_fin_sf$var1.pred)
hist(N15_iso_liz__krig_fin_sf$var1.stdev)
N15_observed_range <- range(iso_liz$d15N, na.rm = TRUE)
N15_predicted_range <- range(N15_iso_liz__krig_fin$krige_output$var1.pred, na.rm = TRUE)
N15_plot_Krige_liz <- ggplot() +
  geom_sf(data =  st_as_sf(land_liz), fill = "darkgrey", lwd = 0.05) +
  geom_sf(data = N15_iso_liz__krig_fin_sf, aes(color = var1.pred)) +
  geom_sf(data =  st_as_sf(land_liz), fill = "darkgrey", lwd = 0.05) +
  scale_color_viridis_c(limits = N15_observed_range, option = "viridis", direction = -1) +
  theme_classic() +
  labs(x = "Longitude", y = "Latitude", color = expression(""* delta^15 * "C" * " (Krige)"))

# Combine INLA and Krige Plots
plot_grid(N15_plot_INLA_liz, N15_plot_Krige_liz)

# Kriging for d15N Model
N15_mean_liz_raster_df_reduced_temp <- N15_mean_liz_raster_df %>% 
  filter(!is.na(N15_mean_liz)) %>% 
  filter(N15_mean_liz >= range(iso_liz$d15N, na.rm = TRUE)[1] & N15_mean_liz <= range(iso_liz$d15N, na.rm = TRUE)[2])
N15_mean_liz_raster_df_reduced <- N15_mean_liz_raster_df_reduced_temp %>% sample_frac(1)
N15_iso_liz_krig_model <- N15_mean_liz_raster_df_reduced %>% 
  rename(Longitude = x, Latitude = y, d15N = N15_mean_liz)
coordinates(N15_iso_liz_krig_model) <- ~Longitude + Latitude
proj4string(N15_iso_liz_krig_model) <- CRS("+proj=longlat +datum=WGS84 +no_defs")
plot(N15_iso_liz_krig_model)
N15_iso_liz_krig_projected_model <- spTransform(N15_iso_liz_krig_model, CRS(crs_projected))
N15_iso_liz_krig_var_model_auto <- autofitVariogram(d15N ~ 1, N15_iso_liz_krig_model)
summary(N15_iso_liz_krig_var_model_auto)
plot(N15_iso_liz_krig_var_model_auto)
N15_iso_liz_krig_fin_model_sf <- autoKrige(d15N ~ 1, N15_iso_liz_krig_projected_model, N15_pred_pts_liz_projected_krige)
N15_iso_liz_krig_fin_model_sf <- st_as_sf(N15_iso_liz_krig_fin_model_sf$krige_output)
range(N15_iso_liz_krig_fin_model_sf$var1.pred)
hist(N15_iso_liz_krig_fin_model_sf$var1.stdev)
N15_plot_Krige_model_liz <- ggplot() +
  geom_sf(data =  st_as_sf(land_liz), fill = "darkgrey", lwd = 0.05) +
  geom_sf(data = N15_iso_liz_krig_fin_model_sf, aes(color = var1.pred)) +
  geom_sf(data =  st_as_sf(land_liz), fill = "darkgrey", lwd = 0.05) +
  scale_color_viridis_c(limits = N15_observed_range, option = "viridis", direction = -1) +
  theme_classic() +
  labs(x = "Longitude", y = "Latitude", color = expression(""* delta^15 * "C" * " (INLA & Krige)"))
N15_plot_Krige_model_liz_legend <- ggplot() +
  geom_sf(data =  st_as_sf(land_liz), fill = "darkgrey", lwd = 0.05) +
  geom_sf(data = N15_iso_liz_krig_fin_model_sf, aes(color = var1.pred)) +
  geom_sf(data =  st_as_sf(land_liz), fill = "darkgrey", lwd = 0.05) +
  scale_color_viridis_c(limits = N15_observed_range, option = "viridis", direction = -1) +
  theme_classic() +
  labs(x = "Longitude", y = "Latitude", color = expression(""* delta^15 * "C" * ""))
Powerpoint_N15_plot_Krige_model_liz <- ggplot() +
  geom_sf(data =  st_as_sf(land_liz), fill = "darkgrey", lwd = 0.05) +
  geom_sf(data = N15_iso_liz_krig_fin_model_sf, aes(color = var1.pred)) +
  geom_sf(data =  st_as_sf(land_liz), fill = "darkgrey", lwd = 0.05) +
  scale_color_viridis_c(limits = N15_observed_range, option = "viridis", direction = -1) +
  geom_point(data = iso_liz, aes(x = Longitude, y = Latitude), size = 0.15) +
  theme_classic() +
  labs(x = "Longitude", y = "Latitude", color = "d15N")

# Generate Final Mean Raster
N15_iso_liz_krig_fin_model_sf_transformed <- st_transform(N15_iso_liz_krig_fin_model_sf, crs = crs('+proj=longlat +datum=WGS84 +no_defs'))
coords <- st_coordinates(N15_iso_liz_krig_fin_model_sf_transformed)
values <- N15_iso_liz_krig_fin_model_sf_transformed$var1.pred
template_raster <- raster(Depth)
N15_mean_iso_liz_f <- rasterize(x = coords, y = template_raster, field = values, fun = mean)
plot(N15_mean_iso_liz_f)
save(N15_mean_iso_liz_f, file = "data/isoscape_data/N15_mean_iso_liz_f.RData")

# Extract and Plot Variance Raster
plot(N15_var_raster_liz)
N15_var_liz_raster_df <- as.data.frame(N15_var_raster_liz, xy = TRUE)
N15_var_liz_raster_df_reduced_temp <- N15_var_liz_raster_df %>% filter(!is.na(N15_var_liz))
N15_var_liz_raster_df_reduced <- N15_var_liz_raster_df_reduced_temp %>% sample_frac(1)
N15_iso_liz_krig_var_model <- N15_var_liz_raster_df_reduced %>% rename(Longitude = x, Latitude = y, N15_var = N15_var_liz)
coordinates(N15_iso_liz_krig_var_model) <- ~Longitude + Latitude
proj4string(N15_iso_liz_krig_var_model) <- CRS("+proj=longlat +datum=WGS84 +no_defs")
plot(N15_iso_liz_krig_var_model)
N15_iso_liz_krig_var_model_auto <- autofitVariogram(N15_var ~ 1, N15_iso_liz_krig_var_model)
summary(N15_iso_liz_krig_var_model_auto)
plot(N15_iso_liz_krig_var_model_auto)
crs_projected <- "+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs"
N15_iso_liz_krig_var_projected_model <- spTransform(N15_iso_liz_krig_var_model, CRS(crs_projected))
N15_iso_liz_krig_var_fin_model <- autoKrige(N15_var ~ 1, N15_iso_liz_krig_var_projected_model, N15_pred_pts_liz_projected_krige)
N15_iso_liz_krig_var_fin_model_sf <- st_as_sf(N15_iso_liz_krig_var_fin_model$krige_output)
N15_plot_Krige_var_model_liz <- ggplot() +
  geom_sf(data =  st_as_sf(land_liz), fill = "darkgrey", lwd = 0.05) +
  geom_sf(data = N15_iso_liz_krig_var_fin_model_sf, aes(color = var1.pred)) +
  geom_sf(data =  st_as_sf(land_liz), fill = "darkgrey", lwd = 0.05) +
  scale_color_viridis_c(option = "rocket", direction = -1) +
  theme_classic() +
  labs(x = "Longitude", y = "Latitude", color = expression(""* delta^15 * "C" * " var"))

# Generate Final Variance Raster
N15_iso_liz_krig_var_fin_model_sf_transformed <- st_transform(N15_iso_liz_krig_var_fin_model_sf, crs = crs('+proj=longlat +datum=WGS84 +no_defs'))
coords <- st_coordinates(N15_iso_liz_krig_var_fin_model_sf_transformed)
values <- N15_iso_liz_krig_var_fin_model_sf_transformed$var1.pred
template_raster <- raster(Depth)
N15_var_iso_liz_f <- rasterize(x = coords, y = template_raster, field = values, fun = mean)
plot(N15_var_iso_liz_f)
save(N15_var_iso_liz_f, file = "data/isoscape_data/N15_var_iso_liz_f.RData")
