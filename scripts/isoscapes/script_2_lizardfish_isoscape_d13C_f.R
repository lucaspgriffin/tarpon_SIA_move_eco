# Load required packages using pacman
pacman::p_load(raster, INLA, ggplot2, cowplot, SpatialEpi, splancs, PBSmapping, maps, mapdata, 
               dplyr, sp, sf, inlabru, lubridate, ggpolypath, units, INLAutils, beepr, magma, 
               automap, cowplot)

# 1. Data Input
# Load environmental data for the models
load("data/isoscape_data/env_land_data/land_liz.RData")
load("data/isoscape_data/env_land_data/depth_fin.RData")
load("data/isoscape_data/env_land_data/SST_fin.RData")
load("data/isoscape_data/env_land_data/CHL_fin.RData")
load("data/isoscape_data/env_land_data/iso_liz_sf_utm.RData")
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
  labs(fill = expression("Chlorophyll a (mg/" * m^3 * ")")) +
  theme_classic()
CHL_plot

# Combine and Save Environmental Data Plots
combined_env_plot <- plot_grid(SST_plot, CHL_plot, depth_plot, labels = "auto", ncol = 1, align = "hv")

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
C13_pred_pts_liz <- SpatialPoints(pred_locs_liz)

# Summary of env covariates observed
range(iso_liz$sst_org)
range(iso_liz$chl_org)
range(iso_liz$depth_org)
iso_liz %>% filter(depth_org > 0)

# 4. Associate Observation and Prediction Data Locations with Mesh Vertices
load("data/isoscape_data/iso_liz_sf.RData")
head(iso_liz_sf_utm)
loc.data <- cbind(iso_liz_sf_utm$Easting, iso_liz_sf_utm$Northing)
A_b_liz <- inla.spde.make.A(mesh_liz, loc = loc.data)
spde_b_liz <- inla.spde2.matern(mesh_liz, alpha = 2)
wb.index_liz <- inla.spde.make.index('w', n.spde = spde_b_liz$n.spde)
A_pred_liz <- pred_grid_liz$proj$A
C13_obs_liz <- nrow(iso_liz)
C13_pred_liz <- nrow(pred_locs_liz)

# 4. Models
# Create the stack
C13_stack_fit_liz <- inla.stack(
  tag = "Fit",
  data = list(y = iso_liz$d13C),
  A = list(1, 1, 1, 1, A_b_liz),
  effects = list(
    Intercept = rep(1, C13_obs_liz),
    sst = iso_liz$sst,
    depth = iso_liz$depth,
    chl = iso_liz$chl,
    w = wb.index_liz))

C13_stack_pred_liz <- inla.stack(data = list(y = NA),
                                 A = list(1, 1, 1, 1, A_pred_liz), tag = "Predicted",
                                 effects = list(
                                   Intercept = rep(1, C13_pred_liz),
                                   sst = env_data_pred_liz$sst,
                                   depth = env_data_pred_liz$depth,
                                   chl = env_data_pred_liz$chl,
                                   w = wb.index_liz))

C13_stack_all_liz <- inla.stack(C13_stack_fit_liz, C13_stack_pred_liz)

# Model selection and run best model
C13_iso_liz_model_selection <- iso_liz %>%
  dplyr::select(d13C, sst, depth, chl)
names(C13_iso_liz_model_selection)
set.seed(1)
C13_spatial_f <- y ~ -1 + Intercept + sst + depth + chl + f(w, model = spde_b_liz)
set.seed(1)
C13_spatial_m_liz <- inla(C13_spatial_f, # C13_result$best_formula
                          family = "gaussian",
                          data = inla.stack.data(C13_stack_fit_liz),
                          control.compute = list(dic = TRUE),
                          control.predictor = list(A = inla.stack.A(C13_stack_fit_liz)))
# Model Outputs
summary(C13_spatial_m_liz)
autoplot(C13_spatial_m_liz, which = c(1), CI = TRUE)
C13_Fit1_liz <- C13_spatial_m_liz$summary.fitted.values$mean[1:dim(iso_liz)[1]]
C13_E1_liz <- iso_liz$d13C - C13_Fit1_liz
plot(C13_Fit1_liz, iso_liz$d13C)
cor.test(C13_Fit1_liz, iso_liz$d13C)

# Homogeneity and Normality
par(mfrow = c(1, 1), mar = c(5, 5, 2, 2), cex.lab = 1.5)
plot(x = C13_Fit1_liz, y = C13_E1_liz)
abline(h = 0, v = 0)
hist(C13_E1_liz, breaks = 25)

# 4. Predict
C13_spatial_pred_liz <- inla(C13_spatial_f,
                             family = "gaussian",
                             data = inla.stack.data(C13_stack_all_liz),
                             control.compute = list(dic = TRUE),
                             control.predictor = list(A = inla.stack.A(C13_stack_all_liz)),
                             quantiles = NULL,
                             control.inla = list(int.strategy = "eb"))

# 5. Extract Prediction Nodes and Response
C13_id_pred_liz <- inla.stack.index(C13_stack_all_liz, "Predicted")$data
C13_mean_liz <- C13_spatial_pred_liz$summary.fitted.values$mean[C13_id_pred_liz]
C13_sd_liz <- C13_spatial_pred_liz$summary.fitted.values$sd[C13_id_pred_liz]
C13_res_liz <- data.frame(pred_locs_liz, C13_mean_liz, C13_var_liz = C13_sd_liz^2)
C13_res_liz_sf <- st_as_sf(C13_res_liz, coords = c("X", "Y"), crs = "+proj=utm +zone=17 +datum=WGS84 +units=km")
C13_res_liz_sf_latlon <- st_transform(C13_res_liz_sf, crs = 4326)
C13_lat_lon <- st_coordinates(C13_res_liz_sf_latlon)
C13_res_liz$Latitude <- C13_lat_lon[, 'Y']
C13_res_liz$Longitude <- C13_lat_lon[, 'X']

# Examine the Matern Correlation Function
C13_SpFi.w <- inla.spde2.result(inla = C13_spatial_m_liz, name = "w", spde = spde_b_liz, do.transfer = TRUE)
C13_Kappa <- inla.emarginal(function(x) x, C13_SpFi.w$marginals.kappa[[1]])
C13_sigmau <- inla.emarginal(function(x) sqrt(x), C13_SpFi.w$marginals.variance.nominal[[1]])
C13_r <- inla.emarginal(function(x) x, C13_SpFi.w$marginals.range.nominal[[1]])
C13_D <- as.matrix(dist(mesh_liz$loc[, 1:2]))
C13_d.vec <- seq(0, max(C13_D), length = 100)
C13_Cor.M <- (C13_Kappa * C13_d.vec) * besselK(C13_Kappa * C13_d.vec, 1)
C13_Cor.M[1] <- 1
par(mfrow = c(1, 1), mar = c(5, 5, 2, 2))
plot(x = C13_d.vec, y = C13_Cor.M, pch = 16, type = "l", cex.lab = 1.5, xlab = "Distance (km)", ylab = "Correlation")

# Plot the Spatial Random Field
C13_w.pm <- C13_spatial_m_liz$summary.random$w$mean
C13_w.sd <- C13_spatial_m_liz$summary.random$w$sd
C13_wproj <- inla.mesh.projector(mesh_liz)
C13_w.pm100_100 <- inla.mesh.project(C13_wproj, C13_w.pm)
C13_w.sd100_100 <- inla.mesh.project(C13_wproj, C13_w.sd)
C13_Grid <- expand.grid(Xkm = C13_wproj$x, Ykm = C13_wproj$y)
C13_Grid$C13_w.pm <- as.vector(C13_w.pm100_100)
C13_Grid$C13_w.sd <- as.vector(C13_w.sd100_100)
col.l <- colorRampPalette(c('red', 'green'))(50)
col.r <- rev(rainbow(30, alpha = 0.35))
C13_plot.wpm <- lattice::levelplot(C13_w.pm ~ Xkm * Ykm, data = C13_Grid, aspect = "iso", col.regions = col.r, scales = list(draw = TRUE), xlab = list("X-coordinates (km)", cex = 1.5), ylab = list("Y-coordinates (km)", cex = 1.5), main = list("Posterior mean spatial random field", cex = 1.5))
C13_plot.wpm
lattice::levelplot(C13_w.sd ~ Xkm * Ykm, col.regions = col.l, data = C13_Grid, scales = list(draw = TRUE), xlab = list("X-coordinates (km)", cex = 1.5), ylab = list("Y-coordinates (km)", cex = 1.5), main = list("Posterior sd spatial random field", cex = 1.5))

# Plotting the GMF
C13_Grid$Xm <- C13_Grid$Xkm * 1000
C13_Grid$Ym <- C13_Grid$Ykm * 1000
C13_Grid_sf <- st_as_sf(C13_Grid, coords = c("Xm", "Ym"), crs = 32617)
st_crs(C13_Grid_sf) <- 32617
C13_Grid_geo <- st_transform(C13_Grid_sf, crs = 4326)
coords <- st_coordinates(C13_Grid_geo)
C13_Grid_geo$Longitude <- coords[, 1]
C13_Grid_geo$Latitude <- coords[, 2]
GMRF_C13 <- ggplot(data = C13_Grid_geo %>% filter(!is.na(C13_w.pm))) +
  geom_tile(aes(x = Longitude, y = Latitude, fill = C13_w.pm), width = diff(range(C13_Grid_geo$Longitude)) / sqrt(nrow(C13_Grid_geo)), height = diff(range(C13_Grid_geo$Latitude)) / sqrt(nrow(C13_Grid_geo))) +
  coord_fixed() +
  labs(fill = "Posterior mean values \nof spatial random field", x = "Longitude", y = "Latitude") +
  scale_fill_viridis_c(option = "turbo") +
  theme_minimal() +
  geom_sf(data =  st_as_sf(land_liz), fill = "darkgrey", lwd = .05) +
  theme_classic()
GMRF_C13

# Plot the results as a raster
C13_mean_liz_r <- subset(C13_res_liz, select = c("Longitude", "Latitude", "C13_mean_liz"))
coordinates(C13_mean_liz_r) <- ~Longitude + Latitude
gridded(C13_mean_liz_r) <- TRUE
C13_mean_liz_raster <- raster(C13_mean_liz_r)
plot(C13_mean_liz_raster)
C13_mean_liz_raster_df <- as.data.frame(C13_mean_liz_raster, xy = TRUE)
C13_observed_range <- range(iso_liz$d13C, na.rm = TRUE)
hist(iso_liz$d13C)
C13_predicted_range <- range(C13_mean_liz_raster_df$C13_mean_liz, na.rm = TRUE)
hist(C13_mean_liz_raster_df$C13_mean_liz)

# Plot raster using ggplot
C13_plot_INLA_liz <- ggplot() +
  geom_sf(data =  st_as_sf(land_liz), fill = "darkgrey", lwd = .05) +
  geom_raster(data = C13_mean_liz_raster_df, aes(x = x, y = y, fill = C13_mean_liz), interpolate = TRUE) +
  scale_fill_viridis_c(name = "Lizardfish INLA \nd13C", na.value = NA, option = "magma", direction = -1) +
  geom_point(data = iso_liz, aes(x = Longitude, y = Latitude), size = 0.75) +
  theme_classic() +
  labs(x = "Longitude", y = "Latitude")

# Explore the rasters more
colfunc <- colorRampPalette(c("#a50026", "#d73027", "#f46d43", "#fdae61", "#fee090", "#ffffbf", "#e0f3f8", "#abd9e9", "#74add1", "#4575b4", "#313695"))
plot(C13_mean_liz_raster, col = colfunc(100), useRaster = FALSE)
range(iso_liz$d13C)
C13_mean_liz_raster_2 <- C13_mean_liz_raster
C13_mean_liz_raster_2[C13_mean_liz_raster_2 < min(iso_liz$d13C)] <- NA
C13_mean_liz_raster_2[C13_mean_liz_raster_2 > max(iso_liz$d13C)] <- NA
plot(C13_mean_liz_raster_2, col = colfunc(100), useRaster = FALSE)
plot(land_liz, add = TRUE, col = "white")
points(liz_xy, pch = 16, cex = 0.5, col = "black")
C13_mean_liz_raster_2[C13_mean_liz_raster_2 < min(iso_liz$depth_org)] <- NA
plot(C13_mean_liz_raster_2, col = colfunc(100), useRaster = FALSE)
plot(land_liz, add = TRUE, col = "white")
points(liz_xy, pch = 16, cex = 0.5, col = "black")

# Variance
C13_var_liz_r <- subset(C13_res_liz, select = c("Longitude", "Latitude", "C13_var_liz"))
coordinates(C13_var_liz_r) <- ~Longitude + Latitude
gridded(C13_var_liz_r) <- TRUE
C13_var_raster_liz <- raster(C13_var_liz_r)
plot(C13_var_raster_liz, col = colfunc(100), useRaster = FALSE)
x_liz <- iso_liz$Longitude
y_liz <- iso_liz$Latitude
liz_xy <- cbind(x_liz,y_liz)
points(liz_xy, pch = 16, cex = 0.5, col = "black")
plot(land_liz, add = TRUE, col = "white")

# Kriging
iso_liz_krig <- iso_liz
coordinates(iso_liz_krig) <- ~Longitude + Latitude
proj4string(iso_liz_krig) <- CRS("+proj=longlat +datum=WGS84 +no_defs")
C13_iso_liz_krig_var_auto <- autofitVariogram(d13C ~ 1, iso_liz_krig)
summary(C13_iso_liz_krig_var_auto)
plot(C13_iso_liz_krig_var_auto)
crs_projected <- "+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs"
iso_liz_krig_projected <- spTransform(iso_liz_krig, CRS(crs_projected))
env_data_pred_liz_krige <- na.omit(as.data.frame(env_data_pred_temp1_liz, xy = TRUE))
pred_locs_liz_krige <- cbind(env_data_pred_liz_krige$x, env_data_pred_liz_krige$y)
colnames(pred_locs_liz_krige) <- c("x", "y")
C13_pred_pts_liz_krige <- SpatialPoints(pred_locs_liz_krige)
proj4string(C13_pred_pts_liz_krige) <- CRS("+proj=longlat +datum=WGS84 +no_defs")
C13_pred_pts_liz_projected_krige <- spTransform(C13_pred_pts_liz_krige, CRS(crs_projected))
C13_iso_liz__krig_fin <- autoKrige(d13C ~ 1, iso_liz_krig_projected, C13_pred_pts_liz_projected_krige)
C13_iso_liz__krig_fin_sf <- st_as_sf(C13_iso_liz__krig_fin$krige_output)
range(C13_iso_liz__krig_fin_sf$var1.pred)
hist(C13_iso_liz__krig_fin_sf$var1.stdev)
C13_observed_range <- range(iso_liz$d13C, na.rm = TRUE)
C13_predicted_range <- range(C13_iso_liz__krig_fin$krige_output$var1.pred, na.rm = TRUE)
C13_plot_Krige_liz <- ggplot() +
  geom_sf(data =  st_as_sf(land_liz), fill = "darkgrey", lwd = 0.05) +
  geom_sf(data = C13_iso_liz__krig_fin_sf, aes(color = var1.pred)) +
  geom_sf(data =  st_as_sf(land_liz), fill = "darkgrey", lwd = 0.05) +
  scale_color_viridis_c(limits = C13_observed_range, option = "magma", direction = -1) +
  theme_classic() +
  labs(x = "Longitude", y = "Latitude", color = expression(""* delta^13 * "C" * " (Krige)"))

# Combine INLA and Krige Plots
plot_grid(C13_plot_INLA_liz, C13_plot_Krige_liz)

# Kriging for d13C Model
C13_mean_liz_raster_df_reduced_temp <- C13_mean_liz_raster_df %>% 
  filter(!is.na(C13_mean_liz)) %>% 
  filter(C13_mean_liz >= range(iso_liz$d13C, na.rm = TRUE)[1] & C13_mean_liz <= range(iso_liz$d13C, na.rm = TRUE)[2])
C13_mean_liz_raster_df_reduced <- C13_mean_liz_raster_df_reduced_temp %>% sample_frac(1)
C13_iso_liz_krig_model <- C13_mean_liz_raster_df_reduced %>% 
  rename(Longitude = x, Latitude = y, d13C = C13_mean_liz)
coordinates(C13_iso_liz_krig_model) <- ~Longitude + Latitude
proj4string(C13_iso_liz_krig_model) <- CRS("+proj=longlat +datum=WGS84 +no_defs")
plot(C13_iso_liz_krig_model)
C13_iso_liz_krig_projected_model <- spTransform(C13_iso_liz_krig_model, CRS(crs_projected))
C13_iso_liz_krig_var_model_auto <- autofitVariogram(d13C ~ 1, C13_iso_liz_krig_model)
summary(C13_iso_liz_krig_var_model_auto)
plot(C13_iso_liz_krig_var_model_auto)
C13_iso_liz_krig_fin_model_sf <- autoKrige(d13C ~ 1, C13_iso_liz_krig_projected_model, C13_pred_pts_liz_projected_krige)
C13_iso_liz_krig_fin_model_sf <- st_as_sf(C13_iso_liz_krig_fin_model_sf$krige_output)
range(C13_iso_liz_krig_fin_model_sf$var1.pred)
hist(C13_iso_liz_krig_fin_model_sf$var1.stdev)
C13_plot_Krige_model_liz <- ggplot() +
  geom_sf(data =  st_as_sf(land_liz), fill = "darkgrey", lwd = 0.05) +
  geom_sf(data = C13_iso_liz_krig_fin_model_sf, aes(color = var1.pred)) +
  geom_sf(data =  st_as_sf(land_liz), fill = "darkgrey", lwd = 0.05) +
  scale_color_viridis_c(limits = C13_observed_range, option = "magma", direction = -1) +
  theme_classic() +
  labs(x = "Longitude", y = "Latitude", color = expression(""* delta^13 * "C" * " (INLA & Krige)"))
C13_plot_Krige_model_liz_legend <- ggplot() +
  geom_sf(data =  st_as_sf(land_liz), fill = "darkgrey", lwd = 0.05) +
  geom_sf(data = C13_iso_liz_krig_fin_model_sf, aes(color = var1.pred)) +
  geom_sf(data =  st_as_sf(land_liz), fill = "darkgrey", lwd = 0.05) +
  scale_color_viridis_c(limits = C13_observed_range, option = "magma", direction = -1) +
  theme_classic() +
  labs(x = "Longitude", y = "Latitude", color = expression(""* delta^13 * "C" * ""))
Powerpoint_C13_plot_Krige_model_liz <- ggplot() +
  geom_sf(data =  st_as_sf(land_liz), fill = "darkgrey", lwd = 0.05) +
  geom_sf(data = C13_iso_liz_krig_fin_model_sf, aes(color = var1.pred)) +
  geom_sf(data =  st_as_sf(land_liz), fill = "darkgrey", lwd = 0.05) +
  scale_color_viridis_c(limits = C13_observed_range, option = "magma", direction = -1) +
  geom_point(data = iso_liz, aes(x = Longitude, y = Latitude), size = 0.15) +
  theme_classic() +
  labs(x = "Longitude", y = "Latitude", color = "d13C")

# Generate Final Mean Raster
C13_iso_liz_krig_fin_model_sf_transformed <- st_transform(C13_iso_liz_krig_fin_model_sf, crs = crs('+proj=longlat +datum=WGS84 +no_defs'))
coords <- st_coordinates(C13_iso_liz_krig_fin_model_sf_transformed)
values <- C13_iso_liz_krig_fin_model_sf_transformed$var1.pred
template_raster <- raster(Depth)
C13_mean_iso_liz_f <- rasterize(x = coords, y = template_raster, field = values, fun = mean)
plot(C13_mean_iso_liz_f)
save(C13_mean_iso_liz_f, file = "data/isoscape_data/C13_mean_iso_liz_f.RData")

# Extract and Plot Variance Raster
plot(C13_var_raster_liz)
C13_var_liz_raster_df <- as.data.frame(C13_var_raster_liz, xy = TRUE)
C13_var_liz_raster_df_reduced_temp <- C13_var_liz_raster_df %>% filter(!is.na(C13_var_liz))
C13_var_liz_raster_df_reduced <- C13_var_liz_raster_df_reduced_temp %>% sample_frac(1)
C13_iso_liz_krig_var_model <- C13_var_liz_raster_df_reduced %>% rename(Longitude = x, Latitude = y, C13_var = C13_var_liz)
coordinates(C13_iso_liz_krig_var_model) <- ~Longitude + Latitude
proj4string(C13_iso_liz_krig_var_model) <- CRS("+proj=longlat +datum=WGS84 +no_defs")
plot(C13_iso_liz_krig_var_model)
C13_iso_liz_krig_var_model_auto <- autofitVariogram(C13_var ~ 1, C13_iso_liz_krig_var_model)
summary(C13_iso_liz_krig_var_model_auto)
plot(C13_iso_liz_krig_var_model_auto)
crs_projected <- "+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs"
C13_iso_liz_krig_var_projected_model <- spTransform(C13_iso_liz_krig_var_model, CRS(crs_projected))
C13_iso_liz_krig_var_fin_model <- autoKrige(C13_var ~ 1, C13_iso_liz_krig_var_projected_model, C13_pred_pts_liz_projected_krige)
C13_iso_liz_krig_var_fin_model_sf <- st_as_sf(C13_iso_liz_krig_var_fin_model$krige_output)
C13_plot_Krige_var_model_liz <- ggplot() +
  geom_sf(data =  st_as_sf(land_liz), fill = "darkgrey", lwd = 0.05) +
  geom_sf(data = C13_iso_liz_krig_var_fin_model_sf, aes(color = var1.pred)) +
  geom_sf(data =  st_as_sf(land_liz), fill = "darkgrey", lwd = 0.05) +
  scale_color_viridis_c(option = "rocket", direction = -1) +
  theme_classic() +
  labs(x = "Longitude", y = "Latitude", color = expression(""* delta^13 * "C" * " var"))

# Generate Final Variance Raster
C13_iso_liz_krig_var_fin_model_sf_transformed <- st_transform(C13_iso_liz_krig_var_fin_model_sf, crs = crs('+proj=longlat +datum=WGS84 +no_defs'))
coords <- st_coordinates(C13_iso_liz_krig_var_fin_model_sf_transformed)
values <- C13_iso_liz_krig_var_fin_model_sf_transformed$var1.pred
template_raster <- raster(Depth)
C13_var_iso_liz_f <- rasterize(x = coords, y = template_raster, field = values, fun = mean)
plot(C13_var_iso_liz_f)
save(C13_var_iso_liz_f, file = "data/isoscape_data/C13_var_iso_liz_f.RData")
