# Assignment model functions
# Sourced by the assignment and sensitivity scripts. A seed is set before the TDF draws, and the distance to
# the highest-probability cell is a great-circle distance. peak_d13C_percentile records where the
# highest-probability cell sits in the d13C distribution of each fish's search area.

library(raster)
library(sf)
library(dplyr)
library(lubridate)
library(mvnmle)
library(parallel)
library(here)

# Settings ####
# Cores for the parallel runs (mclapply forks, which Windows does not support, so one core there)
CORES <- if (.Platform$OS.type == "windows") 1 else max(1, detectCores() - 2, na.rm = TRUE)

# Needed objects ####
# Load isoscapes, search-area layers and the template raster
load_assignment_inputs <- function() {
  e <- new.env()
  load(here("data/assignment_data/depth_raster.RData"), envir = e)
  load(here("data/assignment_data/monthly_sp.RData"), envir = e)
  load(here("data/assignment_data/sat_east_sf.RData"), envir = e)
  load(here("data/assignment_data/sat_west_sf.RData"), envir = e)
  load(here("data/assignment_data/sat_comb_sf.RData"), envir = e)
  load(here("data/move_data/tarpon_latitude/lat_CI_predictions.RData"), envir = e)
  e$iso_liz <- read.csv(here("data/isoscape_data/lizardfish_iso.csv"))

  # Template raster
  full <- extent(-99, -70, 22, 42)
  q <- 0.1
  e$x_seq <- seq(from = full@xmin, to = full@xmax, by = q)
  e$y_seq <- seq(from = full@ymin, to = full@ymax, by = q)
  e$tplate <- raster(nrows = length(e$y_seq), ncols = length(e$x_seq), ext = full, crs = crs(e$depth_raster))

  # Isoscapes (mean and variance)
  load(here("data/isoscape_data/C13_mean_iso_liz_f.RData"), envir = e)
  load(here("data/isoscape_data/C13_var_iso_liz_f.RData"), envir = e)
  load(here("data/isoscape_data/N15_mean_iso_liz_f.RData"), envir = e)
  load(here("data/isoscape_data/N15_var_iso_liz_f.RData"), envir = e)
  lmodC <- resample(e$C13_mean_iso_liz_f, e$tplate)
  lmodCvar <- resample(e$C13_var_iso_liz_f, e$tplate)
  lmodN <- resample(e$N15_mean_iso_liz_f, e$tplate)
  lmodNvar <- resample(e$N15_var_iso_liz_f, e$tplate)
  e$lmodC <- lmodC
  e$lmodN <- lmodN
  e$mu <- cbind(getValues(lmodC), getValues(lmodN))

  # Mask to 0-20 m depth
  depth_mask <- e$depth_raster <= 0 & e$depth_raster >= -20
  depth_mask[depth_mask == 0] <- NA
  e$lmodCvar_masked <- mask(lmodCvar, depth_mask)
  e$lmodNvar_masked <- mask(lmodNvar, depth_mask)

  # Correlation between d13C and d15N across the masked isoscape
  mu_masked <- cbind(getValues(mask(lmodC, depth_mask)), getValues(mask(lmodN, depth_mask)))
  mvn_masked <- mlest(mu_masked)
  e$r <- mvn_masked$sigmahat[1, 2] / (sqrt(mvn_masked$sigmahat[1, 1]) * sqrt(mvn_masked$sigmahat[2, 2]))
  e
}

# Latitude range ####
# Latitude range from the GAM predictions over the incorporation window
calculate_lat_range <- function(given_ordinal_date, given_latitude, month_count, lat_CI_predictions) {
  start_date <- given_ordinal_date - (month_count * 30)
  if (start_date < 1) start_date <- 365 + start_date
  start_date <- max(start_date, 1)
  if (start_date > given_ordinal_date) {
    relevant_predictions <- lat_CI_predictions %>% filter(Date_ord >= start_date | Date_ord <= given_ordinal_date)
  } else {
    relevant_predictions <- lat_CI_predictions %>% filter(Date_ord >= start_date & Date_ord <= given_ordinal_date)
  }
  list(
    range_75 = c(min(min(relevant_predictions$Lower_CI_75, na.rm = TRUE), given_latitude),
                 max(max(relevant_predictions$Upper_CI_75, na.rm = TRUE), given_latitude)),
    range_95 = c(min(min(relevant_predictions$Lower_CI_95, na.rm = TRUE), given_latitude),
                 max(max(relevant_predictions$Upper_CI_95, na.rm = TRUE), given_latitude))
  )
}

# Probability density ####
# Bivariate normal probability density for each cell
# x = value to assign, m = mean vector, v = var vector, r = corr coefficient, ras = raster extent
pdraster <- function(x, m, v, r, ras) {
  pd <- 1 / (2 * pi * sqrt(v[, 1]) * sqrt(v[, 2]) * sqrt(1 - r^2)) * exp(-(1 / (2 * (1 - r^2))) *
    ((x[1] - m[, 1])^2 / v[, 1] + (x[2] - m[, 2])^2 / v[, 2] - (2 * r * (x[1] - m[, 1]) * (x[2] - m[, 2])) /
      (sqrt(v[, 1]) * sqrt(v[, 2]))))
  setValues(ras, pd)
}

# Distance ####
# Great-circle distance (km)
gc_dist_km <- function(lon1, lat1, lon2, lat2) {
  rad <- pi / 180
  a <- sin((lat2 - lat1) * rad / 2)^2 + cos(lat1 * rad) * cos(lat2 * rad) * sin((lon2 - lon1) * rad / 2)^2
  2 * 6371.0088 * asin(pmin(1, sqrt(a)))
}

# Variance terms ####
# Variance terms for one TDF combination. The TDF draws use the seed, so assign_samples() and this
# function give the same values for the same seed. Tarpon SD terms come from all samples passed in.
variance_terms <- function(e, sample, Troph.fract.C, Troph.fract.N, seed = 1) {
  set.seed(seed)

  # Simulate distribution of trophic separations and fractionation
  Troph.diff <- 1
  Troph.diff.dist <- rnorm(10000, mean = Troph.diff, sd = 0.5)
  Troph.fract.dist.C <- rnorm(10000, mean = Troph.fract.C, sd = 0.5)
  Troph.fract.dist.N <- rnorm(10000, mean = Troph.fract.N, sd = 0.5)

  # Calibration errors
  CalibErrorC <- var(Troph.diff.dist * Troph.fract.dist.C)
  CalibErrorN <- var(Troph.diff.dist * Troph.fract.dist.N)

  # Inshore lizardfish variance (individual SD and measurement SD of 0.2)
  C.indL <- round(sd(e$iso_liz$d13C), digits = 2)
  N.indL <- round(sd(e$iso_liz$d15N), digits = 2)
  rep.liz.C.var <- C.indL^2 + 0.2^2
  rep.liz.N.var <- N.indL^2 + 0.2^2

  # Tarpon variance (individual SD and measurement SD of 0.2)
  C.ind <- round(sd(sample$d13C), digits = 2)
  N.ind <- round(sd(sample$d15N), digits = 2)
  var2.C <- 0.2^2 + C.ind^2
  var2.N <- 0.2^2 + N.ind^2

  # Pooled variances: isoscape and lizardfish (sigma_iso), plus tarpon (sigma_assign) and calibration (sigma_off_calib)
  iso_C <- getValues(e$lmodCvar_masked) + rep.liz.C.var
  iso_N <- getValues(e$lmodNvar_masked) + rep.liz.N.var
  list(Troph.diff = Troph.diff, C.indL = C.indL, N.indL = N.indL, C.ind = C.ind, N.ind = N.ind,
       CalibErrorC = CalibErrorC, CalibErrorN = CalibErrorN, iso_C = iso_C, iso_N = iso_N,
       cvar = iso_C + var2.C + CalibErrorC, nvar = iso_N + var2.N + CalibErrorN)
}

# Ranges of the variance terms over the 0-20 m isoscape cells, as reported in the Supplementary Materials
summarise_variance <- function(v) {
  data.frame(tarpon_sd_d13C = v$C.ind, tarpon_sd_d15N = v$N.ind,
             lizardfish_sd_d13C = v$C.indL, lizardfish_sd_d15N = v$N.indL,
             calibration_var_d13C = v$CalibErrorC, calibration_var_d15N = v$CalibErrorN,
             sigma_iso_d13C_min = min(v$iso_C, na.rm = TRUE), sigma_iso_d13C_max = max(v$iso_C, na.rm = TRUE),
             sigma_iso_d15N_min = min(v$iso_N, na.rm = TRUE), sigma_iso_d15N_max = max(v$iso_N, na.rm = TRUE),
             total_d13C_min = min(v$cvar, na.rm = TRUE), total_d13C_max = max(v$cvar, na.rm = TRUE),
             total_d15N_min = min(v$nvar, na.rm = TRUE), total_d15N_max = max(v$nvar, na.rm = TRUE),
             n_cells = sum(!is.na(v$cvar)))
}

# Assignment ####
# Assign one set of samples for one TDF combination.
# sample: data frame of all samples in the run (tarpon SD terms use all of them)
# idx: rows to assign; lat_ci: "95" (fin) or "75" (plasma); month_count: 5 (fin) or 3 (plasma)
assign_samples <- function(e, sample, Troph.fract.C, Troph.fract.N, idx = seq_len(nrow(sample)),
                           lat_ci = "95", month_count = 5, raster_dir = NULL, seed = 1) {
  v <- variance_terms(e, sample, Troph.fract.C, Troph.fract.N, seed = seed)

  # Raw values adjusted for trophic differences and fractionation
  testN <- sample$d15N + v$Troph.diff * Troph.fract.N
  testC <- sample$d13C + v$Troph.diff * Troph.fract.C
  x <- cbind(testC, testN)
  vars <- cbind(v$cvar, v$nvar)

  depth_mask <- e$depth_raster <= 0 & e$depth_raster >= -20
  lat_values <- raster::yFromRow(e$tplate, 1:nrow(e$tplate))

  out <- vector("list", length(idx))
  for (k in seq_along(idx)) {
    i <- idx[k]
    s <- sample[i, ]

    # Monthly utilization distribution for the sampling month (month.name is English in any locale)
    sample_month <- month.name[month(ymd(s$Date))]
    tplate_masked <- mask(setValues(e$tplate, rep(1, ncell(e$tplate))), e$monthly_sp[[sample_month]])

    # Contingent search area
    pt <- st_set_crs(st_as_sf(data.frame(Longitude = s$Longitude, Latitude = s$Latitude),
                              coords = c("Longitude", "Latitude")), 4326)
    intersects_west <- as.numeric(st_distance(pt, e$sat_west_sf)) == 0
    intersects_east <- as.numeric(st_distance(pt, e$sat_east_sf)) == 0
    if (intersects_west & intersects_east) {
      mask_polygon <- e$sat_comb_sf; assigned_polygon <- "Both"
    } else if (intersects_west) {
      mask_polygon <- e$sat_west_sf; assigned_polygon <- "EGOM"
    } else if (intersects_east) {
      mask_polygon <- e$sat_east_sf; assigned_polygon <- "ESeaboard"
    } else {
      next
    }
    tplate_masked <- mask(tplate_masked, as(mask_polygon, "Spatial"))

    # Latitude range predicted for the incorporation window
    lat_range <- calculate_lat_range(yday(s$Date), s$Latitude, month_count, e$lat_CI_predictions)[[paste0("range_", lat_ci)]]
    rows_out <- which(lat_values < lat_range[1] | lat_values > lat_range[2])
    if (length(rows_out) > 0) tplate_masked[cellFromRow(tplate_masked, rows_out)] <- NA

    # 0-20 m depth
    tplate_masked <- mask(depth_mask, tplate_masked)
    tplate_masked[tplate_masked == 0] <- NA
    valid_cells <- !is.na(getValues(tplate_masked))

    mu_mask <- e$mu
    vars_mask <- vars
    mu_mask[!valid_cells, ] <- NA
    vars_mask[!valid_cells, ] <- NA

    # Probability, odds and odds ratio relative to the most likely cell
    ProbCN <- pdraster(x[i, ], mu_mask, vars_mask, e$r, tplate_masked)
    OddTCN <- ProbCN / (1 - ProbCN)
    OddsR.high <- OddTCN / max(getValues(OddTCN), na.rm = TRUE)

    if (!is.null(raster_dir)) {
      writeRaster(OddsR.high, filename = sprintf("%s/probability_raster_%s_C_%.2f_N_%.2f.tif",
                                                 raster_dir, s$ID, Troph.fract.C, Troph.fract.N),
                  format = "GTiff", overwrite = TRUE)
    }

    # Highest-probability cell and distance from the capture location
    peak_cell <- which.max(getValues(OddsR.high))
    Mspot <- xyFromCell(ProbCN, peak_cell)
    # Share of search-area cells with isoscape d13C at or below that of the peak cell
    peak_d13C_percentile <- mean(e$mu[valid_cells, 1] <= e$mu[peak_cell, 1], na.rm = TRUE)
    out[[k]] <- data.frame(ID = s$ID, Latitude_new = Mspot[2], Longitude_new = Mspot[1],
                           dist = gc_dist_km(s$Longitude, s$Latitude, Mspot[1], Mspot[2]),
                           assigned_polygon = assigned_polygon,
                           Troph.fract.C = Troph.fract.C, Troph.fract.N = Troph.fract.N,
                           peak_d13C_percentile = peak_d13C_percentile)
  }
  bind_rows(out)
}

# Parallel runs ####
# Split samples across cores (same seed in each, so the TDF draws match a single run).
# A worker that fails or is killed returns an error object instead of rows, so check every sample came back.
assign_samples_parallel <- function(e, sample, Troph.fract.C, Troph.fract.N, cores = CORES, ...) {
  chunks <- split(seq_len(nrow(sample)), cut(seq_len(nrow(sample)), cores, labels = FALSE))
  res <- mclapply(chunks, function(ix) assign_samples(e, sample, Troph.fract.C, Troph.fract.N, idx = ix, ...), mc.cores = cores)
  stopifnot(all(vapply(res, is.data.frame, logical(1))))
  out <- bind_rows(res)
  stopifnot(nrow(out) == nrow(sample), setequal(out$ID, sample$ID))
  out
}

# TDF grid ####
# Assign the samples for 121 TDF combinations (+/- 0.5 per mil around C_base and N_base, in 0.1 steps)
assign_grid <- function(e, sample, C_base, N_base, lat_ci, month_count) {
  grid <- expand.grid(Troph.fract.N = round(seq(N_base - 0.5, N_base + 0.5, by = 0.1), 2),
                      Troph.fract.C = round(seq(C_base - 0.5, C_base + 0.5, by = 0.1), 2))
  res <- mclapply(seq_len(nrow(grid)), function(g)
    assign_samples(e, sample, grid$Troph.fract.C[g], grid$Troph.fract.N[g], lat_ci = lat_ci, month_count = month_count, seed = 1),
    mc.cores = CORES)
  stopifnot(all(vapply(res, is.data.frame, logical(1))))
  out <- bind_rows(res)
  stopifnot(nrow(out) == nrow(sample) * nrow(grid))
  out
}

# Regions ####
# Regions used in the manuscript
categorize_location <- function(lat, lon) {
  ifelse(lat > 29 & lon < -83.3, "N. Gulf of Mexico",
    ifelse(lat > 26.11 & lat < 30.73 & lon < -81.79, "CW FL",
      ifelse(lat > 26 & lat < 30.73 & lon > -81.4, "E FL",
        ifelse(lat < 26.11, "S FL",
          ifelse(lat > 29.4 & lon > -81.5, "Mid Atlantic", "Unlabeled")))))
}
