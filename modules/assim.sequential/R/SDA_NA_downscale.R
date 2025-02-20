#' @description
#' This function helps to average the ERA5 data based on the start and end dates, and convert it to the GeoTIFF file.
#' @title Average.ERA5.2.GeoTIFF
#' 
#' @param start.date character: start point of when to average the data (e.g., 2012-01-01).
#' @param end.date character: end point of when to average the data (e.g., 2021-12-31).
#' @param in.path character: the directory where your ERA5 data stored (they should named as ERA5_YEAR.nc).
#' @param outdir character: the output directory where the averaged GeoTIFF file will be generated.
#'
#' @return character: path to the exported GeoTIFF file.
#' 
#' @examples
#' @export
#' @author Dongchen Zhang
Average.ERA5.2.GeoTIFF <- function (start.date, end.date, in.path, outdir) {
  # open ERA5 nc file as geotiff format for referencing crs and ext.
  ERA5.tiff <- terra::rast(file.path(in.path, paste0("ERA5_", lubridate::year(start.date), ".nc")))
  dates <- seq(start.date, end.date, "1 year")
  if (length(dates) < 2) {
    PEcAn.logger::logger.info("There is no time range to be calculated!")
    return(NA)
  }
  # initialize final outcomes.
  temp.all <- precip.all <- srd.all <- dewpoint.all <- c()
  for (i in 2:length(dates)) {
    # initialize start and end dates for the current period
    if (i == 1) {
      start <- start.date
    } else {
      start <- as.Date(paste0(lubridate::year(dates[i]), "-01-01"))
    }
    if (i == length(dates)) {
      end <- end.date
    } else {
      end <- as.Date(paste0(lubridate::year(dates[i]), "-12-31"))
    }
    # loop over years.
    for (j in seq_along(dates)) {
      # open ERA5 nc file.
      met.nc <- ncdf4::nc_open(file.path(in.path, paste0("ERA5_", lubridate::year(dates[i]), ".nc")))
      # find index for the date.
      times <- as.POSIXct(met.nc$dim$time$vals*3600, origin="1900-01-01 00:00:00", tz = "UTC")
      time.inds <- which(lubridate::date(times) >= start & lubridate::date(times) <= end)
      # extract temperature.
      PEcAn.logger::logger.info("entering temperature.")
      temp.all <- abind::abind(temp.all, apply(ncdf4::ncvar_get(met.nc, "t2m")[,,,time.inds], c(1,2,4), mean), along = 3)
      # extract precipitation.
      PEcAn.logger::logger.info("entering precipitation.")
      precip.all <- abind::abind(precip.all, apply(ncdf4::ncvar_get(met.nc, "tp")[,,,time.inds], c(1,2,4), mean), along = 3)
      # extract shortwave solar radiation.
      PEcAn.logger::logger.info("entering solar radiation.")
      srd.all <- abind::abind(srd.all, apply(ncdf4::ncvar_get(met.nc, "ssrd")[,,,time.inds], c(1,2,4), mean), along = 3)
      # extract dewpoint.
      PEcAn.logger::logger.info("entering dewpoint.")
      dewpoint.all <- abind::abind(dewpoint.all, apply(ncdf4::ncvar_get(met.nc, "d2m")[,,,time.inds], c(1,2,4), mean), along = 3)
      # close the NC connection.
      ncdf4::nc_close(met.nc)
    }
  }
  # aggregate across time.
  # temperature.
  temp <- apply(temp.all, c(1, 2), mean)
  temp <- PEcAn.utils::ud_convert(temp, "K", "degC")
  # precipitation.
  precip <- apply(precip.all, c(1, 2), mean)
  # solar radiation.
  srd <- apply(srd.all, c(1, 2), mean)
  # dewpoint.
  dewpoint <- apply(dewpoint.all, c(1, 2), mean)
  dewpoint <- PEcAn.utils::ud_convert(dewpoint, "K", "degC")
  # convert dew point to relative humidity.
  beta <- (112 - (0.1 * temp) + dewpoint) / (112 + (0.9 * temp))
  relative.humidity <- beta ^ 8
  VPD <- PEcAn.data.atmosphere::get.vpd(100*relative.humidity, temp)
  # combine together.
  PEcAn.logger::logger.info("Aggregate maps.")
  met.rast <- c(terra::rast(matrix(temp, nrow = dim(temp)[2], ncol = dim(temp)[1], byrow = T)),
                terra::rast(matrix(precip, nrow = dim(precip)[2], ncol = dim(precip)[1], byrow = T)),
                terra::rast(matrix(srd, nrow = dim(srd)[2], ncol = dim(srd)[1], byrow = T)),
                terra::rast(matrix(VPD, nrow = dim(VPD)[2], ncol = dim(VPD)[1], byrow = T)))
  # adjust crs and extents.
  terra::crs(met.rast) <- terra::crs(ERA5.tiff)
  terra::ext(met.rast) <- terra::ext(ERA5.tiff)
  names(met.rast) <- c("temp", "prec", "srad", "vapr")
  # write into geotiff file.
  terra::writeRaster(met.rast, file.path(outdir, paste0("ERA5_met_", lubridate::year(end.date), ".tiff")))
  # end.
  gc()
  return(file.path(outdir, paste0("ERA5_met_", lubridate::year(end.date), ".tiff")))
}



# assemble covariates from different spatial scales/resolutions and crs.
# Here is an example of the `cov.tif.file.list` object:
# cov.tif.file.list <- list(LC = list(dir = "/projectnb/dietzelab/dongchen/anchorSites/downscale/MODIS_NLCD_LC.tif",
#                                     var.name = "LC"),
#                           year_since_disturb = list(dir = "/projectnb/dietzelab/dongchen/anchorSites/downscale/MODIS_LC/outputs/age.tif",
#                                                     var.name = "year_since_disturb"),
#                           agb = list(dir = "/projectnb/dietzelab/dongchen/anchorSites/downscale/AGB/agb.tif",
#                                      var.name = "agb"),
#                           twi = list(dir = "/projectnb/dietzelab/dongchen/anchorSites/downscale/TWI/TWI_resample.tiff",
#                                      var.name = "twi"),
#                           met = list(dir = paths[i],
#                                      var.name = c("temp", "prec", "srad", "vapr")),
#                           soil = list(dir = "/projectnb/dietzelab/dongchen/anchorSites/downscale/SoilGrids.tif",
#                                       var.name = c("PH", "N", "SOC", "Sand")))
# This function helps to stack target data layers from various GeoTIFF maps to a single map 
# cropped and projected to the `base.map`. It also enables the normalization feature to facilitate the ML process.


#' @description
#' This function helps to stack target data layers from various GeoTIFF maps (with different extents, CRS, and resolutions) to a single map.
#' @title stack.covariates.2.geotiff
#' 
#' @param outdir character: the output directory where the stacked GeoTIFF file will be generated.
#' @param base.map.dir character: path to the GeoTIFF file within which the extents and CRS will be used to generate the final map.
#' @param cov.tif.file.list list: a list contains sub-lists with each including path to the corresponding map and the variables to be extracted (e.g., list(LC = list(dir = "path/to/landcover.tiff", var.name = "LC")).
#' @param normalize boolean: decide if we want to normalize each data layer, the default is TRUE.
#' @param cores numeric: how many CPus to be used in the calculation, the default is the total CPU number you have.
#'
#' @return path to the exported GeoTIFF file.
#' 
#' @examples
#' @author Dongchen Zhang
stack.covariates.2.geotiff <- function(outdir, base.map.dir, cov.tif.file.list, normalize = T, cores = parallel::detectCores()) {
  # create the folder if it doesn't exist.
  if (!file.exists(outdir)) {
    dir.create(outdir)
  }
  # parallel loop.
  # register parallel nodes.
  if (cores > length(cov.tif.file.list)) {
    cores <- length(cov.tif.file.list)
  }
  cl <- parallel::makeCluster(as.numeric(cores))
  doSNOW::registerDoSNOW(cl)
  #progress bar.
  pb <- utils::txtProgressBar(min=1, max=length(cov.tif.file.list), style=3)
  progress <- function(n) utils::setTxtProgressBar(pb, n)
  opts <- list(progress=progress)
  # foreach loop.
  paths <- foreach::foreach(f = cov.tif.file.list, 
                            .packages=c("Kendall", "terra"),
                            .options.snow=opts) %dopar% {
                              # load the base map.
                              base.map <- terra::rast(base.map.dir)
                              # read geotif file.
                              temp.rast <- terra::rast(f$dir)
                              # normalize.
                              if (normalize & !"LC" %in% f$var.name) {
                                nx <- terra::minmax(temp.rast)
                                temp.rast <- (temp.rast - nx[1,]) / (nx[2,] - nx[1,])
                              }
                              # set name to layers if we set it up in advance.
                              # otherwise the original layer name will be used.
                              if (!is.null(f$var.name)) {
                                names(temp.rast) <- f$var.name
                              }
                              # raster operations.
                              terra::crs(temp.rast) <- terra::crs(base.map)
                              temp.rast <- terra::crop(temp.rast, base.map)
                              temp.rast <- terra::resample(temp.rast, base.map)
                              # write the raster into disk.
                              file.name <- paste0(f$var.name, collapse = "_")
                              path <- file.path(outdir, paste0(file.name, ".tiff"))
                              terra::writeRaster(temp.rast, path)
                              return(path)
                            } %>% unlist
  # stop parallel.
  parallel::stopCluster(cl)
  foreach::registerDoSEQ()
  gc()
  # combine rasters.
  all.rast <- terra::rast(paths)
  # write all covariates into disk.
  terra::writeRaster(all.rast, file.path(outdir, "covariates.tiff"), overwrite = T)
  # remove previous tiff files.
  unlink(paths)
  # return results.
  return(file.path(outdir, "covariates.tiff"))
}

#' @description
#' convert settings to geospatial points in terra.
#' @title pecan.settings.2.pts
#' 
#' @param settings PEcAn settings: either a character that points to the settings or the actual settings object will be accepted.
#'
#' @return terra spatial points object.
#' 
#' @examples
#' @author Dongchen Zhang
pecan.settings.2.pts <- function(settings) {
  if (is.character(settings)) {
    # read settings.
    settings <- PEcAn.settings::read.settings(settings)
  }
  # grab lat/lon.
  site.locs <- settings$run %>% purrr::map('site') %>% 
    purrr::map_dfr(~c(.x[['lon']],.x[['lat']]) %>% as.numeric)%>% 
    t %>% `colnames<-`(c("Lon","Lat")) %>% as.data.frame()
  # convert lat/lon to terra::vect.
  pts <- terra::vect(site.locs, geom = c("Lon", "Lat"), crs = "EPSG:4326")
  return(pts)
}

#' @description
#' This function helps to build the data frame (pixels by data columns) for only vegetated pixels to improve the efficiency.
#' Note that the `LC` field using the `MODIS land cover` observations (MCD12Q1.061) must be supplied in the covariates to make this function work.
#' @title stack.covariates.2.df
#' 
#' @param rast.dir character: a character that points to the covariates raster file generated by the `stack.covariates.2.geotiff` function.
#' @param cores numeric: how many CPus to be used in the calculation, the default is the total CPU number you have.
#'
#' @return list containing the data frame of covariates for vegetated pixels and the corresponding index of the pixels.
#' 
#' @examples
#' @author Dongchen Zhang
stack.covariates.2.df <- function(rast.dir, cores = parallel::detectCores()) {
  # load maps.
  all.rast <- terra::rast(rast.dir)
  # parallel loop.
  layer.names <- names(all.rast)
  # register parallel nodes.
  if (cores > length(layer.names)) {
    cores <- length(layer.names)
  }
  cl <- parallel::makeCluster(as.numeric(cores))
  doSNOW::registerDoSNOW(cl)
  #progress bar.
  pb <- utils::txtProgressBar(min=1, max=length(layer.names), style=3)
  progress <- function(n) utils::setTxtProgressBar(pb, n)
  opts <- list(progress=progress)
  # foreach loop.
  vecs <- foreach::foreach(r = seq_along(layer.names), 
                           .packages=c("Kendall", "terra"),
                           .options.snow=opts) %dopar% {
                             all.rast <- terra::rast(rast.dir)
                             temp.vec <- matrix(all.rast[[r]], byrow = T)
                             na.inds <- which(is.na(temp.vec))
                             # if it's LC layer.
                             if ("LC" == names(all.rast)[r]) {
                               non.veg.inds <- which(! temp.vec %in% 1:8)
                               na.inds <- unique(c(na.inds, non.veg.inds))
                             }
                             return(list(vec = temp.vec,
                                         na.inds = na.inds))
                           }
  # stop parallel.
  parallel::stopCluster(cl)
  foreach::registerDoSEQ()
  gc()
  # grab uniqued NA index.
  na.inds <- vecs %>% purrr::map("na.inds") %>% unlist %>% unique
  # remove NA from each covariate.
  cov.vecs <- vecs %>% purrr::map(function(v){
    return(v$vec[-na.inds])
  }) %>% dplyr::bind_cols() %>% `colnames<-`(layer.names) %>% as.data.frame()
  non.na.inds <- seq_along(matrix(all.rast[[1]]))[-na.inds]
  return(list(df = cov.vecs, non.na.inds = non.na.inds))
}

#' @description
#' This function helps to create the training dataset of specific variable type and locations for downscaling.
#' TODO: There will be a ratio argument between training and testing samples to testify the ML regression accuracy.
#' @title prepare.train.dat
#' 
#' @param settings character: physical path that points to the pecan settings XML file.
#' @param analysis numeric: data frame (rows: ensemble member; columns: site*state_variables) of updated ensemble analysis results from the `sda_enkf` function.
#' @param covariates.dir character: path to the exported covariates GeoTIFF file.
#' @param variable character: name of state variable. It should match up with the column names of the analysis data frame. 
#'
#' @return matrix within which the first sets of columns contain values of state variables for each ensemble mebers of every site, and the rest columns contain the corresponding covariates.
#' 
#' @examples
#' @author Dongchen Zhang
prepare.train.dat <- function(settings, analysis, covariates.dir, variable) {
  # convert settings into geospatial points.
  pts <- pecan.settings.2.pts(settings)
  # read covariates.
  cov.rast <- terra::rast(covariates.dir)
  # extract covariates by locations.
  predictors <- as.data.frame(terra::extract(cov.rast, pts, ID = FALSE))
  covariate_names <- names(predictors)
  if ("ID" %in% covariate_names) {
    rm.ind <- which("ID" %in% covariate_names)
    covariate_names <- covariate_names[-rm.ind]
    predictors <- predictors[,-rm.ind]
  }
  # grab carbon data.
  var.dat <- analysis[,which(colnames(analysis) == variable)] %>% t %>% 
    as.data.frame() %>% `colnames<-`(paste0("ensemble", seq(nrow(analysis))))
  # combine carbon and predictor.
  full_data <- cbind(var.dat, predictors)
  return(full_data)
}

#' @description
#' This function helps to train the ML model across ensemble members in parallel.
#' @title parallel.rf.train
#' 
#' @param full_data numeric: the matrix generated using the `prepare.train.dat` function.
#' @param cores numeric: how many CPus to be used in the calculation, the default is the total CPU number you have.
#'
#' @return list of trained models across ensemble members.
#' 
#' @examples
#' @author Dongchen Zhang
parallel.rf.train <- function(full_data, cores = parallel::detectCores()) {
  # grab ensemble and predictor index.
  col.names <- colnames(full_data)
  ensemble.inds <- which(grepl("ensemble", col.names, fixed = TRUE))
  predictor.inds <- seq_along(col.names)[-ensemble.inds]
  # parallel train.
  # register parallel nodes.
  if (cores > length(ensemble.inds)) {
    cores <- length(ensemble.inds)
  }
  cl <- parallel::makeCluster(as.numeric(cores))
  doSNOW::registerDoSNOW(cl)
  #progress bar.
  pb <- utils::txtProgressBar(min=1, max=length(ensemble.inds), style=3)
  progress <- function(n) utils::setTxtProgressBar(pb, n)
  opts <- list(progress=progress)
  # foreach loop.
  models <- foreach::foreach(i = ensemble.inds, 
                             .packages=c("Kendall", "stats", "randomForest"),
                             .options.snow=opts) %dopar% {
                               ensemble_col <- col.names[ensemble.inds[i]]
                               formula <- stats::as.formula(paste(ensemble_col, "~", paste(col.names[predictor.inds], collapse = " + ")))
                               randomForest::randomForest(formula,
                                                          data = full_data,
                                                          ntree = 1000,
                                                          na.action = stats::na.omit,
                                                          keep.forest = TRUE,
                                                          importance = TRUE)
                             }
  # stop parallel.
  parallel::stopCluster(cl)
  foreach::registerDoSEQ()
  gc()
  return(models)
}

#' @description
#' This function helps to predict the target variable observations based on the covariates.
#' The prediction is working in parallel across vegetated pixels.
#' @title parallel.prediction
#' 
#' @param base.map.dir character: path to the GeoTIFF file within which the extents and CRS will be used to generate the ensemble maps.
#' @param models list: trained models across ensemble members generated by the `parallel.rf.train` function.
#' @param cov.vecs: numeric: data frame containing covaraites across vegetated pixels generated from the `stack.covariates.2.df` function.
#' @param non.na.inds numeric: the corresponding index of vegetated pixels generated from the `stack.covariates.2.df` function.
#' @param outdir character: the output directory where the downscaled maps will be stored.
#' @param name list: containing the time and variable name to create the final GeoTIFF file name.
#' @param cores numeric: how many CPus to be used in the calculation, the default is the total CPU number you have.
#'
#' @return paths to the ensemble downscaled maps.
#' 
#' @examples
#' @author Dongchen Zhang
parallel.prediction <- function(base.map.dir, models, cov.vecs, non.na.inds, outdir, name, cores = parallel::detectCores()) {
  # load base map.
  base.map <- terra::rast(base.map.dir)
  dims <- dim(base.map)
  # setup progress bar for ensemble members.
  pb <- utils::txtProgressBar(min = 0, max = length(models), style = 3)
  paths <- c()
  # loop over ensemble members.
  for (i in seq_along(models)) {
    # update progress bar.
    utils::setTxtProgressBar(pb, i)
    # go to the next if the current file has already been generated.
    file.name <- paste0(c("ensemble", i, name$time, name$variable), collapse = "_")
    if (file.exists(file.path(outdir, paste0(file.name, ".tiff")))) {
      next
    }
    # register parallel nodes.
    cl <- parallel::makeCluster(cores)
    doSNOW::registerDoSNOW(cl)
    # foreach parallel.
    model <- models[[i]]
    output <- foreach::foreach(d=itertools::isplitRows(cov.vecs, chunks=cores),
                               .packages=c("stats", "randomForest")) %dopar% {
                                 stats::predict(model, d)
                               } %>% unlist
    # export to geotiff map.
    vec <- rep(NA, dims[1]*dims[2])
    vec[non.na.inds] <- output
    map <- terra::rast(matrix(vec, dims[1], dims[2], byrow = T))
    terra::ext(map) <- terra::ext(base.map)
    terra::crs(map) <- terra::crs(base.map)
    terra::writeRaster(map, file.path(outdir, paste0(file.name, ".tiff")))
    paths <- c(paths, file.path(outdir, paste0(file.name, ".tiff")))
    # stop parallel.
    parallel::stopCluster(cl)
    foreach::registerDoSEQ()
    gc()
  }
  return(paths)
}

#' @description
#' This is the main function to execute the RF training and prediction. 
#' Note it will be deployed by each node you requested if the qsub feature is enabled below.
#' @title downscale.rf.main
#' 
#' @param settings character: physical path that points to the pecan settings XML file.
#' @param analysis numeric: data frame (rows: ensemble member; columns: site*state_variables) of updated ensemble analysis results from the `sda_enkf` function.
#' @param covariates.dir character: path to the exported covariates GeoTIFF file.
#' @param time: character: the time tag used to differentiate the outputs from others.
#' @param variable: character: name of state variable. It should match up with the column names of the analysis data frame. 
#' @param outdir character: the output directory where the downscaled maps will be stored.
#' @param base.map.dir character: path to the GeoTIFF file within which the extents and CRS will be used to generate the ensemble maps.
#' @param cores numeric: how many CPus to be used in the calculation, the default is the total CPU number you have.
#'
#' @return paths to the ensemble downscaled maps.
#' 
#' @examples
#' @author Dongchen Zhang
downscale.rf.main <- function(settings, analysis, covariates.dir, time, variable, outdir, base.map.dir, cores = parallel::detectCores()) {
  # create folder specific for the time and carbon type.
  folder.name <- file.path(outdir, paste0(c(variable, time), collapse = "_"))
  if (!file.exists(folder.name)) {
    dir.create(folder.name)
  }
  # prepare training data.
  PEcAn.logger::logger.info("Preparing training data.")
  full_data <- prepare.train.dat(settings = settings, 
                                 analysis = analysis, 
                                 covariates.dir = covariates.dir, 
                                 variable = variable)
  # convert LC into factor.
  if ("LC" %in% colnames(full_data)) {
    full_data[,"LC"] <- factor(full_data[,"LC"])
  }
  # parallel train.
  PEcAn.logger::logger.info("Parallel training.")
  models <- parallel.rf.train(full_data = full_data, cores = cores)
  # save trained models for future analysis.
  # saveRDS(models, file.path(folder.name, "rf_models.rds"))
  save(models, file = file.path(folder.name, "rf_models.Rdata"))
  # convert stacked covariates geotiff file into data frame.
  PEcAn.logger::logger.info("Converting geotiff to df.")
  cov.df <- stack.covariates.2.df(rast.dir = covariates.dir, cores = cores)
  # reconstruct LC because of the computation accuracy.
  cov.df$df$LC[which(cov.df$df$LC < 1)] <- 0
  # convert LC into factor.
  if ("LC" %in% colnames(cov.df$df)) {
    cov.df$df[,"LC"] <- factor(cov.df$df[,"LC"])
  }
  # parallel prediction.
  PEcAn.logger::logger.info("Parallel prediction.")
  paths <- parallel.prediction(base.map.dir = base.map.dir, 
                               models = models, 
                               cov.vecs = cov.df$df, 
                               non.na.inds = cov.df$non.na.inds, 
                               outdir = folder.name, 
                               name = list(time = time, variable = variable), 
                               cores = cores)
  # calculate mean and std.
  PEcAn.logger::logger.info("Calculate mean and std.")
  ras.all <- terra::rast(paths)
  mean <- terra::app(ras.all, "mean")
  std <- terra::app(ras.all, "std")
  # write into geotiff files.
  image.base.name <- paste0(time, "_", variable, ".tiff")
  terra::writeRaster(mean, filename = file.path(folder.name, paste0("mean_", image.base.name)))
  terra::writeRaster(std, filename = file.path(folder.name, paste0("std_", image.base.name)))
  return(list(ensemble.prediction.files = paths,
              mean.prediction.file = file.path(folder.name, paste0("mean_", image.base.name)),
              std.prediction.file = file.path(folder.name, paste0("std_", image.base.name))))
}

#' @description
#' This qsub function helps to run the submitted qsub jobs for running the downscale.rf.main function.
#' @title downscale.qsub.main
#' 
#' @param folder.path Character: physical path to which the job file is located.
#' 
#' @examples
#' @export
#' @author Dongchen Zhang
downscale.qsub.main <- function(folder.path) {
  dat <- readRDS(file.path(folder.path, "dat.rds"))
  out <- downscale.rf.main(dat$settings, dat$analysis.yr, dat$covariates.dir, lubridate::year(dat$time), dat$variable, dat$outdir, dat$base.map.dir, dat$cores)
  saveRDS(out, file.path(folder.path, "res.rds"))
}