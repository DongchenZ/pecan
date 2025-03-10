
#' Temporal downscaling of daily or subdaily CF met data
#'
#' @param cfmet data frame with CF variables generated by \code{\link{load.cfmet}}
#' @param output.dt time step (hours) for output
#' @param lat latitude (for calculating solar radiation)
#' @param ... ignored
#'
#' @return downscaled result
#' @export
#' @author David LeBauer
cfmet.downscale.time <- function(cfmet, output.dt = 1, lat = lat, ...) {
  ## time step
  dt_hr <- as.numeric(round(difftime(cfmet$date[2], cfmet$date[1],  units = "hours")))

  if (dt_hr == output.dt) {
    downscaled.result <- cfmet
  }

#   if("specific_humidity" %in% colnames(cfmet) & (!"relative_humidity" %in% colnames(cfmet))){
#     cfmet$relative_humidity <- cfmet[,list(qair2rh(qair = specific_humidity,
#                                                    temp = PEcAn.utils::ud_convert(air_temperature, "Kelvin", "Celsius"),
#                                                    press = PEcAn.utils::ud_convert(air_pressure, "Pa", "millibar"))]
#   }

  if(dt_hr > output.dt & dt_hr <= 6) {
    downscaled.result <- cfmet.downscale.subdaily(subdailymet = cfmet, output.dt = output.dt)
  } else if (dt_hr > 6 & dt_hr < 24) {
    # cfmet <- cfmet[,list(air_temperature_max = max(air_temperature), air_temperature_min =
    # min(air_temperature), ), by = 'year,doy']) dt_hr <- 24
   PEcAn.logger::logger.error("timestep of input met data is between 6 and 24 hours.\n", "PEcAn will automatically convert this to daily data\n",
                 "you should confirm validity of downscaling, in particular that min / max temperatures are realistic")
  }

  if (dt_hr == 24) {
    if (all(c("tmax", "tmin") %in% colnames(cfmet))) {
      nm <- colnames(cfmet)
      nm[nm == "tmax"] <- "air_temperature_max"
      nm[nm == "tmin"] <- "air_temperature_min"
      colnames(cfmet) <- nm
    }
    downscaled.result <- cfmet.downscale.daily(dailymet = cfmet, output.dt = output.dt, lat = lat)
  } else if (dt_hr > 24) {
   PEcAn.logger::logger.error("only daily and sub-daily downscaling supported")
  }

  return(downscaled.result)
} # cfmet.downscale.time


##' Subdaily to hourly (or less) downscaling
##'
##' Uses simple spline to interpolate variables with diurnal variability, otherwise uses averaging or repeating
##' for variables with no clear diurnal pattern. For all variables except temperature, negative values are set to zero.
##'
##' @param subdailymet data frame with climate variables queried from \code{\link{load.cfmet}}
##' @param output.dt output timestep. default is one hour
##' @export
##' @return weather file with subdaily met variables rescaled to output time step
##' @author David LeBauer
cfmet.downscale.subdaily <- function(subdailymet, output.dt = 1) {
  ## converting surface_downwelling_shortwave_flux_in_air from W/m2 avg to PPFD
  new.date <- subdailymet %>%
    dplyr::group_by(.data$year, .data$month, .data$day, .data$doy) %>%
    dplyr::group_modify(~data.frame(hour = 0:(23 / output.dt) / output.dt))

  new.date$date <- lubridate::ymd_h(paste(new.date$year, new.date$month, new.date$day, new.date$hour))

  downscaled.result <- list()
  tint <- nrow(new.date)/ nrow(subdailymet)
  if(all(c("eastward_wind", "northward_wind") %in% colnames(subdailymet))){
    if(!"wind_speed" %in% colnames(subdailymet)){
      subdailymet$wind_speed <- sqrt(subdailymet$northward_wind^2 + subdailymet$eastward_wind^2)
    }
    downscaled.result[["northward_wind"]] <- rep(subdailymet$northward_wind, each = tint)
    downscaled.result[["eastward_wind"]]  <- rep(subdailymet$eastward_wind, each = tint)
  } else if (!'wind_speed' %in% colnames(subdailymet)){
   PEcAn.logger::logger.error("no wind speed data")
  }
  downscaled.result[["wind_speed"]] <- rep(subdailymet$wind_speed, each = tint)

  solarMJ <- PEcAn.utils::ud_convert(subdailymet$surface_downwelling_shortwave_flux_in_air, paste0("W ", tint, "h"), "MJ" )
  PAR <- 0.486 * solarMJ ## Cambell and Norman 1998 p 151, ch 10
  subdailymet$ppfd <- PEcAn.utils::ud_convert(PAR, "mol s", "micromol h")
  downscaled.result[["ppfd"]] <- subdailymet$ppfd

  downscaled.result[["surface_downwelling_shortwave_flux_in_air"]] <- subdailymet$surface_downwelling_shortwave_flux_in_air


  for(var in c("air_pressure", "specific_humidity",
               "precipitation_flux", "air_temperature",
               "surface_downwelling_shortwave_flux_in_air", "ppfd", "relative_humidity")){
    if(var %in% colnames(subdailymet)){
      ## convert units from subdaily to hourly
      hrscale <- ifelse(var %in% c("surface_downwelling_shortwave_flux_in_air", "precipitation_flux"),
                        output.dt, 1)

      f <- stats::splinefun(as.double(subdailymet$date), (subdailymet[[var]] / hrscale), method = "monoH.FC")
      downscaled.result[[var]] <- f(as.double(new.date$date))
      downscaled.result[[var]][downscaled.result[[var]] < 0] <- 0
      if (var == "relative_humidity") {
        downscaled.result[[var]][downscaled.result[[var]] > 100] <- 100
      }
    }
  }

  downscaled.result <- cbind(new.date, as.data.frame(downscaled.result))
} # cfmet.downscale.subdaily


#' Internal helper to downscale a single row from a daily file
#'
#' @param df one row from dailymet
#' @param tseq vector of hours at which to estimate
#' @param lat latitude
#'
#' @return df with one row for each hour in `tseq`
#'
downscale_one_cfmet_day <- function(df, tseq, lat) {
  if (nrow(df) != 1) {
    PEcAn.logger::logger.error("df must be a one-row dataframe")
  }
  if (length(unique(diff(tseq))) != 1) {
    PEcAn.logger::logger.error("tseq has uneven intervals")
  }

  n <- length(tseq)

  light <- lightME(DOY = df$doy, t.d = tseq, lat = lat) %>%
    as.data.frame() %>%
    dplyr::mutate(
      Itot = .data$I.dir + .data$I.diff,
      resC2 = (.data$Itot - min(.data$Itot)) / max(.data$Itot))

  rhscale <- (cos(2 * pi * (tseq - 10) / n) + 1) / 2

  data.frame(
    hour = tseq,
    solarR = rep(
        df$surface_downwelling_shortwave_flux_in_air * 2.07 * 10^5 / 36000,
        each = n)
      * light$resC2,
    # Note: When dt >= 6, all times get the same T prediction
    Temp = df$tmin +
      (sin(2 * pi * (tseq - 10) / n) + 1) /
      2 * (df$tmax - df$tmin),
    relative_humidity = df$rhmin + rhscale * (df$rhmax - df$rhmin),
    # TODO: Why do we divide by n?
    #   isn't precipitation_flux already an intensity?
    precipitation_flux = rep(df$precipitation_flux / n, each = n),
    wind = rep(df$wind_speed, each = n)) %>%
  dplyr::mutate(
    # TODO: Computation of solarR above already multiplies by resC2.
    #   Is multiplying it again here really correct?
    #   That's how the old data.table version did it
    #   (once when computing `solarR` and again when computing `SolarR`),
    #   so keeping it until proven wrong.
    downwelling_photosynthetic_photon_flux = .data$solarR * light$resC2)
}


##' Simple, Fast Daily to Hourly Climate Downscaling
##'
##' Based on weach family of functions but 5x faster than weachNEW,
##' and requiring metric units (temperature in Kelvins on input and celsius on
##'   output, windspeed in kph, precip in mm, relative humidity as fraction).
##' Derived from the weachDT function in the BioCro package.
##'
##' @param dailymet data frame with climate variables
##' @param output.dt output timestep
##' @param lat latitude (for calculating solar radiation)
#'
#' @importFrom dplyr %>%
#' @importFrom rlang .data
##' @export
##' @return weather file with subdaily timesteps
##' @author David LeBauer
cfmet.downscale.daily <- function(dailymet, output.dt = 1, lat) {

  tint <- 24/output.dt
  tseq <- seq(from = 0, to = 23, by = output.dt)

  if (all(c("air_temperature_max", "air_temperature_min") %in% colnames(dailymet))) {
    nm <- colnames(dailymet)
    nm[nm == "air_temperature_max"] <- "tmax"
    nm[nm == "air_temperature_min"] <- "tmin"
    colnames(dailymet) <- nm
  }

  if (! "wind_speed" %in% colnames(dailymet)) {
    if (all(c("eastward_wind", "northward_wind") %in% colnames(dailymet))) {
      dailymet$wind_speed <- sqrt(dailymet$northward_wind^2 + dailymet$eastward_wind^2)
    } else {
     PEcAn.logger::logger.error("no wind_speed found in daily met dataset")
    }
  }

  return(dailymet %>%
    dplyr::mutate(
      qmin = rh2qair(rh = .data$relative_humidity / 100, T = .data$tmin),
      qmax = rh2qair(rh = .data$relative_humidity / 100, T = .data$tmax),
      pressure = PEcAn.utils::ud_convert(.data$air_pressure, "Pa", "millibar"),
      rhmin = qair2rh(.data$qmin, .data$air_temperature, .data$pressure),
      rhmax = qair2rh(.data$qmax, .data$air_temperature, .data$pressure)) %>%
    dplyr::group_by(.data$year, .data$doy) %>%
    dplyr::group_modify(~downscale_one_cfmet_day(.x, tseq, lat), .keep = TRUE) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      air_temperature = PEcAn.utils::ud_convert(.data$Temp, "kelvin", "celsius")) %>%
    dplyr::select(
      "year", "doy", "hour",
      "downwelling_photosynthetic_photon_flux",
      "air_temperature",
      "relative_humidity",
      "wind",
      "precipitation_flux"))
} # cfmet.downscale.daily


##' Get time series vector from netCDF file
##'
##' internal convenience function for
##' streamlining extraction of data from netCDF files
##' with CF-compliant variable names
##'
##' @param var name of variable to extract
##' @param lati,loni latitude and longitude to extract
##' @param run.dates data frame of dates to read
##' @param met.nc netcdf file with CF variable names
##'
##' @return numeric vector
##' @export
##' @author David Shaner LeBauer
get.ncvector <- function(var, lati = lati, loni = loni, run.dates = run.dates, met.nc) {

  start.idx <- c(latitude = lati, longitude = loni, time = run.dates$index[1])
  count.idx <- c(latitude = 1, longitude = 1, time = nrow(run.dates))
  dim.order <- sapply(met.nc$var$air_temperature$dim, function(x) x$name)
  ncvar_get2 <- function(var) {
    ans <- ncdf4::ncvar_get(nc = met.nc, varid = var, start = start.idx[dim.order], count = count.idx[dim.order])
    return(as.numeric(ans))
  } # ncvar_get2

  if (var %in% attributes(met.nc$var)$names) {
    ans <- ncvar_get2(var)
  } else if (var == "air_pressure") {
    ans <- 1013.25
  } else if (var == "wind") {
    ans <- sqrt(ncvar_get2("northward_wind")^2 + ncvar_get2("eastward_wind")^2)
  } else {
    ans <- NULL
  }

  if (var == "precipitation_flux") {
    precip_units <- met.nc$var[["precipitation_flux"]]$units
    precip_units <- gsub("kg m-2", "mm", precip_units)
    precip_units <- gsub("kg/m2", "mm", precip_units)
    precip_units <- gsub("kg/m\\^2", "mm", precip_units)
    ans <- PEcAn.utils::ud_convert(ans, precip_units, "mm s-1")
  }
  return(ans)
} # cfmet.downscale.time
