library(purrr)
library(foreach)
setwd("/projectnb/dietzelab/dongchen/anchorSites/NA_runs/downscale_maps/")
# average ERA5 to climatic covariates.
outdir <- "/projectnb/dietzelab/dongchen/anchorSites/NA_runs/GridMET"
in.path <- "/projectnb/dietzelab/dongchen/anchorSites/ERA5/"
dates <- c(as.Date("2012-01-01"), seq(as.Date("2012-07-15"), as.Date("2021-07-15"), "1 year"))
start.dates <- dates[1:10]
end.dates <- dates[2:11]
paths <- c()
for (i in 1:10) {
  paths <- c(paths, PEcAnAssimSequential:::Average.ERA5.2.GeoTIFF(start.dates[i], end.dates[i], in.path, outdir))
  print(i)
}
# setup.
base.map.dir <- "/projectnb/dietzelab/dongchen/anchorSites/downscale/MODIS_NLCD_LC.tif"
load("/projectnb/dietzelab/dongchen/anchorSites/NA_runs/SDA_25ens_2024_11_25/sda.all.forecast.analysis.Rdata")
variables <- c("AbvGrndWood", "LAI", "SoilMoistFrac", "TotSoilCarb")
settings <- "/projectnb/dietzelab/dongchen/anchorSites/NA_runs/SDA/pecanIC.xml"
outdir <- "/projectnb/dietzelab/dongchen/anchorSites/NA_runs/downscale_maps/"
cores <- 28
date <- seq(as.Date("2012-07-15"), as.Date("2021-07-15"), "1 year")
# loop over years.
for (i in seq_along(date)) {
  # setup covariates paths and variable names.
  cov.tif.file.list <- list(LC = list(dir = "/projectnb/dietzelab/dongchen/anchorSites/downscale/MODIS_NLCD_LC.tif",
                                      var.name = "LC"),
                            year_since_disturb = list(dir = "/projectnb/dietzelab/dongchen/anchorSites/downscale/MODIS_LC/outputs/age.tif",
                                                      var.name = "year_since_disturb"),
                            agb = list(dir = "/projectnb/dietzelab/dongchen/anchorSites/downscale/AGB/agb.tif",
                                       var.name = "agb"),
                            twi = list(dir = "/projectnb/dietzelab/dongchen/anchorSites/downscale/TWI/TWI_resample.tiff",
                                       var.name = "twi"),
                            met = list(dir = paths[i],
                                       var.name = c("temp", "prec", "srad", "vapr")),
                            soil = list(dir = "/projectnb/dietzelab/dongchen/anchorSites/downscale/SoilGrids.tif",
                                        var.name = c("PH", "N", "SOC", "Sand")))
  # Assemble covariates.
  if (file.exists(paste0(outdir, "covariates_", lubridate::year(date[i]), ".tiff"))) {
    covariates.dir <- paste0(outdir, "covariates_", lubridate::year(date[i]), ".tiff")
  } else {
    covariates.dir <- create.covariates.geotiff(outdir = outdir, 
                                                year = lubridate::year(date[i]), 
                                                base.map.dir = base.map.dir, 
                                                cov.tif.file.list = cov.tif.file.list, 
                                                normalize = T, 
                                                cores = cores)
  }
  # grab analysis.
  analysis.yr <- analysis.all[[i]]
  time <- date[i]
  # loop over carbon types.
  for (j in seq_along(variables)) {
    # setup folder.
    variable <- variables[j]
    folder.path <- file.path(outdir, paste0(variables[j], "_", date[i]))
    dir.create(folder.path)
    save(list = c("settings", "analysis.yr", "covariates.dir", "time", "variable", "folder.path", "base.map.dir", "cores", "outdir"),
         file = file.path(folder.path, "dat.Rdata"))
    # prepare for qsub.
    jobsh <- c("#!/bin/bash -l", 
               "module load R/4.1.2", 
               "echo \"require (PEcAnAssimSequential)", 
               "      require (foreach)",
               "      require (purrr)",
               "      downscale.qsub.main('@FOLDER_PATH@')", 
               "    \" | R --no-save")
    jobsh <- gsub("@FOLDER_PATH@", folder.path, jobsh)
    writeLines(jobsh, con = file.path(folder.path, "job.sh"))
    # qsub command.
    qsub <- "qsub -l h_rt=6:00:00 -l buyin -pe omp @CORES@ -V -N @NAME@ -o @STDOUT@ -e @STDERR@ -S /bin/bash"
    qsub <- gsub("@CORES@", cores, qsub)
    qsub <- gsub("@NAME@", paste0("ds_", i, "_", j), qsub)
    qsub <- gsub("@STDOUT@", file.path(folder.path, "stdout.log"), qsub)
    qsub <- gsub("@STDERR@", file.path(folder.path, "stderr.log"), qsub)
    qsub <- strsplit(qsub, " (?=([^\"']*\"[^\"']*\")*[^\"']*$)", perl = TRUE)
    cmd <- qsub[[1]]
    out <- system2(cmd, file.path(folder.path, "job.sh"), stdout = TRUE, stderr = TRUE)
  }
}
