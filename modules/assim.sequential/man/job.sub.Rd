% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sda.qsub.functions.R
\name{job.sub}
\alias{job.sub}
\title{job.sub}
\usage{
job.sub(settings, rm.file, only.nc)
}
\arguments{
\item{settings}{PEcAn settings object.}

\item{rm.file}{Boolean: if we want to delete files.
the model execution will be submitted if it's `FALSE`.}

\item{only.nc}{Boolean: if we want to only delete NC files.}
}
\description{
This function provides an option for executing model or removing files using cluster.
}
\author{
Dongchen Zhang.
}
