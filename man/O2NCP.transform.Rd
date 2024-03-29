% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ncp.R
\name{O2NCP.transform}
\alias{O2NCP.transform}
\title{transform data to O2NCP format}
\usage{
O2NCP.transform(x)
}
\arguments{
\item{x}{data frame of observations}
}
\value{
data frame in correct format for O2NCP
}
\description{
transform data to O2NCP format
}
\details{
Function to take a typical time series of observations and convert to a start and end type data.frame for use with `O2NCP.mean`.
expects column headings to match the following format (order does not matter):
dateTime = POSIXct, `T` = temperature (oC), `S` = salinity, `C` = oxygen concentration (mmol m-3), `u` = wind speed (m s-1),
`Pslp` = pressure at sea level (mbar), `h` = mixed layer depth (m).
Optionally `Cb` = bottom oxygen (mmol m-3) if entrainment is to be calculated.
Optionally `hmin` for minimum mld for uncertainty analysis
}
