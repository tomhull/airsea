% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/solubility.R
\name{Csat}
\alias{Csat}
\title{Equilibrium Oxygen saturation concentration}
\usage{
Csat(temp, salinity)
}
\arguments{
\item{temp}{numeric vector of water temperature in degrees Celsius}

\item{salinity}{numeric vector of salinity (PSU)}
}
\value{
vector of saturation concentration in mmol m-3
}
\description{
Calculates oxygen saturation concentration in equilibrium with the atmosphere
as per Garcia & Gordon, 1992 (Benson & Kraus data)
}
\details{
TODO
}
\examples{
Csat(10, 35)  # saturation concentration at 10 degrees and 35 salinity
}
\keyword{oxygen}

