#' transform data to O2NCP format
#'
#' @details TODO
#' expects column headings to match the following format (order does not matter):
#' dateTime = POSIXct, T = temperature (oC), S = salinity, C = oxygen concentration (mmol m-3), u = wind speed (m s-1),
#' Pslp = pressure at sea level (mbar), h = mixed layer depth (m).
#' 
#' @param x data frame of observations
#'
#' @return data frame in correct format for O2NCP
#' @export
transform_data <- function(x){
    # transforms timeseries data to inital and post variables in wide format, i.e. like a lhs
    dat = with(x,
               data.frame(dateTime = dateTime,
                          T0 = T, S0 = S,
                          C0 = C, u0 = u,
                          h0 = h, Pslp0 = Pslp))
    dat$T1 = dat$T0 + c(diff(dat$T0), 0)
    dat$S1 = dat$S0 + c(diff(dat$S0), 0)
    dat$u1 = dat$u0 + c(diff(dat$u0), 0)
    dat$C1 = dat$C0 + c(diff(dat$C0), 0)
    dat$h1 = dat$h0 + c(diff(dat$h0), 0)
    dat$Pslp1 = dat$Pslp0 + c(diff(dat$Pslp0), 0)
    dat$timePeriod = c(diff(as.numeric(dat$dateTime)), 0)
    dat$end = dat$timePeriod + dat$dateTime
    return(dat[-nrow(dat),])
}


#' Correct wind speed to 10m height
#' 
#' @details TODO
#'
#' @param U numerical vector of wind speed
#' @param z height at x
#' @param method character string of either "liutang", "COARE" or "hsu". Default is "liutang"
#' @return u10
#' @references Liu W, Tang W. Equivalent neutral wind. Pasadena: Jet Propulsion Laboratory. 1996. Available: http://128.149.33.88/publication/paper/Liu-Tang-1996-jpl.pdf
#' @references COARE formulation courtesy of Tom Bell (Personal communication). 2013.
#' @references Hsu SA, Meindl EA, Gilhousen DB. Determining the Power-Law Wind-Profile Exponent under Near-Neutral Stability Conditions at Sea. J Appl Meteorol. 1994;33: 757-765
#' @export
convert_u10 <- function(U, z, method = 'liutang'){
  if(z < 10){stop("formulations only work for z > 10")}
  #
  liutang <- function(U, z){
    U10 = U / (1 + 2.5 * sqrt(dragCoef(U) * log(z / 10) ))
    return(U10)
  }
  COARE <- function(U, z){
    u10 = U * log(10 / 1E-4) / log(z / 1E-4)
    return(u10)
  }
  hsu <- function(U, z){
    u10 = U * (10 / z)^0.11
    return(u10)
  }
    switch(method,
           liutang = liutang(U, z),
           COARE = COARE(U, z),
           hsu = hsu(U, z)
           )
}