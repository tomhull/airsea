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