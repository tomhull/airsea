#' Wanninkhov et al (2009) Gas transfer velocity
#'
#' Implements the gas transfer velocity parametrisation of Wanninkhov et al (2009).
#'
#' @details TODO
#' @param compound character string of compound of interest
#' @param T vector of temperature in degrees Centigrade
#' @param u vector of wind speed in meters per second at 10 meters height
#' @param S vector of salinity
#' @return vector of gas transfer velocity in meters per second
#' @keywords kw gas exchange
#' @references Wanninkhof, R., Asher, W. E., Ho, D. T., Sweeney, C., & McGillis, W. R. (2009).  'Advances in Quantifying Air-Sea Gas Exchange and Environmental Forcing*. Annual Review of Marine Science, 1(1), 213–244. doi:10.1146/annurev.marine.010908.163742
#' @examples
#' Wann09Kw('O2', 10, 7, 35)  # gas transfer velocity for oxygen
Wann09Kw <- function(compound, T, u, S, normalize=0, schmidt_formulation=0){
    # Kw parametrisation of Wanninkhof2009
    # output in meters per second
    if (normalize!=0) schmidt_number <- normalize else
        schmidt_number <- schmidt_out(compound, T, S, schmidt_formulation)
    (3 + (0.1 * u) + (0.064 * u^2) + (0.011 * u^3)) * ((schmidt_number / 660)^(-0.5)) / (100 * 3600)
}

#' Woolf and Thorpe (1991) equilibrium bubble supersaturation
#'
#' Implements the bubble supersaturation parametrisation of Woolf and Thorpe (1991)
#' Note: parametrisation is only applicable for Oxygen (O2), Nitrogen (N2), Argon (Ar) and Carbon Dioxide (CO2).
#'
#' @details TODO
#' @param compound character string of compound of interest ('O2' ,'N2', 'Ar', 'CO2')
#' @param u vector of wind speed in meters per second at 10 meters height
#' @return vector of equilibrium bubble supersaturation effects (no units)
#' @keywords kw gas exchange bubble
#' @references Woolf, D. K., & Thorpe, S. A. (1991). Bubbles and the air-sea exchange of gases in near-saturation conditions. Journal of Marine Research, 49(3), 435–466. doi:10.1357/002224091784995765
#' @examples
#' WoolfB('Ar', 7)  # bubble supersaturation effect for Argon at 7 m/s winds
WoolfB <- function(compound = 'O2', u){
    # bubble saturation effects from Woolf1991
    # 9m/s/ = 1% supersat
    U = NA
    if(compound == 'N2'){U = 7.2}
    if(compound == 'CO2'){U = 49}
    if(compound == 'Ar'){U = 9.6}
    if(compound == 'O2'){U = 9}
    B = 0.01 * (u / U)^2
    # B is a scaling term and has no units
    return(B)
}

#' Smith (1980) Drag Coefficient
#'
#' Calculates the drag coefficent of Smith (1980)
#' @details TODO
#' @param u vector of wind speed in meters per second at 10 meters height
#' @return vector of drag coefficent
#' @references Smith, S. D. (1980). Wind stress and heat flux over the ocean in gale force winds. J Phys. Oceanogr. 10, 709-726.
#' @keywords drag coefficent
SmithCD <- function(u){
    1e-4*(6.1+0.63*u)
}

#' Saturated Large and Pond (1981) Drag Coefficent
#'
#' Calculates the drag coefficent of Large and Pond with the high wind saturation as per Sullivan et al, (2012)
#' @details TODO
#' @param u vector of wind speed in meters per second at 10 meters height
#' @return vector of drag coefficent
#' @references Sullivan, P. P., Romero, L., McWilliams, J. C., & Melville, W. K. (2012). Transient Evolution of Langmuir Turbulence in Ocean Boundary Layers Driven by Hurricane Winds and Waves. Journal of Physical Oceanography, 42(11), 1959–1980. doi:10.1175/JPO-D-12-025.1
#' @keywords drag coefficent
LargePondCD <- function(u){
    # returns drag coefficent as per Large & Pond, 1981
    # with high wind saturation as per Sullivan et al, 2012
    cd <- function(x){
               if(is.na(x)){return(NA)}
               if(x < 11){return(0.0012)}
               if(x > 20){return(1.8 * 10^-3)}
               else{return((0.49 + 0.065 * x) * 10^-3)}
    }
    sapply(u, cd)
}
