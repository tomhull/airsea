#' Equilibrium bubble supersaturation
#'
#' Implements the bubble supersaturation parametrisation of Woolf and Thorpe (1991)
#' Note: parametrisation is only applicable for Oxygen (O2), Nitrogen (N2), Argon (Ar) and Carbon Dioxide (CO2).
#'
#' @details TODO
#' @param compound character string of compound of interest ('O2' ,'N2', 'Ar' or 'CO2')
#' @param u vector of wind speed in meters per second at 10 meters height
#' @return vector of equilibrium bubble supersaturation effects (no units)
#' @keywords kw gas exchange bubble
#' @references Woolf, D. K., & Thorpe, S. A. (1991). Bubbles and the air-sea exchange of gases in near-saturation conditions. Journal of Marine Research, 49(3), 435-466. doi:10.1357/002224091784995765
#' @examples
#' bubbleSat('Ar', 7)  # bubble supersaturation effect for Argon at 7 m-1 s-1 winds
#' @export
bubbleSat <- function(compound = 'O2', u){
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


#' Drag Coefficient
#'
#' Calculates the drag coefficent
#' @details TODO, discuss methods
#' @param u vector of wind speed in meters per second at 10 meters height
#' @param method optional string determining parametrisation, valid strings are "satLargePond", "smith". default is "satLargePond"
#' @return vector of drag coefficent
#' @references Sullivan, P. P., Romero, L., McWilliams, J. C., & Melville, W. K. (2012). Transient Evolution of Langmuir Turbulence in Ocean Boundary Layers Driven by Hurricane Winds and Waves. Journal of Physical Oceanography, 42(11), 1959-1980. doi:10.1175/JPO-D-12-025.1
#' @references Smith, S. D. (1980). Wind stress and heat flux over the ocean in gale force winds. J Phys. Oceanogr. 10, 709-726.
#' @keywords drag coefficent
#' @export
dragCoef <- function(u, method = 'satLargePond'){
    # Saturated Large and Pond (1981) Drag Coefficent
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

    SmithCD <- function(u){
    # Smith (1980) Drag Coefficient
        1e-4*(6.1+0.63 * u)
    }

    switch(method,
           satLargePond = LargePondCD(u),
           smith = SmithCD(u)
           )
}