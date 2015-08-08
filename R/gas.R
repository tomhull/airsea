# gas transfer velocities

#' water side gas transfer velocity
#'
#' Calculates water side gas transfer velocity (kw).
#'
#' @details TODO - K_calcs scheme
#' By default implements the gas transfer velocity parametrisation of Wanninkhov, 2014 (WA14).
#' utilising XXX'
#' Other options are:
#' Wanninkhov et al, 2009 (WA09),
#' Nighingale et al, 2000 empirical fit to dual tracer data (NG00).
#' Wanninkhov, 1992 (WA92).
#' set normalize to schmidt number value e.g. 660 to compare with other kw curves.
#' @param compound character string of compound of interest
#' @param T vector of temperature in degrees Centigrade
#' @param u vector of wind speed in meters per second at 10 meters height
#' @param S vector of salinity
#' @param method optional string determining kw parameterisation to use, see XXX for list default is 'WA09' (Wanninkhof et al, 2009)
#' @param normalize optional integer Schmitt number, if not = 0 kw values are normalised to specified Schmitt number. Default is 0
#' @param schmidt_method, see `Sch` documentation for details, Default is 'mean'
#' @return vector of gas transfer velocity in meters per second
#' @keywords kw gas exchange transfer velocity
#' @references Johnson, M. T. A numerical scheme to calculate temperature and salinity dependent air-water transfer velocities for any gas. Ocean Sci. Discuss. 7, 251-290 (2010).
#' @references Wanninkhof, R. Relationship between wind speed and gas exchange over the ocean revisited. Limnol. Oceanogr. Methods 12, 351-362 (2014).
#' @references Wanninkhof, R., Asher, W. E., Ho, D. T., Sweeney, C. & McGillis, W. R. Advances in Quantifying Air-Sea Gas Exchange and Environmental Forcing. Ann. Rev. Mar. Sci. 1, 213-244 (2009).
#' @references Nightingale, P. D. et al. In situ evaluation of air-sea gas exchange parameterizations using novel conservative and volatile tracers. Global Biogeochem. Cycles 14, 373-387 (2000).
#' @references Wanninkhof, R. Relationship between wind speed and gas exchange over the ocean. J. Geophys. Res. 97, 7373 (1992).
#' @examples
#' kw('O2', 10, 7, 35)  # gas transfer velocity for oxygen at 10oC, 7 m-1 s-1 winds and 35 salinity.
#' @export

kw <- function(compound, T, u, S, method = 'WA14', normalize = 0, schmidt_method = 'mean'){
    Wann09 <- function(compound, T, u, S, normalize=0, schmidt_method){
        # Kw parametrisation of Wanninkhof2009
        # output in meters per second
        if (normalize!=0) schmidt_number <- normalize else
            schmidt_number <- Sch(compound, T, S, schmidt_method)
        (3 + (0.1 * u) + (0.064 * u^2) + (0.011 * u^3)) * ((schmidt_number / 660)^(-0.5)) / (100 * 3600)
    }
    Wann14 <- function(compound, T, u, S, normalize=0, schmidt_method){
        # Kw parametrisation of Wanninkhof2014
        # output in meters per second
        if (normalize!=0) schmidt_number <- normalize else
            schmidt_number <- Sch(compound, T, S, schmidt_method)
        (0.251 * u^2) * ((schmidt_number / 660)^(-0.5)) / (100 * 3600)
    }
    Nightingale00 <- function(compound, T, u, S, normalize, schmidt_method){
        # empirical fit to dual tracer data by Nightingale et al 2000 (GBC)
        # note k600 not k660 for their study
        if (normalize!=0) schmidt_number <- normalize else
            schmidt_number <- Sch(compound, T, S, schmidt_method)
        (((0.222*u^2)+0.333*u)*(schmidt_number/600)^(-0.5))/(100*3600)
    }
    Wann92 <- function (compound, T, u, S, normalize, schmidt_method){
        #calculate kw transfer velocity in m/s according to Wanninkhof 1992
        if (normalize!=0) schmidt_number <- normalize else
            schmidt_number <- Sch(compound, T, S, schmidt_method)
        (0.31*u^2*((schmidt_number/660)^(-0.5)))/(100*3600)
    }
    switch(method,
           WA09 = Wann09(compound, T, u, S, normalize, schmidt_method),
           WA14 = Wann14(compound, T, u, S, normalize, schmidt_method),
           NG00 = Nightingale00(compound, T, u, S, normalize, schmidt_method),
           WA92 = Wann92(compound, T, u, S, normalize, schmidt_method)
           )
}


#' Woolf and Thorpe (1991) equilibrium bubble supersaturation
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
           smith = SmithCD(u),
           )
}


#' Equilibrium Oxygen saturation concentration
#'
#' Calculates oxygen saturation concentration in equilibrium with the atmosphere
#' as per Garcia & Gordon, 1992 (Benson & Kraus data)
#' @details TODO
#' @param temp numeric vector of water temperature in degrees Celsius
#' @param salinity numeric vector of salinity (PSU)
#' @return vector of saturation concentration in mmol m-3
#' @keywords oxygen
#' @examples
#' Csat(10, 35)  # saturation concentration at 10 degrees and 35 salinity
#' @export
Csat <- function(temp, salinity){
    # Coefficents
    A0 = 2.00907
    A1 = 3.22014
    A2 = 4.05010
    A3 = 4.94457
    A4 = -0.256847
    A5 = 3.88767
    B0 = -0.00624523
    B1 = -0.00737614
    B2 = -0.0103410
    B3 = -0.00817083
    C0 = -4.88682E-07

    Ts = log((298.15-temp)/(273.15+temp))

    O2.sat = A0+(A1*Ts)+(A2*Ts^2)+
    (A3*Ts^3)+(A4*Ts^4)+(A5*Ts^5)+
    salinity*(B0+(B1*Ts)+(B2*Ts^2)+(B3*Ts^3))+
    (C0*salinity^2)

    return((exp(O2.sat))* 44.6608)     # output in uMol/L
}
