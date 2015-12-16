#' water side gas transfer velocity
#'
#' Calculates water side gas transfer velocity (kw).
#'
#' @details TODO - K_calcs scheme
#' By default implements the gas transfer velocity parametrisation of Wanninkhov, 2014 (WA14).
#' Other options are:
#' Wanninkhov et al, 2009 (WA09), ?CO2 only
#' Nighingale et al, 2000 empirical fit to dual tracer data (NG00).
#' Wanninkhov, 1992 (WA92).
#' Liss and Merlivat, 1983 (LM83).
#' set normalize to schmidt number value e.g. 660 to compare with other kw curves.
#' @param compound character string of compound of interest
#' @param T vector of temperature in degrees Centigrade
#' @param u vector of wind speed in meters per second at 10 meters height
#' @param S vector of salinity
#' @param method optional string determining kw parameterisation to use, see details for list default is 'WA14' (Wanninkhof et al, 2009)
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
    lissmer83 <- function (compound, T, u, S, normalize, schmidt_method){
        #calculate kw transfer velocity in m/s according to Liss and Merlivat, 1983
    		k600<-ifelse(u<3.6,0.17*u,
     			ifelse(u<13,(2.85*u)-9.65,(5.9*u)-49.3))
     	    	if (normalize!=0) schmidt_number<-normalize else
           schmidt_number <- Sch(compound, T, S, schmidt_method)
     		ifelse(u<3.6,(k600*((schmidt_number/600)^(-0.66)))/(100*3600),(k600*((schmidt_number/600)^(-0.5)))/(100*3600))
    }
    switch(method,
           WA09 = Wann09(compound, T, u, S, normalize, schmidt_method),
           WA14 = Wann14(compound, T, u, S, normalize, schmidt_method),
           NG00 = Nightingale00(compound, T, u, S, normalize, schmidt_method),
           LM83 = lissmer83(compound, T, u, S, normalize, schmidt_method),
           WA92 = Wann92(compound, T, u, S, normalize, schmidt_method)
           )
}


