#' air side gas transfer velocity
#'
#' Calculates air side gas transfer velocity (ka).
#'
#' @details 
#' By default implements the  gas transfer velocity parametrisation proposed by Johnson 2010, which is a modified form of that presented by Jeffrey et al., 2010 (Ja10)
#' Other options are:
#' Liss 1973, windspeed only parameterisation (L73)
#' Duce et al, 1991 (D91), theoretical ka 
#' Mackay and Yeun 1983 (MY83), empirical fit to wind tunnel study
#' Jeffrey et al 2010 (JE10), unmodified, a parameterisation of NOAA COARE
#' 
#' For full details of parameterisations, choice of drag coefficient and diffusivity calculations see Johnson 2010

#' @param compound character string of compound of interest
#' @param T vector o temperature in degrees Centigrade
#' @param u vector of wind speed in meters per second at 10 meters height
#' @param method optional string determining kw parameterisation to use, see XXX for list default is 'JO10' (Johnson 2010)
#' @return vector of gas transfer velocity in meters per second
#' @keywords ka air side  transfer velocity gas exchange
#' @references Johnson, M. T. A numerical scheme to calculate temperature and salinity dependent air-water transfer velocities for any gas. Ocean Sci. Discuss. 7, 251-290 (2010).
#' @references Jeffery, C., Robinson, I., and Woolf, D.: Tuning a physically-based model of the air-sea gas transfer velocity, Ocean Modell., 31,28–35, doi:10.1016/j.ocemod.2009.09.001, 2010..
#' @references Liss, P.: Processes of gas exchange across an air-water interface, Deep Sea Res., 20, 221–238, doi:10.1016/0011-7471(73)90013-2, 1973.
#' @references Mackay, D. and Yeun, A. T. K.: Mass transfer coefficient correlations for volatilization of organic solutes from water, Environ. Sci. Technol., 17, 211–217, doi:10.1021/es00110a006, 1983
#' @references Duce, R. A., Liss, P., Merrill, J. T., Atlas, E. L., Buat-Menard, P.,Hicks, B. B., Miller, J. M., Prospero, J. M., Arimoto, R., Church, T. M., Ellis, W., Galloway, J. N., Hansen, L., Kickells, T., Knap, A. H., Reinhardt, K. H., Schneider, B., Soudine, A., Tokos, J. J., Tsunogai, S., Wollast, R., and Zhou, M.: The atmospheric input of trace species to the World ocean, Global Biogeochem. Cy., 5, 193–259, 1991
#' @examples
#' kw('O2', 10, 7, 35)  # gas transfer velocity for oxygen at 10oC, 7 m-1 s-1 winds and 35 salinity.
#' @export

ka <- function(compound, T, u, method = 'JO10'){
  Liss1973 <- function(u){
    #calculate ka according to Liss 1973 (non compound-specific)
    (0.005+(0.21*u))/100
  }
  
  
  MackayYeun1983 <- function(compound,u,T)
    #ka according to Mackay and Yeun 1983
    1e-3 + (46.2e-5*sqrt(6.1+(0.63*u))*u*(Sc_air(compound,T))^-0.67)
  
  Duce1991 <- function(compound,u){
    #calculate ka gas phase transfer velocity according to Duce et al 1991
    u/(770+(45*((compounds[compound,"mw"])^(1/3))))
  }
  
  Jeffrey2010<-function(compound,u,T){
    #using smith(1980) Cd
    von_kar<-0.4
    Cd<-(1e-4*(6.1+0.63*u))
    Sc<-Sc_air(compound,T)
    ra<-13.3*sqrt(Sc) + (Cd^(-0.5)) - 5 + log(Sc)/(2*von_kar)
    u_star<-u*sqrt(C_D(u))
    u_star/ra
  }
  
  Johnson2010<-function(compound,u,T){
    (1e-3+Jeffrey_ka(compound,u,T))
  }
  
  switch(method,
         JO10 = Johnson2010(compound,u,T),
         L73 = Liss73(u),
         MY83 = MackayYeun1983(compound,u,T),
         D91 = Duce1991(compound,u),
         JE10 = Jeffrey2010(compound,u,t)
  
}


