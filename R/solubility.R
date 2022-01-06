#' Equilibrium Oxygen saturation concentration
#'
#' Calculates oxygen saturation concentration (C_sat) in equilibrium with the atmosphere
#' as per Garcia & Gordon, 1992, using the Benson and Kraus data.
#'
#' This is also known as C* in some literature.
#'
#' conversions via SCOR WG 142
#'
#' check values from Garcia and Gordon
#' 10oC 35 salinity
#' 6.315 ml/l
#' 274.610 umol kg-1
#'
#' @param temp numeric vector of water temperature in degrees Celsius
#' @param salinity numeric vector of salinity (PSU)
#' @param unit "mmolm" for mmol m-3 (default), "mgl" for mg l-1 or "umolkg" for umol kg-1.
#' @param p_atm atmospheric (air) pressure in hPa (default = 1013.25)
#' @return vector of saturation concentration in mmol m-3
#' @keywords oxygen
#' @references
#' 1. Bittig, H., Körtzinger, A., Johnson, K., Claustre, H., Emerson, S., Fennel, K., Garcia, H., Gilbert, D., Gruber, N., Kang, D.-J., Naqvi, W., Prakash, S., Riser, S., Thierry, V., Tilbrook, B., Uchida, H., Ulloa, O., Xing, X., 2018. SCOR WG 142: Quality Control Procedures for Oxygen and Other Biogeochemical Sensors on Floats and Gliders. Recommendations on the conversion between oxygen quantities for Bio-Argo floats and other autonomous sensor platforms. <https://doi.org/10/ggzjj3>
#' 1. Garcia, H.E., Gordon, L.I., 1992. Oxygen solubility in seawater: Better fitting equations. Limnology and Oceanography 37, 1307–1312. <https://doi.org/10/dxf339>
#' @examples
#' oxygen.sat(10, 35)  # 282.015 mmol m-3, saturation concentration at 10 degrees and 35 salinity
#' oxygen.sat(10, 35, "mll")  # 6.314767 ml/l kg-1 saturation concentration at 10 degrees and 35 salinity
#' oxygen.sat(10, 35, "umolkg")  # 274.6095 umol kg-1 saturation concentration at 10 degrees and 35 salinity
#' oxygen.sat(20, 35, p_atm = 1020.5)  # 284.06 true saturation concentration at 20 degrees and 35 salinity in mmol m-3 when local air pressure is 1020.5 hPa
#' @export
Csat <- function(temp, salinity, unit = "mmolm", p_atm = 1013.25){
  
  if(unit == "umolkg"){
    # umol kg coefficents
    A0 = 5.80871;
    A1 = 3.20291;
    A2 = 4.17887;
    A3 = 5.10006;
    A4 = -9.86643-2;
    A5 = 3.80369;
    B0 = -7.01577e-3;
    B1 = -7.70028e-3;
    B2 = -1.13864e-2;
    B3 = -9.51519e-3;
    C0 = -2.75915e-7;
  }else{
    # cm3 dm-3 coefficents (ml/l)
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
  }
  Ts = log((298.15-temp)/(273.15+temp))
  
  Csat = A0+(A1*Ts)+(A2*Ts^2)+
    (A3*Ts^3)+(A4*Ts^4)+(A5*Ts^5)+
    salinity*(B0+(B1*Ts)+(B2*Ts^2)+(B3*Ts^3))+
    (C0*salinity^2)
  
  # adjust for in-situ pressure # as per SCOR WG 142
  # using vapour pressure correction, but still always at 1 atm so ignore the exp(Vm*P...
  pH2Osat = 1013.25 * (exp(24.4543-(67.4509*(100/(temp+273.15)))-(4.8489*log(((273.15+temp)/100)))-0.000544*salinity)) # saturated water vapour in mbar
  Csat = exp(Csat) * (p_atm-pH2Osat)/(1013.25-pH2Osat)
  
  if(unit == "mmolm"){
    return(Csat * 44.6596)     # convert ml/l to mmol m-3  as per SCOR WG 142
  }
  if(unit == "mll"){
    return(Csat) # no conversion
  }
  if(unit == "mgl"){
    return(Csat / 0.699745)     # convert ml/l to mg/l
  }
  if(unit == "umolkg"){
    return(Csat) # no conversion
  }
  
  # 1 μmol O2 = .022391 ml at sea surface pressure
  # 1 mg/l = 22.391 ml/31.998 = 0.699745 ml/l
  
  else{
    stop("unit not recognised")
  }
}


#' Equilibrium Oxygen saturation concentration (combined fit)
#'
#' Calculates oxygen saturation concentration in equilibrium with the atmosphere
#' as per Garcia & Gordon, 1992, using the combined fit, which is not recommended.
#' This is however the formulation used on Aanderaa optodes.
#'
#' @param temp numeric vector of water temperature in degrees Celsius
#' @param salinity numeric vector of salinity (PSU)
#' @param unit "mmolm" for mmol m-3 (default), "mgl" for mg l-1 or "umolkg" for umol kg-1.
#' @return vector of saturation concentration in mmol m-3
#' @keywords oxygen
#'
Csat.combined <- function(temp, salinity, unit = "mmolm"){
  
  if(unit == "umolkg"){
    # umol kg coefficents
    A0 = 5.80818
    A1 = 3.20684
    A2 = 4.11890
    A3 = 4.93845
    A4 = 1.01567
    A5 = 1.41575
    B0 = -7.01211e-3
    B1 = -7.25958e-3
    B2 = -7.93334e-3
    B3 = -5.54491e-3
    C0 = -1.32412e-7
  }else{
    # cm3 dm-3 coefficents (ml/l)
    A0 = 2.00856
    A1 = 3.22400
    A2 = 3.99063
    A3 = 4.80299
    A4 = 9.78188e-1
    A5 = 1.71069
    B0 = -6.24097e-3
    B1 = -6.93498e-3
    B2 = -6.90358e-3
    B3 = -4.29155e-3
    C0 = -3.11680e-7
  }
  Ts = log((298.15-temp)/(273.15+temp))
  
  Csat = A0+(A1*Ts)+(A2*Ts^2)+
    (A3*Ts^3)+(A4*Ts^4)+(A5*Ts^5)+
    salinity*(B0+(B1*Ts)+(B2*Ts^2)+(B3*Ts^3))+
    (C0*salinity^2)
  
  if(unit == "mmolm"){
    return(exp(Csat) * 44.6596)     # convert ml/l to mmol m-3  as per SCOR WG 142
  }
  if(unit == "mll"){
    return(exp(Csat)) # no conversion
  }
  if(unit == "mgl"){
    return(exp(Csat) / 0.699745)     # convert ml/l to mg/l
  }
  if(unit == "umolkg"){
    return(exp(Csat)) # no conversion
  }
  
  else{
    stop("unit not recognised")
  }
}


#' Gas solubility (Henry's law air-water partioning coefficients)
#' 
#'
#' @param compound e.g. "O2"
#' @param T temperature in Celcius
#' @param S salinity
#'
#' @details Applies polynomial fits to empirical data where available (O2,CO2,N2O etc), or where this data is not available, uses the synthesis of Henry's law contstants and temperature dependencies along with the salinity relationship from Johnson (2010) which in turn uses Rolf Sander's collection of Henry's law constants
#' @details Where available, gas-specific polynomial fits following the model of Weiss 1970 (equation 4) should be used in preference to the Johnson/Sander scheme. This is partcularly important for long-lived (i.e. close-to equilibrium) gases such as N2, O2, Ar, CO2 and N2O. For shorter lived / far from equilibrium gases the simpler scheme will generally to a good job.
#' @references Johnson, M. T. A numerical scheme to calculate temperature and salinity dependent air-water transfer velocities for any gas. Ocean Sci. Discuss. 7, 251-290 (2010).
#' @references Sander, R. Compilation of Henry's Law constants for inorganic and organic species of potential importance in environmental chemistry (Version 3), http://www.henrys-law.org, 1999.
#' @references Weiss, R, The solubility of nitrogen, oxygen and argon in water and seawater, Deep Sea Research (17), 721-735, 1970
#' @return partitioning coefficient in moles/litre/atmosphere
#' @export
#'
#' 
GasSolubility<-function(compound,T, S){
  
  # Johnson_Sander_GasSolubility<-function(compound,T,S)
    KH_Molar_per_atmosphere <- function(compound,T,S){
      # calculate the Henry's law constant in M/atm from Sander data
      # applying the salting out factor from Johnson 2010
      (compounds[compound,"KH"]*exp((compounds[compound,"tVar"])*((1/(T+273.15))-(1/298.15))))/K_H_factor(compound, S)
    }
  
  # TODO add marelac data method
  return(KH_Molar_per_atmosphere(compound, T, S))
}
  

###############################################################################################
### Under the hood of the Johnson method ######################################################
##############################################################################################

  #Calculate Henry's law constant at a given T in pure water according to Sander (1999)
KH0 <- function(compound,T=25){
  12.2/((273.15+T)*(compounds[compound,"KH"])*exp((compounds[compound,"tVar"])*((1/(T+273.15))-(1/298.15))))
}

## Salting out (or in for very small compounds) factor from empirical parameterisation by Johnson 2010
Ks <- function(compound){
  theta = (7.3353282561828962e-04 + (3.3961477466551352e-05*log(KH0(compound))) + (-2.4088830102075734E-06*(log(KH0(compound)))^2) + (1.5711393120941302E-07*(log(KH0(compound)))^3))
  theta*log(Vb(compound))
}

K_H_factor <-function(compound, S){
  #calculate salinity-dependent salting-out scaling factor for KH0 (see manuscript)
  10^(Ks(compound)*S)
}

