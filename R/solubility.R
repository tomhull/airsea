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
    # TODO refer to solubility for non - "O2"
}


#' Gas solubility (Henry's law air-water partioning coefficients)
#' 
#'
#' @param compound 
#' @param S salinity
#' @param T temperature in Celcius
#'
#' @details Applies polynomial fits to empirical data where available (O2,CO2,N2O etc), or where this data is not available, uses the synthesis of Henry's law contstants and temperature dependencies along with the salinity relationship from Johnson (2010) which in turn uses Rolf Sander's collection of Henry's law constants
#' @details Where available, gas-specific polynomial fits following the model of Weiss 1970 (equation 4) should be used in preference to the Johnson/Sander scheme. This is partcularly important for long-lived (i.e. close-to equilibrium) gases such as N2, O2, Ar, CO2 and N2O. For shorter lived / far from equilibrium gases the simpler scheme will generally to a good job.
#' @references Johnson, M. T. A numerical scheme to calculate temperature and salinity dependent air-water transfer velocities for any gas. Ocean Sci. Discuss. 7, 251-290 (2010).
#' @references Sander, R. Compilation of Henry’s Law constants for inorganic and organic species of potential importance in environmental chemistry (Version 3), http://www.henrys-law.org, 1999.
#' @references Weiss, R, The solubility of nitrogen, oxygen and argon in water and seawater, Deep Sea Research (17), 721-735, 1970
#' @return partitioning coefficient in moles/litre/atmosphere
#' @export
#'
#' 
GasSolubility<-function(compound,T,S){
  
  Johnson_Sander_GasSolubility<-function(compound,T,S){
    KH_Molar_per_atmosphere <- function(compound,T,S){
      # calculate the Henry's law constant in M/atm from Sander data
      # applying the salting out factor from Johnson 2010
      (compounds[compound,"KH"]*exp((compounds[compound,"tVar"])*((1/(T+273.15))-(1/298.15))))/K_H_factor(compound,S)
    }
  
  }
  # TODO add marelac data method
}
  
  



###############################################################################################
### Under the hood of the Johnson method ######################################################
##############################################################################################

## Salting out (or in for very small compounds) factor from empirical parameterisation by Johnson 2010
Ks<-function(compound){
  theta = (7.3353282561828962e-04 + (3.3961477466551352e-05*log(KH0(compound))) + (-2.4088830102075734E-06*(log(KH0(compound)))^2) + (1.5711393120941302E-07*(log(KH0(compound)))^3))
  theta*log(Vb(compound))
}

K_H_factor <-function(compound,S){
  #calculate salinity-dependent salting-out scaling factor for KH0 (see manuscript)
  10^(Ks(compound)*S)
}

