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
