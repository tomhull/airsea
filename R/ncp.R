#' Oxygen concentration saturation
#'
#' Calculates oxygen saturation concentration in equilibrium with the atmosphere
#' as per Garcia & Gordon, 1992 (Benson & Kraus, rather than combined)
#' @param temp numeric vector of water temperature in degrees Celsius
#' @param salinity numeric vector of salinity (PSU)
#' @return vector of saturation concentration in micro mol per litre
#' @keywords oxygen
#' @examples
#' Csat(10, 35)  # saturation at 10 degrees and 35 salinity
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

#' NCP model
#'
#' Solves for predicted oxygen for use in NCP estimation
#' Called by run_ncp functions.
#' @param t a integer vector timing sequence
#' @param y the initial oxygen concentration i.e. O2 at t0
#' @param parms a list of parameter values
#' @return vector of predicted oxygen concentrations at solver timesteps
#' @keywords oxygen NCP
ncpModel <- function(t, y, parms){
    with(as.list(c(y, parms)),{
         Csat <- Csat(temperature(t), salinity(t))
         o = C / Cstar
         Prs = Pslp(t) / 1000  # surface pressure scaling
         GasExchange <- (kw(t) / mixedLayer(t)) * Csat * ((1 + B(t)) * Prs - o) # positive gas exchange increases the concentration in the mixed layer
         Entrainment <- dhdt * (Cbottom(t) - C) # positive entrainment increases the concentration in the mixed layer
         dC <- GasExchange + Entrainment
        list(c(dC))
    })
}

#' Rate of change in mixed layer
#'
#' Calculates rate of change in mixed layer depth for NCP estimation
#' returns positive value for MLD deepening i.e. Entrainment, otherwise returns 0
#' called by run_ncp functions
#' @param h0 the mixed layer depth at time 0
#' @param h1 the mixed layer depth at time 1
#' @param time a vector containing the maximum time interval
#' @return the rate of change in mixed layer depth, or 0 if h1 < h0
#' @keywords mixed layer depth entrainment
calc_dhdt <- function(h0, h1, time, entrainmentOff){
    # calculate rate of change in mixed layer depth over time
    # returns 0 for negative change
    dhdt = (h1-h0)/max(time)
    if((dhdt < 0) | entrainmentOff == T){
        return(0)
    }
    else{
        return(dhdt)
    }
}
