#       Viscosity and density of water

#' LaLiberte viscosity
#'
#' Viscosity of pure water according to LaLiberte 2007
#'
#' @details TODO
#' @param T vector of temperature in degrees Centigrade
#' @return viscosity of pure water in cP (mPa s-1)
#' @keywords viscosity
#' @references LaLiberte 2007
n_0 <- function(T){
        (T+246)/(137.37+(5.2842*T)+(0.05594*(T^2)))
}

n_sw <- function(T,S,return_relative_viscotiy=0){
    #calculate the dynamic viscosity of seawater in cP using the viscosity model / mixing rule approach of LaLiberte 2007
    #to return relative viscosity (i.e. viscosity scaling factor) for comparison with other studies, pass n_sw a final argument of "1"
    #after T and S e.g. n_sw(10,35,1)
        #read in a data file containing the mass fraction of each constituent solute in seawater per salinity unit in per mil and the coefficients determined by LaLiberte 2007 for each salt
        # sw_cmf<-read.table("sw_constituent_mass_fractions.dat", header = TRUE, row.names = 1)
        #w_i_ln_n_i_tot is the sum of the product of the mass fraction and the natural log of the viscosities contributed by each solute individually (see LaLiberte 2007, equation 8)
        w_i_ln_n_i_tot <-0
        #sum up the mass fractions to get water mass fraction
        w_i_tot<-0
        for (salt in row.names(sw_cmf)){
            w_i <- sw_cmf[salt,"massfraction"]*S/1000
            w_i_tot<-w_i_tot + w_i
        }
        for (salt in row.names(sw_cmf)){
            w_i <- sw_cmf[salt,"massfraction"]*S/1000
            v1 <- sw_cmf[salt,"v1"]
            v2 <- sw_cmf[salt,"v2"]
            v3 <- sw_cmf[salt,"v3"]
            v4 <- sw_cmf[salt,"v4"]
            v5 <- sw_cmf[salt,"v5"]
            v6 <- sw_cmf[salt,"v6"]
            #wi_tot is used here as eq 12 in LaLiberte requires (1-w_w), which is equivalent. Using this term seems initially counterintuitive - one might expect to use w_i here for each solute individually. However,  "at any given solute concentration, the solute viscosity increases with the total solute content" so 1-w_w (or w_i_tot) is the correct term - pers. comm. from Marc LaLiberte
            n_i<-(exp(((v1*w_i_tot^v2)+v3)/((v4*T) + 1)))/((v5*(w_i_tot^v6))+1)
            w_i_ln_n_i_tot <- w_i_ln_n_i_tot + (w_i*log(n_i))
        }
        ln_n_m<- (1-w_i_tot)*log(n_0(T))+w_i_ln_n_i_tot
        if (return_relative_viscotiy==0) exp(ln_n_m) else
            (exp(ln_n_m)/n_0(T))
}

n_sw_hardy <- function(T,n_zero=1){
    #Non-salinity dependent Hardy  1953 method for calculating viscosity - not used in calculations, purely for comparison with Laliberte method.
    # η = η0 * K/[1+0.03338.T + 0.00018325.T2]
    # use n_zero = 1 to use Laliberte n_0 or n_zero=0 for Hardy n_0
        if (n_zero==0) n0<- 1.787 else
            n0 <- n_0(T)
        (n0*1.052/(1+(0.03338*T)+(0.00018325*(T^2))))
}

#' Seawater density
#'
#' density of seawater according to Millero and Poisson (1981)
#'
#' @details TODO
#' @param T vector of temperature in degrees Centigrade
#' @return density of seawater in kg m-3
#' @keywords seawater density
#' @references Millero and Poisson (1981)
#' @export
p_sw <- function(T,S){
        #coefficients for millero calc
        A <- 0.824493-(0.0040899*T)+(0.000076438*(T^2))-(0.00000082467*(T^3))+(0.0000000053875*(T^4))
        B <- -0.00572466+(0.00010277*T)-(0.0000016546*(T^2))
        C <- 0.00048314
        # density of pure water
        p_o <- 999.842594+(0.06793952*T)-(0.00909529*(T^2))+(0.0001001685*(T^3))-(0.000001120083*(T^4))+(0.000000006536332*(T^5))
        # return salinity-and-temperature-dependent density of seawater (at atmospheric pressure)
        (p_o+(A*S)+(B*(S^(3/2)))+(C*S))
}

v_sw <- function(T,S) {
    #calculate kinmatic viscosity of seawater in cm/s for Schmidt number calculation
        # dynamic viscosity - convert to S.I. (Kg/m/s)
        n = n_sw(T,S)/1000
        # density already in S.I.
        p = p_sw(T,S)
        #multiply by 10000 to go from m2/s to cm2/s
        10000*n/p
}
