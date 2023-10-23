#' molar volume
#'
#' Calculates the molar volume at boiling point using the Schroeder method, or takes overruling value from compounds data.
#' note: for compounds containing elements other than C,H,O,N,S, Cl, Br, F, and I, the Schroeder method won't work so Vb must be specified
#' db, tb and rings should be the number of double and triple bonds, and ring features in the molecule respectively.
#'
#' @details TODO
#' @param compound character string
#' @return vector of molar volume at boiling point
#' @keywords volume
#' @references TODO
#' @export
Vb <- function(compound){
    ringval<-ifelse(compounds[compound,"rings"]>0,-7,0)
    ifelse(compounds[compound,"Vb"]>0,compounds[compound,"Vb"],7*(compounds[compound,"C"]+compounds[compound,"H"]+compounds[compound,"O"]+compounds[compound,"N"]+compounds[compound,"db"])+14*compounds[compound,"tb"]+31.5*compounds[compound,"Br"]+24.5*compounds[compound,"Cl"]+10.5*compounds[compound,"F"]+38.5*compounds[compound,"I"]+21*compounds[compound,"S"]+32*compounds[compound,"Se"]+ringval)	
}

#' Schmidt number
#'
#' Schmidt number
#'
#' @details Calculates Schmidt number using one of several methods:
#' \itemize{
#'  \item JS - Mean of WC and HM as per Johnson, 2000, the default
#'  \item WA - Wanningkov, 2014 (assumes salinity = 35)
#'  \item HL - Hayduk and Laudie, 1974
#'  \item HM - Hayduk and Minhas, 1982
#'  \item WC - Wilkie and Chang, 1955
#' }
#' @param compound character string
#' @param TEMP vector of temperature in degrees Celsius
#' @param SAL vector of salinity
#' @param method string matching Schmidt calculation method, default is 'JS'
#' @return vector of molar volume at boiling point
#' @keywords volume
#' @references TODO
#' @export
Sch <- function(compound, TEMP, SAL, method = 'JS'){
    diff_HM <- function(compound, TEMP, SAL){
        #Hayduk and Minhas (1982) diffusion coefficient calculation
            EpsilonStar <- (9.58 / Vb(compound))-1.12
            1.25e-8 * (Vb(compound)^(-0.19) - 0.292)*((TEMP + 273.15)^1.52)*((n_sw(TEMP, SAL)) ^ EpsilonStar)
    }
    schmidt_HM <- function(compound, TEMP, SAL){
        #calculate schmidt number from HM diffusivity
        (v_sw(TEMP, SAL)) / diff_HM(compound, TEMP, SAL)
    }
    diff_HL <- function(compound, TEMP, SAL){
        #calculate diffusivity by Hayduk and Laudie (1974) method
        #NOTE - only T dependence from n_sw calculation
        13.26e-5 / ((n_sw(TEMP, SAL)^1.4)*(Vb(compound)^0.589))
    }
    schmidt_HL <- function(compound, TEMP, SAL){
        #calculate schmidt number from HL diffusivity
            (v_sw(TEMP, SAL)) / diff_HL(compound, TEMP, SAL)
    }

    diff_WC <- function(compound, TEMP, SAL){
        #Wilkie and Chang (1955) diffusion coefficient
            # association factor of solvent (2.6 in the case of water according to Poling 2001; although Wanninkhof suggests 2.26)
            phi <- 2.6
            ((TEMP + 273.15) * 7.4e-8 * (phi * 18.01)^ 0.5) / ((n_sw(TEMP, SAL)) * (Vb(compound)^ 0.6))
    }
    schmidt_WC <- function(compound, TEMP, SAL){
        #calculate schmidt number from WC diffusivity
            (v_sw(TEMP, SAL)) / diff_WC(compound, TEMP, SAL)
    }

    schmidt_JS <- function(compound, TEMP, SAL){
        #calculate mean schmidt number in water
            mean_diff <- 0.5 * (diff_WC(compound, TEMP, SAL) + diff_HM(compound, TEMP, SAL))
            v_sw(TEMP, SAL) / mean_diff
    }
    
    schmidt_WA <- function(compound, TEMP){
        # polynomial fit as per Wanningkov 2014, assumes salinity = 35
        # which is based on Hayduk and Laudie for O2
        coef = wann14_coef[wann14_coef$Gas == compound,]
        coef$A + TEMP * coef$B + TEMP^2 * coef$C + TEMP^3 * coef$D + TEMP^4 * coef$E
    }

    switch(method,
           HL = schmidt_HL(compound, TEMP, SAL),
           HM = schmidt_HM(compound, TEMP, SAL),
           WC = schmidt_WC(compound, TEMP, SAL),
           WA = schmidt_WA(compound, TEMP),
           JS = schmidt_JS(compound, TEMP, SAL)
           )
}
