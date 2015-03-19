#           Schmidt number calculations

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
    data(compounds, envir=environment())
    ringval<-ifelse(compounds[compound,"rings"]>0,-7,0)
    ifelse(compounds[compound,"Vb"]>0,compounds[compound,"Vb"],7*(compounds[compound,"C"]+compounds[compound,"H"]+compounds[compound,"O"]+compounds[compound,"N"]+compounds[compound,"db"])+14*compounds[compound,"tb"]+31.5*compounds[compound,"Br"]+24.5*compounds[compound,"Cl"]+10.5*compounds[compound,"F"]+38.5*compounds[compound,"I"]+21*compounds[compound,"S"]+32*compounds[compound,"Se"]+ringval)	
}


#' Schmidt number
#'
#' Schmidt number
#'
#' @details TODO
#' @param compound character string
#' @param T vector of temperature in degrees Centigrade
#' @param S vector of salinity
#' @param method string matching Schmidt caculation method, default is 'mean'
#' @return vector of molar volume at boiling point
#' @keywords volume
#' @references TODO
#' @export
Sch <- function(compound, T, S, method = 'mean'){
    diff_HL <- function(compound,T,S){
        #calculate diffusivity by Hayduk and Laudie (1974) method
        #NOTE - only T dependence from n_sw calculation
        13.26e-5 / ((n_sw(T, S) ^ 1.4) * (Vb(compound) ^ 0.589))
    }
    schmidt_HL <- function(compound, T, S){
        #calculate schmidt number from HL diffusivity
            (v_sw(T, S)) / diff_HL(compound, T, S)
    }

    diff_HM <- function(compound, T, S){
        #Hayduk and Minhas (1982) diffusion coefficient calculation
            EpsilonStar <- (9.58 / Vb(compound))-1.12
            1.25e-8 * (Vb(compound)^(-0.19) - 0.292)*((T + 273.15)^1.52)*((n_sw(T, S)) ^ EpsilonStar)
    }
    schmidt_HM <- function(compound,T,S){
        #calculate schmidt number from HM diffusivity
        (v_sw(T, S)) / diff_HM(compound, T, S)
    }

    diff_WC <- function(compound, T, S){
        #Wilkie and Chang (1955) diffusion coefficient
            # association factor of solvent (2.6 in the case of water according to Poling 2001; although Wanninkhof suggests 2.26)
            phi <- 2.6
            ((T + 273.15) * 7.4e-8 * (phi * 18.01)^ 0.5) / ((n_sw(T, S)) * (Vb(compound)^ 0.6))
    }
    schmidt_WC <- function(compound, T, S){
        #calculate schmidt number from WC diffusivity
            (v_sw(T, S)) / diff_WC(compound, T, S)
    }

    mean.schmidt <- function(compound, T, S){
        #calculate mean schmidt number in water
            mean_diff <- 0.5 * (diff_WC(compound, T, S) + diff_HM(compound, T, S))
            v_sw(T, S) / mean_diff
    }

    switch(method,
           HL = schmidt_HL(compound, T, S),
           HM = schmidt_HM(compound, T, S),
           WC = schmidt_WC(compound, T, S),
           mean = mean.schmidt(compound, T, S)
           )
}
