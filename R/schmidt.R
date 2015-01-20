# schmidt number calculations

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
