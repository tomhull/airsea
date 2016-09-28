# functions for oxygen based NCP calculations

#' Net community production (NCP) from O2
#'
#' calculate net community production based on mean O2 observations.
#'
#' @details Full explanation of the method and methods are found in Hull et al,
#'   2015 Entrainment will only be calculated where an
#'   \code{entrainment} variable is found in the input data and it's value is
#'   \code{TRUE}. When being using within a monte-carlo analysis you can include
#'   bias estimates in the following variables \code{kw_error, B_error,
#'   Csat_error}. Where these variables are found these parameters will be
#'   scaled accordingly.
#' @param dat data frame matching the format outlined in \link[airsea]{O2NCP.transform}
#' @param kw_method character string passed to kw, default is 'NG00'
#' @param bubbles logical value to enable bubble supersaturation param, default is TRUE (on)
#' @param entrainment logical value to enable entrainment calculation, default is FALSE (off)
#' @param output_conc logical value to enable expressing NCP as just a biological concentration change rather than a flux.
#' @param output_conc logical value to print debug messages
#'  (units of mmol m-3), default is FALSE (off)
#' @return a vector of NCP in mmol per m-3 per supplied time interval
#' @references Hull et al, 2015 http://www.biogeosciences-discuss.net/12/15611/2015/bgd-12-15611-2015.html
#' @export
O2NCP.mean <- function(dat, kw_method = 'NG00', bubbles = T, entrainment = F, output_conc = F, debug = F){
  # expects single row of LHS style data.frame
  # works with all factors constant
  if(!"kw_error" %in% colnames(dat)){kw_error = 0}
  if(!"B_error" %in% colnames(dat)){B_error = 0}
  if(!"Csat_error" %in% colnames(dat)){Csat_error = 0}
      # if entrainment state variables found use them
  if(!("Cb0" %in% colnames(dat)) & entrainment == T){entrainment = F; if(debug){print("Cb not found, no entrainment calculated")}}

    with(dat, {
        ti = timePeriod
        # calculate averages
        k = (kw('O2', T0, u0, S0, method = kw_method) + kw('O2', T1, u1, S1, method = kw_method))/2
        k = k + ((k / 100) * kw_error) # apply kw error
        h = (h0 + h1)/2
        Prs = ((Pslp0 + Pslp1)/2) / 1013.25  # surface pressure scaling
        B = (bubbleSat('O2', u0) + bubbleSat('O2', u1))/2
        B = B + ((B / 100) * B_error) # apply bubble error
        S = (Csat(T0, S0) + Csat(T1, S1))/2
        S = S + ((S / 100) * Csat_error) # apply bubble error

        # are we going to calculate entrainment?
        if(entrainment == T){
            Cb = (Cb0 + Cb1)/2
            dhdt = (h1 - h0)/ti # calculate entrainment
            dhdt[dhdt < 0] = 0
            if(debug){print(paste("calculating entrainment, dhdt=", dhdt, "Cb=", Cb))}
            dhdt[is.na(Cb)] = 0
            Cb[is.na(Cb)] = 0
        }else{
            # if not set to 0 for no entrainment
            dhdt = 0
            Cb = 0
        }

        if(bubbles == F){
          if(debug){print("bubbles off")}
          B = 0
          Prs = 1
        }
        r = (k / h) + ((1 / h) * dhdt) #  = everything that multiples C
        if(debug){print(paste("r =", r))}
        Q = (k/h)*S*(1 + B) * Prs + ((1/h) * dhdt * Cb) # Q = everything except J that doesn't multiply C
        if(debug){print(paste("Q =", Q))}
        J = r * h * ((C1 - C0) / (1 - exp(-r * ti)) + C0) - Q * h
        if(debug){print(paste("J =", J))}
        if(output_conc == T){
          J = J / h
        }
        return(J * ti) # return as mmol m-2 per supplied time
        })
}

#' transform data to O2NCP format
#'
#' @details TODO
#' expects column headings to match the following format (order does not matter):
#' dateTime = POSIXct, T = temperature (oC), S = salinity, C = oxygen concentration (mmol m-3), u = wind speed (m s-1),
#' Pslp = pressure at sea level (mbar), h = mixed layer depth (m).
#' Optionally Cb = bottom oxygen (mmol m-3) if entrainment is to be calculated.
#' Optionally hmin for minimum mld for uncertainty analysis
#' 
#' @param x data frame of observations
#'
#' @return data frame in correct format for O2NCP
#' @export
O2NCP.transform <- function(x){
    # transforms timeseries data to inital and post variables in wide format, i.e. like a lhs
    dat = with(x,
               data.frame(dateTime = dateTime,
                          T0 = T, S0 = S,
                          C0 = C, u0 = u,
                          h0 = h, Pslp0 = Pslp))
    dat$T1 = dat$T0 + c(diff(dat$T0), 0)
    dat$S1 = dat$S0 + c(diff(dat$S0), 0)
    dat$u1 = dat$u0 + c(diff(dat$u0), 0)
    dat$C1 = dat$C0 + c(diff(dat$C0), 0)
    if("Cb" %in% colnames(x)){
      dat$Cb0 = x$Cb
      dat$Cb1 = dat$Cb0 + c(diff(dat$Cb0), 0)
    }
    if("hmin" %in% colnames(x)){
      dat$hmin0 = x$hmin
      dat$hmin1 = dat$hmin0 + c(diff(dat$hmin0), 0)
    }
    dat$h1 = dat$h0 + c(diff(dat$h0), 0)
    dat$Pslp1 = dat$Pslp0 + c(diff(dat$Pslp0), 0)
    dat$timePeriod = c(diff(as.numeric(dat$dateTime)), 0)
    dat$end = dat$timePeriod + dat$dateTime
    return(dat[-nrow(dat),])
}
