# functions for oxygen based NCP calculations


#' O2 NCP
#'
#' calculate NCP based on O2 observations, assumes linear change in observed values between
#' timesteps
#'
#' @details TODO
#' @param dat data frame matching the format outlined in XXX
#' @param entrainment if True (default) calculate NCP with entrainment, Cb0 and Cb1 must be supplied
#' @param kw_method character string passed to kw, default is 'WA13'
#' @param asVolume if True NCP is not multipled by average mixed layer depth
#' @return a vector of NCP in mmol per m-3 per supplied time interval, or mmol per m-3 if asVolume is True.
#' @export
O2NCP.linear <- function(dat, entrainment = T, kw_method = 'WA13', asVolume = F){
    #subfunctions
    Csat.t <- function(tx){
        # calculate oxygen saturation at tx
             S = Csat(temp(tx), sal(tx))
             S = S + ((S / 100) * Csat_error) # apply bubble error
        return(S)
    }
    rtx <- function(tx){
        # returns r at tx (any point between t0 and t1)
        k = kw('O2', temp(tx), ws(tx), sal(tx), kw_method)
        k = k + ((k / 100) * kw_error) # apply kw error
        k / mld(tx) + dhdt
    }
    Rt <- function(tx){
        # vectorised, integrates rtx up to suppled time
        return(sapply(tx, function(y) integrate(rtx, lower = 0, upper = y)$value))
    }
    Q2.f <- function(tx){
        # exponential of integrated r at supplied time
        return(exp(Rt(tx)))
    }
    Gtx <- function(tx){
            # kw at tx
        kwx = kw('O2', temp(tx), ws(tx), sal(tx), method = kw_method)
        kwx = kwx + ((kwx / 100) * kw_error) # apply kw error
            # return flux at tx
        Prs = Pslp(tx) / 1000  # surface pressure scaling
        B = bubbleSat('O2', ws(tx))
        B = B + ((B / 100) * B_error) # apply bubble error
        Gx = (kwx / mld(tx)) * Csat.t(tx) * (1 + B ) * Prs + (dhdt * Cb(tx))
        return(Gx)
    }
    Q1.f <- function(tx){
        return(Gtx(tx) * exp(Rt(tx)))
    }
    calc_dhdt <- function(h0, h1, times, entrainment){
        # calculate rate of change in mixed layer depth over time
        # returns 0 for negative change
        dhdt = (h1-h0)/max(times)
        if(entrainment & dhdt > 0){
            return(dhdt)
        }
        else{
            return(0)
        }
    }

    ncp = vector()
    for(i in 1:nrow(dat)){
    # variable interpolators, can't be vectorised
    ti = dat$timePeriod[i]
    dhdt = calc_dhdt(dat$h0[i], dat$h1[i], ti, entrainment)
        mld  = approxfun(c(0, ti), c(dat$h0[i], dat$h1[i]))
        temp = approxfun(c(0, ti), c(dat$T0[i], dat$T1[i]))
        sal  = approxfun(c(0, ti), c(dat$S0[i], dat$S1[i]))
        ws   = approxfun(c(0, ti), c(dat$u0[i], dat$u1[i]))
        Pslp   = approxfun(c(0, ti), c(dat$Pslp0[i], dat$Pslp1[i]))
        if("Cb0" %in% colnames(dat)){
            Cb   = approxfun(c(0, ti), c(dat$Cb0[i], dat$Cb1[i]))
        }else{
            Cb   = approxfun(c(0, ti), c(0, 0))
        }
            # error vectors
        if("kw_error" %in% colnames(dat)){
            kw_error = dat$kw_error[i]
        }else{ kw_error = 0 }
        if("B_error" %in% colnames(dat)){
            B_error = dat$B_error[i]
        }else{ B_error = 0 }
        if("Csat_error" %in% colnames(dat)){
            Csat_error = dat$Csat_error[i]
        }else{ Csat_error = 0 }

    Q2 = integrate(Q2.f, lower = 0, upper = ti)$value
    Q1 = integrate(Q1.f, lower = 0, upper = ti)$value

        # calculate NCP (J)

    J = (dat$C1[i] * exp(Rt(ti)) - dat$C0[i] - Q1) / Q2
    if(asVolume == T){
        ncp = c(ncp, J * ti) # return as mmol m-3 per supplied time
        }
    else{
        ncp = c(ncp, J * ti * ((dat$h0[i] + dat$h1[i])/2)) # return as mmol m-2 per supplied time
        }
    }
    return(ncp)
}

#' O2 NCP simplified (mean)
#'
#' calculate NCP based on O2 observations, using mean values of observations
#'
#' @details TODO
#' @param dat data frame matching the format outlined in XXX
#' @param asVolume if True NCP is not multipled by average mixed layer depth
#' @param kw_method character string passed to kw, default is 'WA13'
#' @return a vector of NCP in mmol per m-3 per supplied time interval, or mmol per m-3 if asVolume is True.
#' @export
O2NCP.mean <- function(dat, asVolume = F, kw_method = 'WA13'){
    # expects single row of LHS style data.frame
    # works with all factors constant
    if(!"kw_error" %in% colnames(dat)){kw_error = 0}
    if(!"B_error" %in% colnames(dat)){B_error = 0}
    if(!"Csat_error" %in% colnames(dat)){Csat_error = 0}
        # if entrainment state variables found use them
    if("entrainment" %in% colnames(dat)){use_entrainment = T}else{use_entrainment = F}

    with(dat, {

             # calculate averages
             k = (kw('O2', T0, u0, S0) + kw('O2', T1, u1, S1))/2
             k = k + ((k / 100) * kw_error) # apply kw error
             h = (h0 + h1)/2
             Prs = ((Pslp0 + Pslp1)/2) / 1000  # surface pressure scaling
             B = (bubbleSat('O2', u0) + bubbleSat('O2', u1))/2
             B = B + ((B / 100) * B_error) # apply bubble error
             S = (Csat(T0, S0) + Csat(T1, S1))/2
             S = S + ((S / 100) * Csat_error) # apply bubble error
             if('Cb0' %in% names(dat)){
                 Cb = (Cb0 + Cb1)/2
             }
             else{Cb = 0} # check if bottom o2 available
             ti = timePeriod

            # are we going to calculate entrainment?
             if(use_entrainment == T){
                dhdt = (h1 - h0)/ti # calculate entrainment
                dhdt[dhdt < 0] = 0
                dhdt[entrainment == F] = 0
             }else{
                 # if not set to 0 for no entrainment
                 dhdt = 0
             }

             r = (k / h) + ((1 / h) * dhdt) #  = everything that multiples C (residence time)
             f. = (k/h)*S*(1 + B) * Prs + (dhdt * Cb) # q = everything except J that doesn't multiply C

             J = r * h * ((C1 - C0) / (1 - exp(-r * ti)) + C0) - f. * h
             if(asVolume == T){
                 return((J / h) * ti) # return as mmol m-3 per supplied time
             }
             else{
                 return(J * ti) # return as mmol m-2 per supplied time
             }
           })
}
