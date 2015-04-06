# functions for oxygen based NCP calculations


#' O2 NCP
#'
#' calculate NCP based on O2 observations
#'
#' @details TODO
#' @param dat data frame matching the format outlined in XXX
#' @param entrainment if True (default) calculate NCP with entrainment, Cb0 and Cb1 must be supplied
#' @param kw_method character string passed to kw, default is 'WA09'
#' @return a vector of NCP in uMol/l per m-2 per supplied time interval
#' @export
O2NCP <- function(dat, entrainment = T, kw_method = 'WA09'){
    #subfunctions
    Csat.t <- function(tx){
        # calculate oxygen saturation at tx
        return(Csat(temp(tx), sal(tx)))
    }
    rtx <- function(tx){
        # returns r at tx (any point between t0 and t1)
        kw('O2', temp(tx), ws(tx), sal(tx), kw_method) / mld(tx) + dhdt
    }
    Rt <- function(tx){
        # vectorised, integrates rtx up to suppled time
        return(sapply(tx, function(y) integrate(rtx, lower = 0, upper = y)$value))
    }
    Q2.f <- function(tx){
        # exponential of integrated r at supplied time
        return(exp(Rt(tx)))
    }
    Gtx <- function(tx, kw.method){
            # kw at tx
        kwx = kw('O2', temp(tx), ws(tx), sal(tx), method = kw_method)
            # return flux at tx
        Prs = Pslp(tx) / 1000  # surface pressure scaling
        Gx = (kwx / mld(tx)) * Csat.t(tx) * (1 + bubbleSat('O2', ws(tx))) * Prs + (dhdt * Cb(tx))
        return(Gx)
    }
    Q1.f <- function(tx){
        Gtx(tx) * exp(Rt(tx))
    }
    calc_dhdt <- function(h0, h1, times, entrainment){
        # calculate rate of change in mixed layer depth over time
        # returns 0 for negative change
        dhdt = (h1-h0)/max(times)
        if(entrainment & dhdt > 0){
            print('Entraining')
            return(dhdt)
        }
        else{
            return(0)
        }
    }

    ncp = vector()
    for(i in 1:nrow(dat)){
    # variable interpolators, can't be vectorised
    ti = dat[i,]$timePeriod
    dhdt = calc_dhdt(dat[i,]$h0, dat[i,]$h1, ti, entrainment)
        mld  = approxfun(c(0, ti), c(dat[i,]$h0, dat[i,]$h1))
        temp = approxfun(c(0, ti), c(dat[i,]$T0, dat[i,]$T1))
        sal  = approxfun(c(0, ti), c(dat[i,]$S0, dat[i,]$S1))
        ws   = approxfun(c(0, ti), c(dat[i,]$u0, dat[i,]$u1))
        Pslp   = approxfun(c(0, ti), c(dat[i,]$Pslp0, dat[i,]$Pslp1))
        Cb   = approxfun(c(0, ti), c(dat[i,]$Cb0, dat[i,]$Cb1))

    Q2 = integrate(Q2.f, lower = 0, upper = ti)$value
    Q1 = integrate(Q1.f, lower = 0, upper = ti)$value

        # calculate NCP (J)

    J = (dat[i,]$C1 * exp(Rt(ti)) - dat[i,]$C0 - Q1) / Q2
    ncp = c(ncp, J * ti) # return as uMol per supplied time per unit area
    }
    return(ncp)
}

#' O2 NCP simplified (mean)
#'
#' calculate NCP based on O2 observations, using mean variables
#'
#' @details TODO
#' @param dat data frame matching the format outlined in XXX
#' @param entrainment if True (default) calculate NCP with entrainment
#' @param kw_error
#' @param B_error
#' @param Csat_error
#' @return a vector of NCP in uMol/l per m-2 per supplied time interval
#' @export
O2NCP.simple <- function(dat, entrainment = T, kw_error = 0, B_error = 0, Csat_error = 0){
    # expects single row of LHS style data.frame
    # works with all factors constant
    with(dat, {
             if('Cb0' %in% names(dat)){Cb = (Cb0 + Cb1)/2}else{Cb = 0} # check if bottom o2 available
             # calculate averages
             k = (kw('O2', T0, u0, S0) + kw('O2', T1, u1, S1))/2
             k = k + ((k / 100) * kw_error) # apply kw error
             h = (h0 + h1)/2
             Prs = ((Pslp0 + Pslp1)/2) / 1000  # surface pressure scaling
             B = (bubbleSat('O2', u0) + bubbleSat('O2', u1))/2
             B = B + ((B / 100) * B_error) # apply bubble error
             S = (Csat(T0, S0) + Csat(T1, S1))/2
             S = S + ((S / 100) * Csat_error) # apply bubble error
             ti = timePeriod
             if(entrainment == T){
                 dhdt = (h1 - h0)/ti # calculate entrainment
                 dhdt[dhdt < 0] = 0
             }else{
                 dhdt = 0
             }

             q. = (k/h)*S*(1 + B) * Prs + (dhdt * Cb) # q = everything except J that doesn't multiply C
             r = (k / h) + dhdt # -r = everything that multiples C (residence time)

             J = r*(((C1 - q./r) * exp(r * ti) - (C0 - q./r))/ (exp(r * ti) - 1))
             return(J * ti)
           })
}


#' NCP test data
#'
#' example data for testing NCP calculations
#'
#' @format a data frame with 53940 rows and 10 variables:
#'\describe{
#' \item{price}{usdollars}
#' \item{price}{usdollars}
#' }
#' @source FIXME Tom
"ncpTest"
