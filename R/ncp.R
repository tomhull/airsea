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
        Gx = (kwx / mld(tx)) * Csat.t(tx) * (1 + bubbleSat('O2', ws(tx))) * 1 + (dhdt * Cb(tx))
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
        Cb   = approxfun(c(0, ti), c(dat[i,]$Cb0, dat[i,]$Cb1))

    Q2 = integrate(Q2.f, lower = 0, upper = ti)$value
    Q1 = integrate(Q1.f, lower = 0, upper = ti)$value

        # calculate NCP (J)

    J = (dat[i,]$C1 * exp(Rt(ti)) - dat[i,]$C0 - Q1) / Q2
    ncp = c(ncp, J * ti) # return as uMol per supplied time per unit area
    }
    return(ncp)
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
