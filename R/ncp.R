# functions for oxygen based NCP calculations


#' timeseries translator
#'
#' divides up timeseries for NCP calculation
#'
#' @details TODO
#' @param x a data frame matching the format outlined in XXX
#' @return data.frame in the format expected by O2NCP()
#' @export
transform_data <- function(x){
    # transforms timeseries data to inital and post variables in wide format, i.e. like a lhs
    dat = data.frame(dateTime = x$dateTime, T0 = x$T, S0 = x$S, C0 = x$C, Cb0 = x$Cb, u0 = x$u, h0 = x$h, Pslp0 = x$Pslp)
    dat$T1 = dat$T0 + c(diff(dat$T0), 0)
    dat$S1 = dat$S0 + c(diff(dat$S0), 0)
    dat$Cb1 = dat$Cb0 + c(diff(dat$Cb0), 0)
    dat$u1 = dat$u0 + c(diff(dat$u0), 0)
    dat$C1 = dat$C0 + c(diff(dat$C0), 0)
    dat$B0 = bubbleSat(dat$u0)
    dat$B1 = bubbleSat(dat$u1)
    dat$kw0 = kw('O2', dat$T0, dat$u0, dat$S0)
    dat$kw1 = kw('O2', dat$T1, dat$u1, dat$S1)
    dat$h1 = dat$h0 + c(diff(dat$h0), 0)
    if("hMin" %in% colnames(x)){
        dat$hMin0 = x$hMin
        dat$hMin1 = x$hMin + c(diff(x$hMin), 0)
        dat$hMax0 = x$hMax
        dat$hMax1 = x$hMax + c(diff(x$hMax), 0)
    }
    dat$Pslp1 = dat$Pslp0 + c(diff(dat$Pslp0), 0)
    dat$timePeriod = c(diff(as.numeric(dat$dateTime)), 0)
    return(dat[-nrow(dat),])
}

#' O2 NCP
#'
#' calculate NCP based on O2 observations
#'
#' @details TODO
#' @param dat data frame matching the format outlined in XXX
#' @param entrainment if True (default) calculate NCP with entrainment, Cb0 and Cb1 must be supplied
#' @param kw.method character string passed to kw, default is 'WA09'
#' @return a vector of NCP in uMol/l
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

    ncp = vector()
    for(i in 1:nrow(dat)){
    # variable interpolators, can't be vectorised
    ti = dat[i,]$timePeriod
    dhdt = 0
        mld  = approxfun(c(0, ti), c(dat[i,]$h0, dat[i,]$h1))
        temp = approxfun(c(0, ti), c(dat[i,]$T0, dat[i,]$T1))
        sal  = approxfun(c(0, ti), c(dat[i,]$S0, dat[i,]$S1))
        ws   = approxfun(c(0, ti), c(dat[i,]$u0, dat[i,]$u1))
        Cb   = approxfun(c(0, ti), c(dat[i,]$Cb0, dat[i,]$Cb1))

    Q2 = integrate(Q2.f, lower = 0, upper = ti)$value
    Q1 = integrate(Q1.f, lower = 0, upper = ti)$value

        # calculate NCP (J)

    J = (dat[i,]$C1 * exp(Rt(ti)) - dat[i,]$C0 - Q1) / Q2
    ncp = c(ncp, J * ti) # return as uMol per supplied time
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
