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
#' @param x data frame matching the format outlined in XXX
#' @param entrainment if True (default) calculate NCP with entrainment, Cb0 and Cb1 must be supplied
#' @return a vector of NCP in uMol/l
#' @export
O2NCP <- function(x, entrainment = T){
    with(x, {

    ti = timePeriod

     if(entrainment == T){
        if('Cb0' %in% names(x) & 'Cb1' %in% names(x)){
            stop('Cb0 or Cb1 not supplied')
        }else{
            # calculate dhdt
            if(h1 - h0 > 0){
                dhdt = (h1 - h0) / ti
            }

        }
    }else{
        dhdt = 0
     }

    r0 = (kw0 / h0) + dhdt # -r = everything that multiples C (residence time + entrainment rate)
    r1 = (kw1 / h1) + dhdt # -r = everything that multiples C

    rtx = approxfun(c(0, ti),c(r0, r1))

    Rt <- function(x){
        # vectorised
        return(sapply(x, function(y) integrate(rtx, lower = 0, upper = y)$value))
    }

    Q2.f <- function(x){
        return(exp(Rt(x)))
    }
    Q2 = integrate(Q2.f, lower = 0, upper = ti)$value

        # q

    qtx <- function(x){
        kw = approxfun(c(0, ti),c(kw0, kw1))
        h = approxfun(c(0, ti),c(h0, h1))
        B = approxfun(c(0, ti),c(B0, B1))
        Cb = approxfun(c(0, ti),c(Cb0, Cb1))
        temp = approxfun(c(0, ti), c(T0, T1))
        sal = approxfun(c(0, ti), c(S0, S1))

        CstarGG.t <- function(x){
            return(CstarGG(sal(x), temp(x)))
        }
        qx = (kw(x) / h(x)) * CstarGG.t(x) * (1 + B(x)) * 1 + (dhdt * Cb(x))
        return(qx)
    }

    Q1.f <- function(x){
        qtx(x) * exp(Rt(x))
    }
    Q1 = integrate(Q1.f, lower = 0, upper = ti)$value
    print(c(Q1, Q2, C0, C1, exp(Rt(ti)), qtx(ti)))

        # J

    J = (C1 * exp(Rt(ti)) - C0 - Q1) / Q2
    return(J * ti)
           })
}
