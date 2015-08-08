# implementation of liang 2013 piston and bubble stuff
# using liang matlab names where possible

solub <- function(T, S, igas){
  
  # Returns the solubility [mol/(m^3 Pa)] as a function of temperature
  # and salinity for the following gases as per Weiss rather than Sarmiento & Gruber, 2006:
  # igas
  # 1 = O2
  # 2 = N2
  # 3 = CO2
  # 4 = Ar
  # 5 = Ne
  # 6 = He
  # 7 = Kr
  # 8 = Xe
  # 9 = DMS
  
  t = (T + 273.15) * 0.01 # t = kelvin/100
  
  # Oxygen
  if(igas == 1){
    soly = -173.4292 + 249.6339 / t + 143.3483 * log(t) - 21.8492 * t + S * (-0.033096 + 0.014259 * t -0.0017 * t^2)
    return( exp(soly) / (22.4*101325*0.20946) )
  }
  # CO2
  if(igas==3){
    soly = -160.7333 + 215.4152 / t + 89.8920 * log(t) - 1.47759 * t^2 + S * (0.029941 -0.027455 * t + 0.0053407 *t^2)
    return( exp(soly) * 1e3 / 101325 )
  }
  
  # TODO remaining gasses
}

solubility <- function(x, gasname){
  # cites wannikhoff 1992
  # x = water temp in oC
  # sol = solubility in mole/m-3/atm
  # al = dimensionless solubility
  # sc schmidt number
  r1 = 8.314
  tsx=(x + 273.16) / 100 # I don't know why he changes is nomeclature here, same as t
  
  if(gasname == 'CO2'){
    a = c(-60.2409, 93.4157, 23.3585)
    b = c(0.023517, -0.023656, 0.0047036)
    a0 = c(2073.1, 125.62, 3.6276, 0.043219)
    s = a[1] + a[2] / tsx + a[3] * log(tsx) + 35 * ( b[1] + b[2] * tsx + b[3] * tsx * tsx)
    sol = exp(s) * 1000
  }
  if(gasname == 'O2'){
    a = c(-58.3877, 85.8079, 23.8439)
    b = c(-0.034892,  0.016241, -0.0020114)
    a0 = c(1953.4, 128.0, 3.9918, 0.050091);
    sol = solub(x, 35, 1) * 101325 * 0.20946
  }
  al = sol / 1E5 * r1 * tsx * 100
  sc = a0[1] - a0[2] * x + a0[3] * x^2 - a0[4] * x^3
  return(c(sol, al, sc))
}

piston_velocity_L12 <- function(u10, sst, gasname){
  
  y  = solubility(sst, gasname) # assumes 35 salinity!
  
  alc = y[2]
  Sc = y[3]
  
  cd10 = dragCoef(u10)

  usr = sqrt(cd10) * u10
             
  lam = 13.3
  A = 1.3
  phi = 1
  rhow = 1022
  rhoa = 1.0
  hw = lam / A / phi
  scwc = Sc
  tkt = 0.01

  rwo = sqrt(rhow / rhoa) * (hw * sqrt(scwc) + (log(0.5 / tkt) / 0.4))
   
  ha = lam
  scac = 0.9  # %air-side schmidt number
   
  ra = ha * sqrt(scac) + 1 / sqrt(cd10) - 5 + 0.5 * log(scac) / 0.4 # %air side
     
  usrw = usr / sqrt(1000)
  #   
  kbb = 1.98e6 * usrw^(2.76) / ((scwc / 660)^(2/3)) / 100 / 3600
  #   
  rw = rwo
  #   
  kw = 1 / (1 / (usr /rw / sqrt(scwc / 660) + kbb) + 1 / (usr / (ra * alc)))
  return(kw)
}

gasex_overpressure_L12 <- function(u10, sst, gasname){
  
  k0 = piston_velocity_L12(u10, sst, gasname)
  
  cd10 = dragCoef(u10)
  
  usr = sqrt(cd10) * u10 / sqrt(1000)
  
  y  = solubility(sst, gasname)
  sol=y[1]
  scwc=y[3]
  
  sigmap = 1.5244 * usr^(1.06)
  finj = 5.56 * usr^(3.86)
  
  kbb = 1.98e6 * usr^(2.76) / ((scwc / 660) ^ (2 / 3)) / 100 / 3600
  
  del = (kbb * sigmap * y[1] + finj) / k0 / y[1]
  
  if(gasname == 'O2'){chi = 0.20946}
  if(gasname == 'N2'){chi=0.78084}
  if(gasname == 'CO2'){chi=0.00036}
  if(gasname == 'Ar'){chi=0.00934}
  
  del = kbb * (sigmap + finj * chi / kbb / y[1]) / k0
  return(del)
}

# comp = data.frame(u10 = seq(0.5, 20, length.out = 100))
# comp$liang_10 = gasex_overpressure_L12(comp$u10, 10, 'O2')
# comp$liang_15 = gasex_overpressure_L12(comp$u10, 15, 'O2')
# comp$liang_20 = gasex_overpressure_L12(comp$u10, 20, 'O2')
# comp$woolf= bubbleSat('O2', comp$u10)
# ggplot( melt(comp, id.var = 'u10') ) + geom_line(aes(u10, value, color = variable)