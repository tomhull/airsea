# implementation of liang 2013 piston and bubble stuff
# using liang matlab names where possible
# require(airsea)

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
  if(igas == "O2"){
      # soly coefs taken from table.2 Weiss (volumetric solubilty constants)
    soly = -173.4292 + 249.6339 / t + 143.3483 * log(t) - 21.8492 * t + S * (-0.033096 + 0.014259 * t -0.0017 * t^2)
    return( exp(soly) / (22.4*101325*0.20946) ) # 22.4 * mean sea level pressure (Pa) * air mole fraction
  }
  # Ar
  if(igas == "Ar"){
      # soly coefs taken from table.2 Weiss (volumetric solubilty constants)
    soly = -173.5146 + 245.4510 / t + 141.8222 * log(t) - 21.8020 * t + S * (-0.034474 + 0.014934 * t -0.0017729 * t^2);
    return( exp(soly) / (22.4*101325*0.00934) )
  }
  # CO2
  if(igas == "CO2"){
    soly = -160.7333 + 215.4152 / t + 89.8920 * log(t) - 1.47759 * t^2 + S * (0.029941 -0.027455 * t + 0.0053407 *t^2)
    return( exp(soly) * 1e3 / 101325 )
  }
  
  # is this equiv to martin KH when converted to SI
}

solubility <- function(x, gasname, S = 35){
  # cites wannikhoff 1992
  # x = water temp in oC
  # sol = solubility in mole/m-3/atm
  # al = dimensionless solubility
  # sc schmidt number
  r1 = 8.314 # gas constant
  tsx=(x + 273.16) / 100 # I don't know why he changes is nomeclature here, same as t
  ts = x + 273.16
  # print(paste("solubility called", gasname, x))
  
  if(gasname == 'CO2'){
    a = c(-60.2409, 93.4157, 23.3585)
    b = c(0.023517, -0.023656, 0.0047036)
    a0 = c(2073.1, 125.62, 3.6276, 0.043219)
    s = a[1] + a[2] / tsx + a[3] * log(tsx) + S * ( b[1] + b[2] * tsx + b[3] * tsx * tsx)
    sol = exp(s) * 1000
  }
  if(gasname == 'O2'){
    a = c(-58.3877, 85.8079, 23.8439) # these are the bunsen solubility constants from weiss
    b = c(-0.034892,  0.016241, -0.0020114) # as above
    a0 = c(1953.4, 128.0, 3.9918, 0.050091)
    sol = solub(x, S, "O2") * 101325 * 0.20946 # convert to atm and * mole fraction
  }
  if(gasname == 'Ar'){
    a = c(-55.6578, 82.0262, 22.5929) # these are the bunsen solubility constants from weiss
    b = c(-0.036267,  0.078374, -0.0033120) # as above
    a0 = c(1909.1, 125.09, 3.9012, 0.048953)
    sol = solub(x, S, "Ar") * 101325 * 0.00934 # convert to atm and then?
  }
  
  al = sol / 1E5 * r1 * ts # convert Hcp to ?dimensionless
  sc = a0[1] - a0[2] * x + a0[3] * x^2 - a0[4] * x^3
  return(data.frame(sol, al, sc))
}

#' Liang et al 2012 Kw param
#'
#' @param u10 
#' @param sst 
#' @param sal 
#' @param gasname 
#'
#' @return vector of kw in meters per second
#' @export
#'
kw_L12 <- function(u10, sst, sal = 35, gasname){
  
  # print(paste("piston called", gasname, u10, sst, sal))
  y  = solubility(sst, gasname, sal) # assumes 35 salinity!
  
  alc = y$al # dimentionless henry's law solubility
  Sc = y$sc # water side schmidtt
  # Sc = Sch(gasname, sst, sal)
  
  cd10 = dragCoef(u10) # drag coeficient

  usr = sqrt(cd10) * u10
             
  lam = 13.3 # 
  A = 1.3 # fairall2011 empirical constant
  phi = 1 # surface bouyancey flux enhancement compensation
  rhow = 1022 # water density
  rhoa = 1.0 # atmospheric density
  hw = lam / A / phi # ?incorrect, fairall defines this as hw = 13.3 / (A * phi)
  kappa = 0.4 # from Fairall2011 eq 3
  tkt = 0.01
  scac = 0.93  # air-side schmidt number, why this isn't calculated i don't know
  # Scac = Sc_air(gasname, sst)  # air-side schmidt number

    # water side resistance to transfer (due to molecular turbulent processes)
    # from fairall2011 eq13 sort of
  rw = sqrt(rhow / rhoa) * (hw * sqrt(Sc) + (log(0.5 / tkt) / kappa))
   
   
    # this is sort of fairall2011 eq 15
  ra = lam * sqrt(scac) + 1 / sqrt(cd10) - 5 + 0.5 * log(scac) / kappa # air side resistance
     
  usrw = usr / sqrt(1000) # frictional velocty ? doesn't use the definition in the paper
  #   
  kbb = 1.98e6 * usrw^(2.76) / ((Sc / 660)^(2/3)) / 100 / 3600
  #   
  kw = 1 / (1 / (usr /rw / sqrt(Sc / 660) + kbb) + 1 / (usr / (ra * alc)))
  return(kw) # in 
}

#' Liang et al 2012 Bubble param
#'
#' @param u10 
#' @param sst 
#' @param gasname 
#'
#' @return bubble sat term (dimensionless)
#' @export
B_L12 <- function(u10, sst, gasname, sal){
  
  k0 = kw_L12(u10, sst, sal, gasname = gasname)
  
  cd10 = dragCoef(u10)
  
  usr = sqrt(cd10) * u10 / sqrt(1000) # 
  
  y  = solubility(sst, gasname, sal)
  sol=y$sol # henrys law volumetric solubility from weiss
  scwc=y$sc # schmidt number
  
  sigmap = 1.5244 * usr^(1.06)
  finj = 5.56 * usr^(3.86)
  
  kbb = 1.98e6 * usr^(2.76) / ((scwc / 660) ^ (2 / 3)) / 100 / 3600
  
  del = (kbb * sigmap * y[1] + finj) / k0 / y[1]
  
  if(gasname == 'O2'){chi=0.20946}
  if(gasname == 'N2'){chi=0.78084}
  if(gasname == 'CO2'){chi=0.00036}
  if(gasname == 'Ar'){chi=0.00934}
  
  del = kbb * (sigmap + finj * chi / kbb / sol) / k0
  return(del)
}

# plot(piston_velocity_L12(0:20, 20, 35, "CO2")*100*60*60)
# points(kw("CO2", 20, 0:20, 35)*100*60*60, col = "red")

