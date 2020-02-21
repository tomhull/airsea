vp_h2o <- function(temp, sal){
  # Vapor pressure of seawater in Pa (Green & Carritt, 1967)
  tk = temp + 273.15
  (1 - 0.000537 * sal) * 101325 * exp(18.1973 * (1 - 373.16 / tk) + 3.1813E-07 * (1 - exp(26.1205 * (1 - tk / 373.16))) - 0.018726 * (1 - exp(8.03945 * (1 - 373.16 / tk))) + 5.02802 * log(373.16 / tk))
}

vp_h2o_hamme <- function(temp, sal){
  # temp: ºC
  # sal: psu
  # vp_h2o_hamme: Pa (0.01 mbar)
  #   Vapour pressure of pure water: D. Ambrose and I.J. Lawrenson (1972)
  #    "The vapour pressure of water"
  #    Journal of Chemical Thermodynamics, v.4, p. 755-671.
  #   Correction for seawater: Frank J. Millero and Wing H. Leung
  #    "The thermodynamics of seawater at one atmosphere" (1976)
  #    American Journal of Science, v. 276, p. 1035-1077.
  tk = temp + 273.15
  Tmod = (2 * tk - (648 + 273)) / (648 - 273)
  
  # Calculate value of Chebyshev polynomial 
  Chebyshev = (2794.0144 / 2) + (1430.6181 * (Tmod)) + (-18.2465 * (2 * Tmod ^ 2 - 1)) + (7.6875 * (4 * Tmod ^ 3 - 3 * Tmod)) + (-0.0328 * (8 * Tmod ^ 4 - 8 * Tmod ^ 2 + 1)) + (0.2728 * (16 * Tmod ^ 5 - 20 * Tmod ^ 3 + 5 * Tmod)) + (0.1371 * (32 * Tmod ^ 6 - 48 * Tmod ^ 4 + 18 * Tmod ^ 2 - 1)) + (0.0629 * (64 * Tmod ^ 7 - 112 * Tmod ^ 5 + 56 * Tmod ^ 3 - 7 * Tmod)) + (0.0261 * (128 * Tmod ^ 8 - 256 * Tmod ^ 6 + 160 * Tmod ^ 4 - 32 * Tmod ^ 2 + 1)) + (0.02 * (256 * Tmod ^ 9 - 576 * Tmod ^ 7 + 432 * Tmod ^ 5 - 120 * Tmod ^ 3 + 9 * Tmod)) + (0.0117 * (512 * Tmod ^ 10 - 1280 * Tmod ^ 8 + 1120 * Tmod ^ 6 - 400 * Tmod ^ 4 + 50 * Tmod ^ 2 - 1)) + (0.0067 * (1024 * Tmod ^ 11 - 2816 * Tmod ^ 9 + 2816 * Tmod ^ 7 - 1232 * Tmod ^ 5 + 220 * Tmod ^ 3 - 11 * Tmod))
  
  # Vapor pressure of pure water in Pa
  vapor_0sal_Pa = 1000 * 10 ^ (Chebyshev / tk)
  
  # Correct vapour pressure for salinity
  return( vapor_0sal_Pa + 101325 / 760 * (sal * (-0.0023311 - 0.00014799 * temp - 7.52E-06 * temp ^ 2 - 5.5185E-08 * temp ^ 3) + sal ^ 1.5 * (-1.132E-05 - 8.7086E-06 * temp + 7.4936E-07 * temp ^ 2 - 2.6327E-08 * temp ^ 3)))
}
