########################################################
#       GAS PHASE TRANSFER VELOCITY CALCULATIONS       #
########################################################




####################################################################
####      all ka paramterisations return values in m/s          ####
####################################################################



Duce_ka_with_Sc<-function(compound,u,C_d_form,T=15){
    # use C_d = 1 > Duce; C_d = 2 > Smith; C_d = 3 > Large Pond
    if (C_d_form==1) u_star<-Duce_u_star(u)
    if (C_d_form==2) u_star<-Smith_u_star(u)
    if (C_d_form==3) u_star<-Large_Pond_u_star(u)
    ra_turb <-u/(u_star)^2
    ra_diff <- (5/u_star)*((Sc_air(compound,T))^(2/3))
    1/(ra_turb + ra_diff)
    
}



ka <- function(compound,u,T){
   #use new formulation, based on Jeffrey et al. (2010), by default
   Jeffrey_modified_ka(compound,u,T)
   
   #new_ka(compound,u,T)
   #Jeffrey_ka(compound,u,T)
   #Shahin_ka(compound,u,T)
   #MackayYeun_ka(compound,u,T)
   #Duce_ka_with_sc(compound,u,T)
   #Duce_ka(compound,u)
   #Liss_ka(u)
}

#####################################################################
#  RESISTANCE TO TRANSFER, RG/RL and TOTAL TRANSFER VECLOCITIES     #
#####################################################################

rl <- function(compound,T,u,S,schmidt_formulation=0){
    #calculate rl (=1/kw) after Liss and Slater 1974
        1/kw(compound,T,u,S,0,schmidt_formulation)
}

rl_prime <- function(compound,T,u,S,schmidt_formulation=0){
    #rl_prime = KH/kw
        KH(compound,T,S)/kw(compound,T,u,S,0,schmidt_formulation)
}

rg <- function(compound,T,u,S){
    #calculate rg (=1/KH*ka) after Liss and Slater
        1/(KH(compound,T,S)*ka(compound,u,T))
}

rg_prime <- function(compound,u){
    #rg_prime = 1/ka 
        1/ka(compound,u,T)
}

rg_by_rl <- function(compound,T,u,S){
    #calculate rg/rl to give relative contributions of the two phases to the total transfer velocity
        rg(compound,T,u,S)/rl(compound,T,u,S)
}

rgp_by_rlp <- function(compound,T,u,S){
    #calculate rg_prime/rl_prime to give relative contributions of the two phases to the total transfer velocity.
    #if this and rg_by_rl aren't equal for a given gas and set of conditions, something's gone very wrong
        rg_prime(compound,u)/rl_prime(compound,T,u,S)
}

Kw <- function(compound,T,u,S){
    #calculate total transfer velocity (Kw) after Liss and Slater
        1/(rg(compound,T,u,S)+rl(compound,T,u,S))
}

Ka <- function(compound,T,u,S){
    #calculate total transfer velocity (Ka) after Liss and Slater
        1/(rg_prime(compound,u)+rl_prime(compound,T,u,S))
}

#Sanity checking
Sanity <- function(){
cat("Sanity checking - make sure the numbers coming out are approximately what would be expected")
cat("They won't be exactly the same as the de-facto values as calculation methods, Henry's law and other input values etc will not be identical")
cat("\n\n")
cat("first of all, what about the Schmidt numbers for CO2 at 20 Celcius in freshwater and seawater? \n - the commonly quoted values are 600 and 660 respectively, but these are very sensitive so within +/- 30 of these values is fine") 
cat(paste("\nfreshwater: ", schmidt("CO2",20,0)))
cat(paste("\nseawater: ", schmidt("CO2",20,35)))
cat("\n\n")
cat("Hopefully that was OK. How about the values of some transfer velocities?")
cat("\n")
cat("for each formulation, the expected values at 4,10, and 16 m/s and schmidt number of 600 are given, and then the values calculated by this model \n")
cat("1: Liss and Merlivat (1986) \n")
cat("             4m/s    10m/s    16m/s \n")
cat("expected:    2       19       45.5  cm/hr \n")
cat(paste("calculated: ", round(360000*Lissmerkw("CO2",20,4,0,normalize=600),1),"   ", round(360000*Lissmerkw("CO2",20,10,0,normalize=600),1), "   ", round(360000*Lissmerkw("CO2",20,16,0,normalize=600),1)))

cat("\n \n 2: Wanninkhof(1992) \n")
cat("             4m/s    10m/s    16m/s \n")
cat("expected:    5.2     32.5     83   cm/hr \n")
cat(paste("calculated: ", round(360000*Wannkw("CO2",20,4,0,normalize=600),1),"   ", round(360000*Wannkw("CO2",20,10,0,normalize=600),1), "   ", round(360000*Wannkw("CO2",20,16,0,normalize=600),1)))

cat("\n \n 3: Nightingale et al (1999) \n")
cat("             4m/s    10m/s    16m/s \n")
cat("expected:    4.9     25.5     62   cm/hr \n")
cat(paste("calculated: ", round(360000*Nightingkw("CO2",20,4,0,normalize=600),1),"   ", round(360000*Nightingkw("CO2",20,10,0,normalize=600),1), "   ", round(360000*Nightingkw("CO2",20,16,0,normalize=600),1),"\n"))



}


#THE END

