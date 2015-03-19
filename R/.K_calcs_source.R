########################################################
#       GAS PHASE TRANSFER VELOCITY CALCULATIONS       #
########################################################


Tucker_D_air <- function(compound,T){
    #calculate diffusivity in air in cm2/sec
    #M_a is molar weight of air
    M_a <- 28.97
    M_b <- compounds[compound,"mw"]
    M_r <- (M_a + M_b)/(M_a*M_b)
    #assume 1ATM
    P <- 1
    #assume molar volume air is 20.1 cm3/mol
    V_a <- 20.1	
    (0.001*((T+273.15)^1.75)*sqrt(M_r))/(P*((V_a^(1/3))+(Vb(compound)^(1/3))))^2
}

D_air <- function(compound,T){
    Tucker_D_air(compound,T)
} 

n_air <- function(T){
    # dynamic viscosity of saturated air according to Tsiligiris 2008
    SV_0 = 1.715747771e-5
    SV_1 = 4.722402075e-8
    SV_2 = -3.663027156e-10
    SV_3 = 1.873236686e-12
    SV_4 = -8.050218737e-14
    
    # in N.s/m^2 (Pa.s)
    u_m = SV_0+(SV_1*T)+(SV_2*T^2)+(SV_3*T^3)+(SV_4*T^4)
    u_m
}

p_air <- function(T){
    # density of saturated air according to Tsiligiris 2008 in kg/m^3
    SD_0 = 1.293393662
    SD_1 = -5.538444326e-3
    SD_2 = 3.860201577e-5
    SD_3 = -5.2536065e-7
    p = SD_0+(SD_1*T)+(SD_2*T^2)+(SD_3*T^3)
    p
}

v_air <- function(T) {
    #calculate kinmatic viscosity of air in cm2/s for Schmidt number calculation
        # dynamic viscosity 
        n = n_air(T)
        # density 
        p = p_air(T)
        #multiply by 10000 to go from m2/s to cm2/s
        10000*n/p
}

Sc_air <- function (compound,T){
    #calculate the schmidt number of a given gas in air
    v_air(T)/D_air(compound,T)
}


####################################################################
####      all ka paramterisations return values in m/s          ####
####################################################################

Duce_ka <- function(compound,u){
    #calculate ka gas phase transfer velocity according to Duce et al 1991
        u/(770+(45*((compounds[compound,"mw"])^(1/3))))
}

Duce_ka_with_Sc<-function(compound,u,C_d_form,T=15){
    # use C_d = 1 > Duce; C_d = 2 > Smith; C_d = 3 > Large Pond
    if (C_d_form==1) u_star<-Duce_u_star(u)
    if (C_d_form==2) u_star<-Smith_u_star(u)
    if (C_d_form==3) u_star<-Large_Pond_u_star(u)
    ra_turb <-u/(u_star)^2
    ra_diff <- (5/u_star)*((Sc_air(compound,T))^(2/3))
    1/(ra_turb + ra_diff)
    
}


Liss_ka <- function(u){
    #calculate ka according to Liss 1973 (non compound-specific)
    (0.005+(0.21*u))/100
}

MackayYeun_ka <- function(compound,u,T)
    #ka according to Mackay and Yeun 1983
    1e-3 + (46.2e-5*sqrt(6.1+(0.63*u))*u*(Sc_air(compound,T))^-0.67)

Shahin_ka <- function(compound,u,T){
    #calculate transfer velocity from gas phase diffusivity in air according to Shahin et al 2002 in m/s
    #checked against shahin data - they find a ka of ~3 cm/s at 6m/s wind - approx 5 times that of Duce!)- see Johnson 2010 (final article) for an explanation
    (sqrt(D_air(compound,T))*((0.98*u)+1.26))/100
}

Jeffrey_ka<-function(compound,u,T){
    #using smith(1980) Cd
    von_kar<-0.4
    Cd<-(1e-4*(6.1+0.63*u))
    Sc<-Sc_air(compound,T)
    ra<-13.3*sqrt(Sc) + (Cd^(-0.5)) - 5 + log(Sc)/(2*von_kar)
    u_star<-u*sqrt(C_D(u))
    u_star/ra
}

new_ka<-function(compound,u,T){
#this is the ka propsed in Johnson 2010 (OS discussions paper), which has been superceded by a modified version of Jeffrey 2010 (scheme_ka below) for the resubmitted  paper which will hopefully be published in OS
 u_star<-u*sqrt(C_D(u))
 ra_turb <-u/(u_star)^2
 ra_diff <-(5/u_star)*((Sc_air(compound,T))^(2/3))
 katest<-1e-3 +1/(ra_turb + ra_diff)
 ifelse(is.na(katest),1e-3,katest)
}

Jeffrey_modified_ka<-function(compound,u,T){
    (1e-3+Jeffrey_ka(compound,u,T))
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

