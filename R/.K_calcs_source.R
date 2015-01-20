
##########################################################
##		DRAG coefficient C_D			##
##########################################################

C_D<-function(u,formulation=1){
    if (formulation==1) {smithCD(u)} else {-999}
}

############################################################
############################################################
##      LIQUID PHASE TRANSFER VELOCITY CALCULATIONS       ##
############################################################
############################################################

#################################################
#all diffusion coefficients calculated in cm^2/s#
#################################################

diff_HM <- function(compound,T,S){
    #Hayduk and Minhas (1982) diffusion coefficient calculation
        EpsilonStar <- (9.58/Vb(compound))-1.12
        1.25e-8*(Vb(compound)^(-0.19)-0.292)*((T+273.15)^1.52)*((n_sw(T,S))^EpsilonStar)
}

schmidt_HM <- function(compound,T,S){
    #calculate schmidt number from HM diffusivity
    (v_sw(T,S))/diff_HM(compound,T,S)
}

diff_HL <- function(compound,T,S){
    #calculate diffudivity by Hayduk and Laudie (1974) method
    #NOTE - only T dependence from n_sw calculation
    13.26e-5/(((n_sw(T,S))^1.4)*(Vb(compound)^0.589))
}

schmidt_HL <- function(compound,T,S){
    #calculate schmidt number from HM diffusivity
        (v_sw(T,S))/diff_HL(compound,T,S)
}

diff_WC <- function(compound,T,S){
    #Wilkie and Chang (1955) diffusion coefficient
        # associaction factor of solvent (2.6 in the case of water according to Poling 2001; although Wanninkhof suggests 2.26)
        phi <- 2.6
        ((T+273.15)*7.4e-8*(phi*18.01)^0.5)/((n_sw(T,S))*(Vb(compound)^0.6))
}

schmidt_WC <- function(compound,T,S){
    #calculate schmidt number from WC diffusivity
        (v_sw(T,S))/diff_WC(compound,T,S)
}

schmidt <- function(compound,T,S){
    #calculate mean schmidt number in water
        mean_diff<-0.5*(diff_WC(compound,T,S)+diff_HM(compound,T,S))
        v_sw(T,S)/mean_diff
}

schmidt_out<-function(compound,T,S,schmidt=0){
    #output schmidt number T dependence for a list of T
    #calculated with mean schmidt number by default, or specify 1 for Hayduck and Minhas or 2 for Wilkie and Chang
        if (schmidt==0) schmidt(compound,T,S) else
            if (schmidt==1) schmidt_HM(compound,T,S) else
                schmidt_WC(compound,T,S)
}



#########################################################
###   liquid phase transfer velcocity calculations    ###
###           all values returned in m/s              ###
#########################################################


Wannkw <- function (compound,T,u,S,normalize=0,schmidt_formulation=0){
    #calculate kw transfer velocity in m/s according to Wanninkhof 1992
    #specifying a 'normalize' value allows the calculation of the transfer velocity for a user-specified value of schmidt number
        if (normalize!=0) schmidt_number<-normalize else
        schmidt_number<-schmidt_out(compound,T,S,schmidt_formulation)
        (0.31*u^2*((schmidt_number/660)^(-0.5)))/(100*3600)
}

Lissmerkw<- function(compound,T,u,S,normalize=0,schmidt_formulation=0){
        k600<-ifelse(u<3.6,0.17*u,
            ifelse(u<13,(2.85*u)-9.65,(5.9*u)-49.3))
            if (normalize!=0) schmidt_number<-normalize else
                schmidt_number<-schmidt_out(compound,T,S,schmidt_formulation)
        ifelse(u<3.6,(k600*((schmidt_number/600)^(-0.66)))/(100*3600),(k600*((schmidt_number/600)^(-0.5)))/(100*3600))
}

Lissmer_approx_kw <- function(compound,T,u,S,normalize=0,schmidt_formulation=0){
    #calculate kw in m/s according to liss and merlivat 
    # note k600 not k660
    # in order to make calculations easier, an exponential relationship has been fitted to the Liss and Merlivat Curve)
        if (normalize!=0) schmidt_number<-normalize else
            schmidt_number<-schmidt_out(compound,T,S,schmidt_formulation)
        (0.075*u^2.25*((schmidt_number/600)^(-0.5))/(100*3600))
}

Sweeneykw <- function (compound,T,u,S,normalize=0,schmidt_formulation=0){
    #calculate kw according to Sweeney 2007
        if (normalize!=0) schmidt_number<-normalize else
            schmidt_number<-schmidt_out(compound,T,S,schmidt_formulation)
        (0.27*u^2*((schmidt_number/660)^(-0.5)))/(100*3600)
}

Nightingkw <- function(compound,T,u,S,normalize=0,schmidt_formulation=0){
    # empirical fit to dual tracer data by Nightingale et al 2000 (GBC)
    # note k600 not k660 for their study
        if (normalize!=0) schmidt_number<-normalize else
            schmidt_number<-schmidt_out(compound,T,S,schmidt_formulation)
        (((0.222*u^2)+0.333*u)*(schmidt_number/600)^(-0.5))/(100*3600)
}

Woolf97kw <- function(compound,T,u,S,normalize=0,wb_formulation=0,schmidt_formulation=0){
    if (normalize!=0) schmidt_number<-normalize else
            schmidt_number<-schmidt_out(compound,T,S,schmidt_formulation)
        (((0.222*u^2)+0.333*u)*(schmidt_number/600)^(-0.5))/(100*3600)
    # whitecapping component included in this physically based model of gas exchange
    # W_b is whitecap coverage
    #different whitecapping paramaterisations could be applied here:
    #0 is woolf default from Monahan..
    if (wb_formulation==0) W_b<-3.84e-6 * (u^(3.41))
    
    # K_b is whitecap-dependent (ostensibly bubble-related) component of total k_w
    # beta is equal to liquid over gas unitless KH using water T for both phases (equal to 1/KH in this scheme).  
        
    beta<-1/KH(compound,T,S)
    denominator<-(beta*(1+((14*beta*(schmidt_number^(-0.5)))^(-1/1.2)))^1.2)	
    K_b = 2540*W_b/denominator
    # K_o is non-whitecapping friction velocity relationship component of k_w
    u_star<-u*sqrt(C_D(u))
    K_o<-1.57e-4*u_star*((600/schmidt_number)^0.5)
    ((K_o*360000) + K_b)/360000
}

McGilliskw<-function(compound,T,u,S,normalize=0,schmidt_formulation=0){
    #empirical fit to GasEx data by McGillis et al (JGR 2001)
    schmidt_number<-schmidt_out(compound,T,S,schmidt_formulation)
    ((3.3+(0.026*(u^3)))*((schmidt_number/600)^(-0.5)))/(100*3600)
}

kw <- function(compound,T,u,S,normalize=0,schmidt_formulation=0) {
    #calculate mean kw. set normalize to schmidt number value e.g. 660 to compare with other kw curves. Doesn't matter what compound or value of T is passed to the function (unless using the Woolf1997 bubble scheme), only u and schmidt value to normalise to, but you have to pass it some values otherwise R will break.
        sf<-schmidt_formulation
        nmlz<-normalize
        Nightingkw(compound,T,u,S,nmlz,sf)
    #Woolf97kw(compound,T,u,S,normalize=nmlz,schmidt_formulation=sf)
    #Wannkw(compoundT,u,S,normalize=nmlz,schmidt_formulation=sf)
}



########################################################
#              HENRY'S LAW COEFFICIENTS                #
########################################################

KH0 <- function(compound,T=25){
    #Calculate Henry's law constant at a given T in pure water according to Sander (1999)
        12.2/((273.15+T)*(compounds[compound,"KH"])*exp((compounds[compound,"tVar"])*((1/(T+273.15))-(1/298.15))))
}


Ks<-function(compound){
    theta = (7.3353282561828962e-04 + (3.3961477466551352e-05*log(KH0(compound))) + (-2.4088830102075734E-06*(log(KH0(compound)))^2) + (1.5711393120941302E-07*(log(KH0(compound)))^3))
    theta*log(Vb(compound))
}

K_H_factor <-function(compound,S){
    #calculate salinity-dependent salting-out scaling factor for KH0 (see manuscript)
        10^(Ks(compound)*S)
}

KH <- function(compound,T,S){
    #Calculate gas-over-liquid unitless Henry's law constant at a given T and Salinity
        KH0(compound,T)*K_H_factor(compound,S)
}

KH_Molar_per_atmosphere <- function(compound,T,S){
    # calculate the Henry's law constant in M/atm from Sander data
    # applying the salting out factor above
        (compounds[compound,"KH"]*exp((compounds[compound,"tVar"])*((1/(T+273.15))-(1/298.15))))/K_H_factor(compound,S)
}





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

