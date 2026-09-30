import numpy as np
from numba import njit 

#############################################
#
# Plant root water uptake functions
#
#############################################
@njit(inline='always')
def rootStressFunction(psie,T,pars):
    # We use loops here, because they are highly readable, and numba will compile them efficiently
    beta=np.zeros(len(psie))
    for i in range(len(psie)):
        if T[i]<0:
            # No transpiration from frozen soils
            beta[i]=0.
        elif psie[i]<pars['psi_wilt']:
            # No transpiration from below wilting point
            beta[i]=0.
        elif psie[i]<pars['psi_opt']:
            # Limited transpiration below optimum psi:
            beta[i]=(psie[i]-pars['psi_wilt'])/(pars['psi_opt']-pars['psi_wilt'])
        elif psie[i]<pars['psi_crit']:
            # Unstressed transpiration below critical psi:
            beta[i]=1.
        else: 
            # No transpiration above critical psi
            beta[i]=0.

    return beta

@njit(inline='always')
def soilEvapStressFunction(psie,T,pars):
    # Only consider the top cell here, i=0:

    if T[0]<0:
        # No evaporation if the soil is frozen
        gamma=0.
    elif psie[0]<pars['psi_soilE_min']:
        # No evaporation if the soil is dry
        gamma=0.
    elif psie[0]<pars['psi_soilE_max']:
        # Limited evaporation below soilE_max
        gamma=(psie[0]-pars['psi_soilE_min'])/(pars['psi_soilE_max']-pars['psi_soilE_min'])
    else:
        # Unstressed evaporation
        gamma=1.

    return gamma


@njit(inline='always')
def rootDensityFunction(z,dz,pars):
    # Get relative root density, gr(z)
    g=np.exp(-z/pars['rootDepth'])
    gr=g/np.sum(g*dz)
    return gr

@njit(inline='always')
def rootUptake(E_PT,psie,T,z,dz,pars):
    beta=rootStressFunction(psie,T,pars)
    gr=rootDensityFunction(z,dz,pars)
    sv=-E_PT*beta*gr
    return sv


@njit(inline='always')
def soilEvaporation(E_PS,psie,T,dz,pars):
    gamma=soilEvapStressFunction(psie,T,pars)
    sv=-E_PS*gamma/dz[0]
    return sv
 
@njit(inline='always')
def getEvapFluxes(E_PT,E_PS,psie,T,z,dz,nt,pars,const):
    E_AT=np.zeros(nt)
    E_AS=np.zeros(nt)
    jE=np.zeros(nt)
    
    for i in range(nt):
        sv=rootUptake(E_PT[i],psie[i,:],T[i,:],z,dz,pars)
        E_AT[i]=-np.sum(sv*dz)
        s_ES=soilEvaporation(E_PS[i],psie[i,:],T[i,:],dz,pars)
        E_AS[i]=-s_ES*dz[0]
        sv[0]+=s_ES
        su=sv*const['cp_liq']*const['rho_liq']*T[i,:]
        jE[i]=-np.sum(su*dz)

    E_AT[1:]=(E_AT[1:]+E_AT[:-1])/2.
    E_AT[0]=0.
    E_AS[1:]=(E_AS[1:]+E_AS[:-1])/2.
    E_AS[0]=0.
    jE[1:]=(jE[1:]+jE[:-1])/2.
    jE[0]=0.

    return E_AT,E_AS,jE




# @njit(inline='always')
def satvappresFun(T): # saturation vapour pressure (Pa)
    es = 0.6108*np.exp((17.27*T)/(T+237.3))*1000
    return es

# @njit(inline='always')
def DeltaFun(T,es):  # Delta (Pa/deg C)
    Delta = (4098*es)/(T+237.3)**2
    return Delta

# @njit(inline='always')
def gammaFun(p,const):  # gamma (Pa/deg c)
    gamma=const['cp_air']*p/const['epsilon']/const['lambda_v']
    return gamma

# @njit(inline='always')
def vappresFun(RH,es): # actual vapour pressure (Pa)
    ea = (RH/100)*es
    return ea

# @njit(inline='always')
def windFun(uz,pars): # calculate wind speed at 2m height (m/s)
    u2 = uz*4.87/np.log(67.8*pars['z_u']-5.42) # uz wind speed m/s measured at hieght z (m)
    return u2

# The Penman Monteith Combination Equation
# @njit(inline='always')
def GetPotentialEvap(SWnet,LWnet,Ta,P,RH,U,G,pars,const,rs):
    
    # Input variables with units:
    Rn=(SWnet+LWnet)       # W/m2
    # P                    # Pa
    es=satvappresFun(Ta)   # Pa
    ea=vappresFun(RH,es)   # Pa
    gamma =gammaFun(P,const)     # Pa/deg C
    Delta= DeltaFun(Ta,es) # Pa/deg C
            
    # Aerodynamic resistance:
    d  = (2/3)*pars['canopyHeight']   # zero displacement height
    z0m = 0.123*pars['canopyHeight']  # roughness length of momentum
    z0h = 0.1*z0m                     # roughness length for heat and vapour transfer 
    k  = 0.41                         # von Karman constant

    ra = (
        np.log((pars['z_u']-d)/z0m)
        *np.log((pars['z_rh']-d)/(z0h))
    )/(k**2)/U 

    # Density of air:
    
    rho_air=P/const['Rd']/(Ta+273.15)
    
    # PM-Combination Eqn to return latent heat flux:
    numerator=(
        Delta*(Rn-G) + 
        const['cp_air']*rho_air*(es-ea)/ra
    )
    
    denominator=(
        Delta+gamma*(1+rs/ra)
    )

    latentHeatFlux=numerator/denominator   # J/m2/s

    # Get evaporation rate in mm/d
    PE=latentHeatFlux/const['lambda_v']    # kg/m2/s
    PE=PE*86400                            # mm/d

    # Remove negative values
    PE[PE<0]=0
    
    return PE
