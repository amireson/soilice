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

@njit(inline='always')
def gamma(p):  # gamma (KPa/deg c)
    g = (0.665/1000)*p # p Atmo. pressure (KPa)
    return g

@njit(inline='always')
def Delta(T,es):  # Delta (KPa/deg c)
    d = (4098*es)/(T+237.3)**2
    return d

@njit(inline='always')
def satvappres(T): # saturation vapour pressure (KPa)
    es = 0.6108*np.exp((17.27*T)/(T+237.3))
    return es

@njit(inline='always')
def vappres(RH,es): # actual vapour pressure (KPa)
    ea = (RH/100)*es
    return ea

@njit(inline='always')
def wind(uz,z): # wind speed at 2m height (m/s)
    u2 = uz*4.87/np.log(67.8*z-5.42) # uz wind speed m/s measured at hieght z (m)
    return u2


@njit(inline='always')
def GetPotentialEvap(SWnet,LWnet,T,p,RH,uz,z,zavg,rs): # input data: SW&LW(W/m2), T(deg c), p(Pa), RH(%), uz(m/s), z(m)
    
    #Equation parameters with the righ units:    
    
    Rn=(SWnet+LWnet)*(0.0864)   # MJ/m2/d (as 30*60 is the number of secondes in period)
    p = p/1000         # KPa
    es=satvappres(T)   # KPa
    ea=vappres(RH,es)  # KPa
    g =gamma(p)        # KPa/deg c
    delta= Delta(T,es)    # KPa/deg c
    
    # u2 (m/s)
    if z==2:
        u2=uz
    else:
        u2=wind(uz,z)
        
    # calculations for the aerodynamic resistance:
    d  = (2/3)*zavg   #zero displacment height
    z0 = 0.123*zavg # roughness length of momentum
    k  = 0.41       # von Karman constant
    
    #1
    ra = (np.log((2-d)/z0)*np.log((2-d)/(z0*0.1)))/(k**2)/u2  # the FAO eqaution (s/m)

    #2
#     Uf = (k*u2)/np.log((2-d)/z0) # Dingman book equation
#     ra=u2/(Uf**2)
    

    # calculate the equation terms:
    a = 0.408*delta*Rn
    b = (g*185396.2*(es-ea)/((T+273)))*(1/ra)
    c = delta+g*(1+(rs/ra))     #*(1+0.34*u2)
    
    PE = (a+b)/c
    
    PE[PE<0]=0
    return PE
