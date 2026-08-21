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
 
