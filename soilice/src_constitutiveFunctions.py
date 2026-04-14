import numpy as np
from numba import njit 

#############################################
#
# USER SPECIFIED CONSTITUITIVE FUNCTIONS
# Maybe editted. 
# Do not change function inputs/outputs
#
#############################################

# Hydraulic properties:
@njit(inline='always')
def thetaFun(psi,pars):
    Se=(1+(psi*-pars['alpha'])**pars['n'])**(-pars['m'])
    Se[psi>0.]=1.0
    return pars['thetaR']+(pars['thetaS']-pars['thetaR'])*Se

@njit(inline='always')
def CFun(psi,pars):
    Se=(1+(psi*-pars['alpha'])**pars['n'])**(-pars['m'])
    Se[psi>0.]=1.0
    dSedh=pars['alpha']*pars['m']/(1-pars['m'])*Se**(1/pars['m'])*(1-Se**(1/pars['m']))**pars['m']
    return Se*pars['Ss']+(pars['thetaS']-pars['thetaR'])*dSedh

@njit(inline='always')
def KFun(psie,psif,pars,const):
    # VG K fun, using psie (call function with psie=psif if you want to overide that)
    Se=(1+(psie*-pars['alpha'])**pars['n'])**(-pars['m'])
    Se[psie>0.]=1.0
    K=pars['Ks']*Se**pars['neta']*(1-(1-Se**(1/pars['m']))**pars['m'])**2

    # Optinally include impedance (set pars['impedance']=0. to remove this)
    # NOTE - do not include impedance if you are using cryoK=1 - it will not work
    thetaL=thetaFun(psif,pars)
    thetaT=thetaFun(psie,pars)
    thetaI=const['rho_liq']/const['rho_ice']*(thetaT-thetaL)
    K=K*10**(-pars['impedance']*thetaI)
     
    return K

# Thermal conductivity function
@njit(inline='always')
def thermalKfun(psie,psif,T,pars,const):
    # Uses geometric mean of each component
    thetaL=thetaFun(psif,pars)
    thetaT=thetaFun(psie,pars)
    thetaI=const['rho_liq']/const['rho_ice']*(thetaT-thetaL)
    theta_mineral=pars['theta_mineral']
    theta_org=pars['theta_org']
    thetaA=pars['thetaS']-thetaT
    kappa = (
        (const['kappa_liq'] ** thetaL) *
        (const['kappa_ice'] ** thetaI) *
        (const['kappa_air'] ** thetaA) *
        (pars['kappa_soil'] ** theta_mineral)*
        (pars['kappa_org'] ** theta_org))
    return kappa

# Bulk heat capacity function
@njit(inline='always')
def CBFun(psie,psif,pars,const):
    # Uses arithmetic mean of each component
    thetaL=thetaFun(psif,pars)
    thetaT=thetaFun(psie,pars)
    thetaI=const['rho_liq']/const['rho_ice']*(thetaT-thetaL)
    theta_mineral=pars['theta_mineral']
    theta_org=pars['theta_org']
    #thetaA=thetaS-thetaT ~ assumed negligible
    CB = (
        (const['cp_ice']*thetaI*const['rho_ice']) +
        (const['cp_liq']*thetaL*const['rho_liq']) +
        (pars['cp_soil']*theta_mineral*pars['rho_soil']) +
        (pars['cp_org']*theta_org*pars['rho_org']))
    return CB

# Slope function dtheta_L/dT
@njit(inline='always')
def SFCslope(psie,psif,pars,const):
    C=const['lambda_f']/const['g']/const['T0']
    dthdpsi=CFun(psif,pars)
    dthdT=C*dthdpsi
    return dthdT
    
# GCE
@njit(inline='always')
def GCEFun(T,pars,const):
    psi=T*const['lambda_f']/(const['g']*const['T0'])
    psi[psi>0]=0.
    return psi


