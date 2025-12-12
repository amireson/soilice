import numpy as np
import pandas as pd
import time
from scipy.integrate import ode
from numba import jit 

from numba import types
from numba.typed import Dict

def MakeDictArray():
    d=Dict.empty(
    key_type=types.unicode_type,
    value_type=types.float64[:],)
    return d

def MakeDictFloat():
    d=Dict.empty(
    key_type=types.unicode_type,
    value_type=types.float64,)
    return d

# Thermal conductivity function
@jit(nopython=True)
def thermalKfun(psie,psif,T,pars,const):
    # Uses geometric mean of each component
    thetaL=thetaFun(psif,pars)
    thetaT=thetaFun(psie,pars)
    thetaI=const['rho_liq']/const['rho_ice']*(thetaT-thetaL)
    thetaG=pars['thetaS']-thetaT
    kappa = (
        (const['kappa_liq'] ** thetaL) *
        (const['kappa_ice'] ** thetaI) *
        (const['kappa_air'] ** thetaG) *
        (pars['kappa_soil'] ** pars['theta_mineral'])*
        (pars['kappa_org'] ** pars['theta_org']))
    return kappa
    
# GCE
@jit(nopython=True)
def GCEfun(T,pars,const):
    psi=T*const['lambda_f']/(const['g']*const['T0'])
    psi[psi>0]=0.
    return psi

# Slope function dtheta_L/dT
@jit(nopython=True)
def gdashfun(T,pars,const):
    psi=GCEfun(T,pars,const)
    
    C=const['lambda_f']/const['g']/const['T0']
    U=1+np.abs(pars['alpha']*C*T)**pars['n']
    dthdSe=(pars['thetaS']-pars['thetaR'])
    dSedu=-pars['m']*U**(-pars['m']-1)  
    dudT=pars['n']*(pars['alpha']*C)**pars['n']*(np.abs(T))**(pars['n']-1)
    dthdT=-dthdSe*dSedu*dudT
    dthdT[T>0]=0. 
    return dthdT

# Bulk heat capacity function
@jit(nopython=True)
def CBfun(psie,psif,pars,const):
    thetaL=thetaFun(psif,pars)
    thetaT=thetaFun(psie,pars)
    thetaI=const['rho_liq']/const['rho_ice']*(thetaT-thetaL)
    thetaS=1-pars['thetaS']
    thetaG=thetaS-thetaT
    CB=(const['cp_ice']*thetaI*const['rho_ice'])+(const['cp_liq']*thetaL*const['rho_liq'])+(pars['cp_soil']*thetaS*pars['rho_soil'])
    return CB

# Hydraulic properties:
@jit(nopython=True)
def thetaFun(psi,pars):
    Se=(1+(psi*-pars['alpha'])**pars['n'])**(-pars['m'])
    Se[psi>0.]=1.0
    return pars['thetaR']+(pars['thetaS']-pars['thetaR'])*Se

@jit(nopython=True)
def fdashFun(psi,pars):
    Se=(1+(psi*-pars['alpha'])**pars['n'])**(-pars['m'])
    Se[psi>0.]=1.0
    dSedh=pars['alpha']*pars['m']/(1-pars['m'])*Se**(1/pars['m'])*(1-Se**(1/pars['m']))**pars['m']
    return Se*pars['Ss']+(pars['thetaS']-pars['thetaR'])*dSedh

@jit(nopython=True)
def KFun(psie,psif,pars,const):
    # Impedance model for K after Taylor and Luthin
    thetaL=thetaFun(psif,pars)
    thetaT=thetaFun(psie,pars)
    thetaI=const['rho_liq']/const['rho_ice']*(thetaT-thetaL)

    Se=(1+(psie*-pars['alpha'])**pars['n'])**(-pars['m'])
    Se[psie>0.]=1.0
    Ke=pars['Ks']*Se**pars['neta']*(1-(1-Se**(1/pars['m']))**pars['m'])**2

    K=Ke*10**(-10*thetaI)

    return K

@jit(nopython=True)
def Richards(t,psif,psie,dz,pars,const,opts,nz,qI):

    # Get hydraulic conductivity
    K=KFun(psie,psif,pars,const)
    Kmid=(K[:-1]+K[1:])/2.0

    # initialize vectors:
    q=np.zeros(nz+1)
    
    # Upper boundary: infiltration rate
    psiTop=np.minimum(psie[0],0)
    # qImax=-pars['Ks']*(psiTop/(dz[0]/2)-1)
    qImax=-K[0]*(psiTop/(dz[0]/2)-1)
    
    # q[0]=qI # 
    q[0]=np.minimum(qI,qImax)
    
    # lower boundary: free (gravity) drainage 
    q[-1]=K[-1]*opts['gravity']*opts['freeDrainage']
    
    # internal nodes
    cryoflow=-Kmid*((psif[1:]-psif[:-1])/((dz[1:]+dz[:-1])/2)-opts['gravity'])
    normflow=-Kmid*((psie[1:]-psie[:-1])/((dz[1:]+dz[:-1])/2)-opts['gravity'])
    q[1:-1]=opts['cryoflow']*cryoflow+(1-opts['cryoflow'])*normflow

    fdash=fdashFun(psie,pars)
    
    # continuity
    dthetaTdt=-(q[1:]-q[:-1])/dz       
    dpsiedt=1/fdash*dthetaTdt
    
    return dthetaTdt,dpsiedt,q

@jit(nopython=True)
def heatbalanceFun(t,psie,psif,T,TTop,TBot,jTopAdv,jTopNonAdv,dz,pars,const,opts,nz,dthetaTdt,q):
    
    # Determine the thermal cond and heat capacity for given temperature
    kappa=thermalKfun(psie,psif,T,pars,const)
    
    # Calculate the conductive heat flux:
    jd=np.zeros(nz+1)
    
    # Internal conduction fluxes using an average thermal conductivity:
    jd[1:-1]=-(kappa[1:]+kappa[:-1])/2*(T[1:]-T[:-1])/((dz[1:]+dz[:-1])/2)

    # Upper conduction boundary - no conduction (note jG comes in as advection, even if it is conduction):
    jd[0]=-kappa[0]*(T[0]-TTop)/(dz[0]/2.)*opts['conductionTop']
    jd[0] += jTopNonAdv
    
    # Lower boundary - no conduction:
    jd[-1]=-kappa[-1]*(TBot-T[-1])/(dz[-1]/2.)*opts['conductionBot']
    
    # Calculate the advective heat flux:
    ja=np.zeros(nz+1)

    # Internal (central difference ~ consider changing this)
    # ja[1:-1]=q[1:-1]*const['rho_liq']*const['cp_liq']*(T[1:]+T[:-1])/2.

    # Internal (upstream T used for advection)
    ja[1:-1]=q[1:-1]*const['rho_liq']*const['cp_liq']*T[:-1]
    ja[1:-1][q[1:-1]<0]=q[1:-1][q[1:-1]<0]*const['rho_liq']*const['cp_liq']*T[1:][q[1:-1]<0]
    
    # Upper boundary:
    ja[0]=jTopAdv 

    # Lower boundary - free drainage:
    ja[-1]=q[-1]*const['rho_liq']*const['cp_liq']*T[-1]
    
    # Putting it all together
    j=jd+opts['withadv']*ja
    
    # Heat balance terms:
    gdash=gdashfun(T,pars,const)
    CB=CBfun(psie,psif,pars,const)

    fluxDiv=-(j[1:]-j[:-1])/dz
    
    # Change in temperature in frozen conditions:
    storageTerm=(const['cp_ice']*T-const['lambda_f'])*const['rho_liq']*dthetaTdt
    denom=const['rho_liq']*gdash*(opts['massflag']*T*(const['cp_liq']-const['cp_ice'])+const['lambda_f'])+CB
    dTdt=(fluxDiv-storageTerm)/denom

    # Change in temperature in unfrozen conditions:
    storageTermUF=const['cp_liq']*T*const['rho_liq']*T*dthetaTdt
    dTdtUF=(fluxDiv-storageTermUF)/CB

    # Combine correctly:
    dTdt[psie<=psif]=dTdtUF[psie<=psif]
    
    return dTdt,j


# Model function wrapper, called by ODE solver
def ODEfun(t,DV,qI,TTop,TBot,TInf,jTopBC,dz,pars,const,opts,nz):
    return ODEfunCall(t,DV,qI,TTop,TBot,TInf,jTopBC,dz,pars,const,opts,nz)

# Model function
@jit(nopython=True)
def ODEfunCall(t,DV,qI,TTop,TBot,TInf,jTopBC,dz,pars,const,opts,nz):

    ind_psi=np.arange(nz)*2+2
    ind_T=np.arange(nz)*2+3
    psie=DV[ind_psi]
    T=DV[ind_T]
    
    psif=GCEfun(T,pars,const)
    psif=np.minimum(psie,psif)

    if opts['simulateFlow']:
        dthetaTdt,dpsiedt,q=Richards(t,psif,psie,dz,pars,const,opts,nz,qI)
    else:
        dthetaTdt=np.zeros(nz)
        dpsiedt=np.zeros(nz)
        q=np.zeros(nz+1)

    if opts['simulateTransport']:
        jTopAdv=q[0]*const['cp_liq']*const['rho_liq']*TInf
        jTopNonAdv=opts['groundHeatFlux']*jTopBC
        dTdt,j=heatbalanceFun(t,psie,psif,T,TTop,TBot,jTopAdv,jTopNonAdv,dz,pars,const,opts,nz,dthetaTdt,q)
    else:
        dTdt=np.zeros(nz)
        j=np.zeros(nz+1)
        
    dDVdt=np.zeros(2*nz+4)
    dDVdt[ind_psi]=dpsiedt
    dDVdt[ind_T]=dTdt
    dDVdt[0]=q[0]
    dDVdt[1]=j[0]
    dDVdt[-2]=q[-1]
    dDVdt[-1]=j[-1]

    return dDVdt
    
# Run model
def run(dt,t,dz,nz,T0,psi0,qI,TTop,TBot,TInf,jTopBC,parsD,const,opts,rtol=1e-7):

    constD=MakeDictFloat()
    for k in const: constD[k]=const[k]
    optsD=MakeDictFloat()
    for k in opts: optsD[k]=opts[k]
        
    nt=len(t)
    T=np.zeros((nt,nz))
    T[0,:]=T0
    psie=np.zeros((nt,nz))
    psie[0,:]=psi0
    fluxes=np.zeros((nt,4))
    
    DV = np.zeros((nt,2*nz+4))
    ind_psi=np.arange(nz)*2+2
    ind_T=np.arange(nz)*2+3
    
    DV[0,ind_psi] = psi0
    DV[0,ind_T] = T0
        
    r = ode(ODEfun)
    r.set_integrator('vode',method='BDF',uband=3,lband=3,rtol=rtol)
    
    tic=time.time()
    for i in range(len(t)-1):
        
        r.set_initial_value(np.hstack([0,0,DV[i,2:-2],0,0]), 0)
        params=(qI[i],TTop[i],TBot[i],TInf[i],jTopBC[i],dz,parsD,constD,optsD,nz)
        
        r.set_f_params(*params)
        r.integrate(dt)
        DV[i+1,:]=r.y

    runtime=time.time()-tic
    print('ode, with jac runtime = %.2f seconds'%(runtime))
    
    psie=DV[:,ind_psi]
    T=DV[:,ind_T]
    qT=DV[:,0]
    jT=DV[:,1]
    qB=DV[:,-2]
    jB=DV[:,-1]
    
    i,j=T.shape
    psif=np.array([GCEfun(T[i,:],parsD,constD) for i in range(nt)])
    psif=np.minimum(psie,psif)
    thetaL=np.array([thetaFun(psif[i,:],parsD) for i in range(nt)])
    thetaT=np.array([thetaFun(psie[i,:],parsD) for i in range(nt)])
    thetaI=const['rho_liq']/const['rho_ice']*(thetaT-thetaL)
    
    return psie,T,thetaL,thetaI,qT,qB,jT,jB
