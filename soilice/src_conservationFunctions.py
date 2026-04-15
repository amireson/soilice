import numpy as np
from numba import njit

try:
    from src_constitutiveFunctions import thetaFun, CFun, KFun, thermalKfun, SFCslope, GCEFun, CBFun
except ImportError:
    from soilice.src_constitutiveFunctions import thetaFun, CFun, KFun, thermalKfun, SFCslope, GCEFun, CBFun

#############################################
#
# CONSERVATION EQUATIONS (edit with 
# extreme caution)
#
#############################################

@njit(inline='always')
def Richards(t,psif,psie,dz,pars,const,opts,nz,upperBC):

    # Get hydraulic conductivity
    K=(
        KFun(psif,psif,pars,const)*opts['cryoK']+
        KFun(psie,psif,pars,const)*(1-opts['cryoK'])
    )

    Kmid=(K[:-1]+K[1:])/2.0
    Ksurf=KFun(np.array([upperBC]),np.array([upperBC]),pars,const)

    # initialize vectors:
    q=np.zeros(nz+1)
    
    # Upper boundary: infiltration rate
    qImax=-K[0]*(psie[0]/(dz[0]/2)-1)
    #psiTop=np.minimum(psie[0],0)
    #qImax=-K[0]*(psiTop/(dz[0]/2)-1)
    
    # q[0]=qI # 
    q[0]=np.minimum(upperBC,qImax)*opts['infiltration']
    #q[0]+= -K[0]*((psie[0]-upperBC)/(dz[0]/2)-1)*(1-opts['infiltration'])
    q[0]+= -(Ksurf[0]+K[0])/2*((psie[0]-upperBC)/(dz[0]/2)-1)*(1-opts['infiltration'])
    
    # lower boundary: free (gravity) drainage 
    q[-1]=K[-1]*opts['gravity']*opts['freeDrainage']
    
    # internal nodes
    cryoflow=-Kmid*((psif[1:]-psif[:-1])/((dz[1:]+dz[:-1])/2)-opts['gravity'])
    normflow=-Kmid*((psie[1:]-psie[:-1])/((dz[1:]+dz[:-1])/2)-opts['gravity'])
    q[1:-1]=opts['cryoGradient']*cryoflow+(1-opts['cryoGradient'])*normflow

    C=CFun(psie,pars)
    
    # continuity
    dthetaTdt=-(q[1:]-q[:-1])/dz       
    dpsiedt=1/C*dthetaTdt
    
    return dthetaTdt,dpsiedt,q

@njit(inline='always')
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
    Fdash=SFCslope(psie,psif,pars,const)
    CB=CBFun(psie,psif,pars,const)

    fluxDiv=-(j[1:]-j[:-1])/dz
    
    # Change in temperature in frozen conditions:
    storageTerm=(const['cp_ice']*T-const['lambda_f'])*const['rho_liq']*dthetaTdt
    denom=const['rho_liq']*Fdash*(T*(const['cp_liq']-const['cp_ice'])+const['lambda_f'])+CB
    dTdt=(fluxDiv-storageTerm)/denom

    # Change in temperature in unfrozen conditions:
    storageTermUF=const['cp_liq']*const['rho_liq']*T*dthetaTdt
    dTdtUF=(fluxDiv-storageTermUF)/CB

    # Combine correctly:
    dTdt[psie<=psif]=dTdtUF[psie<=psif]
    
    return dTdt,j

@njit(inline='always')
def ODEfunCall(t,DV,upperBC,TTop,TBot,TInf,jTopBC,dz,pars,const,opts,nz):

    ind_psi=np.arange(nz)*2+2
    ind_T=np.arange(nz)*2+3
    psie=DV[ind_psi]
    T=DV[ind_T]
    
    psif=GCEFun(T,pars,const)
    psif=np.minimum(psie,psif)

    if opts['simulateFlow']:
        dthetaTdt,dpsiedt,q=Richards(t,psif,psie,dz,pars,const,opts,nz,upperBC)
    else:
        dthetaTdt=np.zeros(nz)
        dpsiedt=np.zeros(nz)
        q=np.full(nz+1,pars['q'][0] if hasattr(pars['q'], "__getitem__") else pars['q'])

    if opts['simulateTransport']:
        jTopAdv=q[0]*const['cp_liq']*const['rho_liq']*TInf
        jTopNonAdv=jTopBC
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

