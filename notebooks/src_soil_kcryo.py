#############################################
#
# USER SPECIFIED CONSTITUITIVE FUNCTIONS
# Maybe editted
#
#############################################
import numpy as np
from numba import njit 

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

### MODIFIED ###
@njit(inline='always')
def KFun(psie,psif,pars,const):
    # thetaL=sm.thetaFun(psif,pars)
    # thetaT=sm.thetaFun(psie,pars)
    # thetaI=1000./918.*(thetaT-thetaL)
        
    # Se=(1+(psie*-pars['alpha'])**pars['n'])**(-pars['m'])
    # Se[psie>0.]=1.0
    # Ke=pars['Ks']*Se**pars['neta']*(1-(1-Se**(1/pars['m']))**pars['m'])**2
    
    # K=Ke*10**(-10*thetaI)

    # Using frozen psi to get K:
    Se=(1+(psif*-pars['alpha'])**pars['n'])**(-pars['m'])
    Se[psif>0.]=1.0
    K=pars['Ks']*Se**pars['neta']*(1-(1-Se**(1/pars['m']))**pars['m'])**2
    
    return K

# Thermal conductivity function
@njit(inline='always')
def thermalKfun(psie,psif,T,pars,const):
    # Uses geometric mean of each component
    thetaL=thetaFun(psif,pars)
    thetaT=thetaFun(psie,pars)
    thetaI=const['rho_liq']/const['rho_ice']*(thetaT-thetaL)
    thetaA=pars['thetaS']-thetaT
    kappa = (
        (const['kappa_liq'] ** thetaL) *
        (const['kappa_ice'] ** thetaI) *
        (const['kappa_air'] ** thetaA) *
        (pars['kappa_soil'] ** pars['theta_mineral'])*
        (pars['kappa_org'] ** pars['theta_org']))
    return kappa

# Slope function dtheta_L/dT
@njit(inline='always')
def SFCslope(T,pars,const):
    psi=GCEFun(T,pars,const)
    
    C=const['lambda_f']/const['g']/const['T0']
    U=1+np.abs(pars['alpha']*C*T)**pars['n']
    dthdSe=(pars['thetaS']-pars['thetaR'])
    dSedu=-pars['m']*U**(-pars['m']-1)  
    dudT=pars['n']*(pars['alpha']*C)**pars['n']*(np.abs(T))**(pars['n']-1)
    dthdT=-dthdSe*dSedu*dudT
    dthdT[T>0]=0. 
    return dthdT
    
# GCE
@njit(inline='always')
def GCEFun(T,pars,const):
    psi=T*const['lambda_f']/(const['g']*const['T0'])
    psi[psi>0]=0.
    return psi

# Bulk heat capacity function
@njit(inline='always')
def CBFun(psie,psif,pars,const):
    thetaL=thetaFun(psif,pars)
    thetaT=thetaFun(psie,pars)
    thetaI=const['rho_liq']/const['rho_ice']*(thetaT-thetaL)
    thetaS=1-pars['thetaS']
    #thetaA=thetaS-thetaT
    CB=(const['cp_ice']*thetaI*const['rho_ice'])+(const['cp_liq']*thetaL*const['rho_liq'])+(pars['cp_soil']*thetaS*pars['rho_soil'])
    return CB


#############################################
#
# MODEL CORE FUNCTIONS - do not edit below
#
#############################################
import time
from scipy.integrate import ode

from numba import types
from numba.typed import Dict

import dill

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


@njit
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

    C=CFun(psie,pars)
    
    # continuity
    dthetaTdt=-(q[1:]-q[:-1])/dz       
    dpsiedt=1/C*dthetaTdt
    
    return dthetaTdt,dpsiedt,q

@njit
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
    Fdash=SFCslope(T,pars,const)
    CB=CBFun(psie,psif,pars,const)

    fluxDiv=-(j[1:]-j[:-1])/dz
    
    # Change in temperature in frozen conditions:
    storageTerm=(const['cp_ice']*T-const['lambda_f'])*const['rho_liq']*dthetaTdt
    denom=const['rho_liq']*Fdash*(opts['massflag']*T*(const['cp_liq']-const['cp_ice'])+const['lambda_f'])+CB
    dTdt=(fluxDiv-storageTerm)/denom

    # Change in temperature in unfrozen conditions:
    storageTermUF=const['cp_liq']*const['rho_liq']*T*dthetaTdt
    dTdtUF=(fluxDiv-storageTermUF)/CB

    # Combine correctly:
    dTdt[psie<=psif]=dTdtUF[psie<=psif]
    
    return dTdt,j


# Model function wrapper, called by ODE solver
def ODEfun(t,DV,qI,TTop,TBot,TInf,jTopBC,dz,pars,const,opts,nz):
    return ODEfunCall(t,DV,qI,TTop,TBot,TInf,jTopBC,dz,pars,const,opts,nz)

# Model function
@njit
def ODEfunCall(t,DV,qI,TTop,TBot,TInf,jTopBC,dz,pars,const,opts,nz):

    ind_psi=np.arange(nz)*2+2
    ind_T=np.arange(nz)*2+3
    psie=DV[ind_psi]
    T=DV[ind_T]
    
    psif=GCEFun(T,pars,const)
    psif=np.minimum(psie,psif)

    if opts['simulateFlow']:
        dthetaTdt,dpsiedt,q=Richards(t,psif,psie,dz,pars,const,opts,nz,qI)
    else:
        dthetaTdt=np.zeros(nz)
        dpsiedt=np.zeros(nz)
        q=np.zeros(nz+1)

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

# class for saving model output
class modelInOut:
    pass


class model:
    def __init__(self,opts=None,rtol=1e-7):
        self.opts=opts
        self.rtol=rtol

    def setOpts(self):
        
        opts={}

        # Solve correct equations (set to zero to remove Cdtheta/dt)
        if (input('Solve correct mass balance equations, y/n? (default = y)').lower() or 'y')=='y':
            opts['massflag']=1.
        else:
            opts['massflag']=0.

        # Horizontal/vertical
        if (input('Horizontal or vertical, h/v? (default = v)').lower() or 'v')=='v':
            opts['gravity']=1.
        else:
            opts['gravity']=0.
        
        # Cryosuction (yes means gradient based on psif, no means gradient based on psie)
        if (input('Include cryosuction, y/n? (default = n)').lower() or 'n')=='n':
            opts['cryoflow']=0.
        else:
            opts['cryoflow']=1.

        # Include advection
        if (input('Include advection, y/n? (default = y)').lower() or 'y')=='y':
            opts['withadv']=1.
        else:
            opts['withadv']=0.

        # Conduction at the upper boundary based on TTop
        if (input('Conduction, based on TTop, at upper boundary, y/n? (default = n)').lower() or 'n')=='n':
            opts['conductionTop']=0.
        else:
            opts['conductionTop']=1.

        # Conduction at the lower boundary based on TBot
        if (input('Conduction, based on TBot, at lower boundary, y/n? (default = n)').lower() or 'n')=='n':
            opts['conductionBot']=0.
        else:
            opts['conductionBot']=1.

        # Simulate flow
        if (input('Simulate flow, y/n? (default = y)').lower() or 'y')=='y':
            opts['simulateFlow']=1.
        else:
            opts['simulateFlow']=0.

        # Simulate heat transport
        if (input('Simulate heat transport, y/n? (default = y)').lower() or 'y')=='y':
            opts['simulateTransport']=1.
        else:
            opts['simulateTransport']=0.

        # Simulate free Drainage
        if (input('Simulate free drainage on lower boundary (n means no flow), y/n? (default = y)').lower() or 'y')=='y':
            opts['freeDrainage']=1.
        else:
            opts['freeDrainage']=0.

        self.opts=opts

    def printOpts(self):
        opts=self.opts
        print('*************************************************************************\n     SUMMARY OF MODEL OPTIONS:')
        
        print(f'     * {'Solving correct mass bal eqns' if opts['massflag']==1  else 'Solving simplified mass bal eqns'}')
        print(f'     * {'Vertical flow' if opts['gravity']==1  else 'Horizontal flow'}')
        print(f'     * {'Includes cryosuction based flow' if opts['cryoflow']==1  else 'No cryosuction'}')
        print(f'     * {'Included heat advection' if opts['withadv']==1  else 'No advection of heat'}')
        print(f'     * {'Conduction on the upper boundary' if opts['conductionTop']==1  else 'No conduction on the upper boundary'}')
        print(f'     * {'Conduction on the lower boundary' if opts['conductionBot']==1  else 'No conduction on the lower boundary'}')
        print(f'     * {'Simulating flow' if opts['simulateFlow']==1  else 'No flow'}')
        print(f'     * {'Simulating heat transport' if opts['simulateTransport']==1  else 'No heat transport'}')
        print(f'     * {'Free draining lower boundary condition' if opts['freeDrainage']==1  else 'No (mass) flow lower boundary condition'}')
        print('*************************************************************************\n')
    def setPars(self,pars):
        # Determine whether the soil parameters are uniform or not:
        uniform_pars=True
        for k in pars:
            if np.ndim(pars[k])>0:
                uniform_pars=False
        
        self.pars=pars
        
        if not(uniform_pars):
            # Non-uniform soil parameters (at least one should be given with dim nz):
            for k in pars:
                if np.ndim(pars[k])==0:
                    # Here we have one parameter given, so duplicate it to make it uniform
                    self.pars[k]=np.full(self.nz,pars[k])
                elif len(pars[k])==1:
                    # Here we have one parameter given, so duplicate it to make it uniform
                    self.pars[k]=np.full(self.nz,pars[k])
                elif len(pars[k])==2:
                    # Here we have two parameters defining a linear relationship of the form
                    self.pars[k]=np.interp(self.z,[self.z0,self.zMax],pars[k]) 
                elif len(pars[k])==3:
                    # Exponential relationship where:
                    # pars[k][0]=parameter value at ground surface
                    # pars[k][1]=parameter value at 1m depth
                    # pars[k][2]=parameter value at infinity
                    c=np.log((pars[k][0]-pars[k][2])/(pars[k][1]-pars[k][2]))
                    self.pars[k]=pars[k][2]+(pars[k][0]-pars[k][2])*np.exp(-c*self.z)
        
    def setConst(self,const):
        self.const=const

    def zGrid(self,bz):
        self.dz=np.diff(bz)
        self.z=bz[:-1]+self.dz/2
        self.nz=len(self.z)
        self.z0=bz[0]
        self.zMax=bz[-1]
        
    def tGrid(self,t0,tMax,dt):
        self.dt=dt
        self.t=np.arange(0,tMax+dt,dt)
        self.nt=len(self.t)
        
    def BCs(self,jTopBC=0,qI=0,TInf=0,TTop=0,TBot=0):
        
        def _to_timeseries(x):
            if np.ndim(x) == 0:        
                return np.zeros(self.nt)+x
            else:
                return x

        self.jTopBC = _to_timeseries(jTopBC)
        self.qI     = _to_timeseries(qI)
        self.TInf   = _to_timeseries(TInf)
        self.TTop   = _to_timeseries(TTop)
        self.TBot   = _to_timeseries(TBot)

    def ICs(self,T0,psi0):

        def _to_array(x):
            if np.ndim(x) == 0:
                # Uniform initial condition based on scalar
                return np.zeros(self.nz)+x
            elif len(x) == 1:
                # Uniform initial condition based on array with just one value
                return np.zeros(self.nz)+x
            elif len(x) == 2:
                # Linear initial condition from top to bottom:
                return np.interp(self.z,[self.z0,self.zMax],x)
            elif len(x)==self.nz:
                # Fully specified 
                return x

        self.T0=_to_array(T0)
        self.psi0=_to_array(psi0)
        
    # Run model
    def run(self):

        uniform_pars=True
        for k in self.pars:
            if np.ndim(self.pars[k])>0:
                uniform_pars=False
        
        if uniform_pars:
            # Uniform soil parameters (more efficient):
            parsD=MakeDictFloat()
        else:
            # Non-uniform soil parameters (at least one should be given with dim nz):
            parsD=MakeDictArray()
            
        for k in self.pars: parsD[k]=self.pars[k]
            
        constD=MakeDictFloat()
        for k in self.const: constD[k]=self.const[k]
        optsD=MakeDictFloat()
        for k in self.opts: optsD[k]=self.opts[k]
        
        T=np.zeros((self.nt,self.nz))
        T[0,:]=self.T0
        psie=np.zeros((self.nt,self.nz))
        psie[0,:]=self.psi0
        fluxes=np.zeros((self.nt,4))
        
        DV = np.zeros((self.nt,2*self.nz+4),dtype=np.float64)
        ind_psi=np.arange(self.nz)*2+2
        ind_T=np.arange(self.nz)*2+3
        
        DV[0,ind_psi] = self.psi0
        DV[0,ind_T] = self.T0
            
        r = ode(ODEfun)
        r.set_integrator('vode',method='BDF',uband=3,lband=3,rtol=self.rtol)
        
        tic=time.time()
        for i in range(self.nt-1):
            
            r.set_initial_value(np.hstack([0,0,DV[i,2:-2],0,0]), 0)
            params=(self.qI[i],self.TTop[i],self.TBot[i],self.TInf[i],self.jTopBC[i],self.dz,parsD,constD,optsD,self.nz)
            
            r.set_f_params(*params)
            r.integrate(self.dt)
            DV[i+1,:]=r.y
    
        runtime=time.time()-tic
        
        print('*************************************************************************')
        print('     soilice ran successfully!')
        print('               runtime:     % .2f seconds'%(runtime))        
        # print(' ode, with jac runtime:     % .2f seconds'%(runtime))

        # Pack up model input/output:
        inOut=modelInOut()
        for i in ['dt','nt','t','dz','nz','z','z0','zMax',
                  'pars','const',
                  'jTopBC','qI','TInf','TTop','TBot',
                  'T0','psi0','opts','rtol']: 
            setattr(inOut, i, getattr(self, i))
        
        inOut.psie=DV[:,ind_psi]
        inOut.T=DV[:,ind_T]
        inOut.qT=DV[:,0]
        inOut.jT=DV[:,1]
        inOut.qB=DV[:,-2]
        inOut.jB=DV[:,-1]
        
        i,j=T.shape
        inOut.psif=np.array([GCEFun(inOut.T[i,:],parsD,constD) for i in range(self.nt)])
        inOut.psif=np.minimum(inOut.psie,inOut.psif)
        inOut.thetaL=np.array([thetaFun(inOut.psif[i,:],parsD) for i in range(self.nt)])
        inOut.thetaT=np.array([thetaFun(inOut.psie[i,:],parsD) for i in range(self.nt)])
        inOut.thetaI=self.const['rho_liq']/self.const['rho_ice']*(inOut.thetaT-inOut.thetaL)

        # Print mass/energy balance errors:
        self.balanceClosure(inOut)
        print('*************************************************************************\n')
        return inOut

    def balanceClosure(self,inOut):

        nt,nz=inOut.thetaL.shape
        ml=inOut.thetaL*inOut.dz*inOut.const['rho_liq']
        mi=inOut.thetaI*inOut.dz*inOut.const['rho_ice']
        ms=np.zeros((nt,nz))+((1-inOut.pars['thetaS'])*inOut.dz*inOut.pars['rho_soil'])
        u=(ml*inOut.const['cp_liq']+mi*inOut.const['cp_ice']+ms*inOut.pars['cp_soil'])*inOut.T-mi*inOut.const['lambda_f']
        ml=np.sum(ml,axis=1)
        mi=np.sum(mi,axis=1)
        u=np.sum(u,axis=1)
        du=u[-1]-u[0]
        m=ml+mi
        dm=m[-1]-m[0]
        qT=inOut.qT.sum()*inOut.const['rho_liq']
        qB=inOut.qB.sum()*inOut.const['rho_liq']
        jT=inOut.jT.sum()
        jB=inOut.jB.sum()
        print(f'    Mass balance error: {(qT-qB)/dm*100-100: .2e} %')
        print(f'                        {(qT-qB-dm): .2e} kg')
        print(f'  Energy balance error: {(jT-jB)/du*100-100: .2e} %')
        print(f'                        {(jT-jB-du): .2e} J')
        
        inOut.u=u
        inOut.m=m

    