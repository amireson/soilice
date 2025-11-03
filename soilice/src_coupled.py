# This script loosely couples the snow and soil models. That is, it solves each
# separately at the time step level.

import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
import time 
from scipy.integrate import ode
from numba import jit 

from .src_soil import MakeDictArray
from .src_soil import MakeDictFloat

from .src_soil import ODEfun as soilODEfun
from .src_soil import GCEfun
from .src_soil import thetaFun

from .src_snow import ODEfun as snowODEfun
from .src_snow import uFun
from .src_snow import TFun
        
def run(df,dz,nz,pars,constIN,optsIN,T0,psi0,msini,Tsini,TBot,rtol=1e-7):

    const=MakeDictFloat()
    for k in constIN: const[k]=constIN[k]
    opts=MakeDictFloat()
    for k in optsIN: opts[k]=optsIN[k]
        
    # Define time grid:
    t=(df.index-df.index[0]).total_seconds() # Time in seconds
    dt=(t[1]-t[0])                           # Timestep in seconds
    nt=len(t)
    P=df['P'].values                         # Precipitation rate (kg/m2/s)
    Ta=df['T'].values                        # Air temperature (deg C)
    SWin=df['SW'].values                     # Radiation (J/m2/s)
    LWin=df['LW'].values                     # Radiation (J/m2/s)
    U=df['U'].values                         # Windspeed (m/s)
    BP=df['BP'].values                       # Barometric pressure (Pa)
    SH=df['SH'].values                       # Specific humidity (kg/kg)

    # Dependent variable array, with initial conditions
    DVsoil = np.zeros((nt,2*nz+4))
    ind_psi=np.arange(nz)*2+2
    ind_T=np.arange(nz)*2+3
    
    DVsoil[0,ind_psi] = psi0
    DVsoil[0,ind_T] = T0
    
    DVsnow=np.zeros((nt,16))
    DVsnow[0,0]=msini          # Initial snow ice mass
    DVsnow[0,1]=0.             # Initial snow liquid mass
    DVsnow[0,2]=uFun(msini,0,Tsini,const) # Initial snow enthalpy
    DVsnow[0,3:]=0.            # Cumulative fluxes
    
    # Snow initial density and depth
    rhoS=np.zeros(len(t))
    rhoS[0]=pars['rho_newSnow']
    zs=msini/pars['rho_newSnow']
    
    potentialInfiltration=np.zeros(len(t))

    thetaLS = 0.   # Used for soil evaporation calculation.
    dz_soil=dz[0]  # Thickness of the top soil layer
        
    rsnow = ode(snowODEfun)
    snowtol = np.concatenate([np.full(3, 1e-6),
                       np.full(13, 1e-1)])
    rsnow.set_integrator('vode',method='BDF',rtol=snowtol,atol=snowtol,first_step=0.01)

    rsoil = ode(soilODEfun)
    rsoil.set_integrator('vode',method='BDF',uband=3,lband=3,rtol=rtol)
    
    tic=time.time()
    for i in range(1,len(t)):
        
        # Solve snow:
        TSS=DVsoil[i-1,3]      # Soil temperature from previous time step
        
        rsnow.set_initial_value(DVsnow[i-1,:], 0)
        params_snow=(pars,const,SWin[i],LWin[i],P[i],Ta[i],U[i],BP[i],SH[i],TSS,rhoS[i-1],zs,dz_soil)
        rsnow.set_f_params(*params_snow)
        
        rsnow.integrate(dt)
        DVsnow[i,:]=rsnow.y
        
        # Calculate new snow density: first average with new snow:
        if DVsnow[i,0]>DVsnow[i-1,0]:             # If the snowpack got deeper, average the fresh snow:
            ms_new=DVsnow[i,0]-DVsnow[i-1,0]
            Delta_zs=ms_new/pars['rho_newSnow']
            zs_old=DVsnow[i-1,0]/rhoS[i-1]
            zs=zs_old+Delta_zs
            rhoS[i]=DVsnow[i,0]/zs
        else:                             # Otherwise keep the snowpack the same
            rhoS[i]=rhoS[i-1]
        
        if DVsnow[i,0]>0.1:                   # If there is snow on the ground, age the snow
            # Calculate new snow density: next age the snowpack:
            rhoS[i]=(rhoS[i]-pars['rho_maxSnow'])*np.exp(-pars['snow_aging_coefficient']*dt)+pars['rho_maxSnow']
        else:                             # Otherwise reset to the density of fresh snow
            rhoS[i]=pars['rho_newSnow']

        # Prepare inputs for soil model:
        qR_ground=(DVsnow[i,8]-DVsnow[i-1,8])/dt            # rain on the ground
        qD=(DVsnow[i,5]-DVsnow[i-1,5])/dt                   # snow drainage
        qI=(qR_ground+qD)/const['rho_liq']  # Potential Infiltration flux

        if qD>0:
            TInf=0.                 # Snow drainage is assumed to have a temperature of zero
        else:
            TInf=Ta[i]              # Rainfall on the ground has temperature equal to air temperature
            
        jG=(DVsnow[i,10]-DVsnow[i-1,10])/dt              # Soil heat flux due to conduction with snow
        jRn_ground=(DVsnow[i,15]-DVsnow[i-1,15])/dt      # Soil heat flux due to ground radiation

        jG=jG+jRn_ground
        
        # Solve soil:
        rsoil.set_initial_value(DVsoil[i-1,:], 0)
        params_soil=(qI,0,TBot[i],TInf,jG,dz,pars,const,opts,nz)
        
        rsoil.set_f_params(*params_soil)
        rsoil.integrate(dt)
        
        DVsoil[i,:]=rsoil.y
        psie0=DVsoil[i,2]
        thetaLS=thetaFun(np.array([psie0]),pars)
            
    runtime=time.time()-tic
    print('ode, with jac runtime = %.2f seconds'%(runtime))

    # Unpack output:

    # Convert soil fluxes from non-cumulative total to timestep totals:
    DVsoil[1:,[0,1,-2,-1]]=np.diff(DVsoil[:,[0,1,-2,-1]],axis=0)
    
    psie=DVsoil[:,ind_psi]
    T=DVsoil[:,ind_T]
    qST=DVsoil[:,0]
    jST=DVsoil[:,1]
    qSB=DVsoil[:,-2]
    jSB=DVsoil[:,-1]

    i,j=T.shape
    psif=np.reshape(GCEfun(T.flatten(),pars,const),(i,j))
    psif=np.minimum(psie,psif)
    thetaL=np.reshape(thetaFun(psif.flatten(),pars),(i,j))
    thetaT=np.reshape(thetaFun(psie.flatten(),pars),(i,j))
    thetaI=const['rho_liq']/const['rho_ice']*(thetaT-thetaL)


    # Convert snow fluxes from non-cumulative total to timestep totals:
    DVsnow[1:,3:]=np.diff(DVsnow[:,3:],axis=0)
    
    # Pack output:
    snow_df=pd.DataFrame(index=df.index)
    snow_df['mi']=DVsnow[:,0]          # Snow ice mass (kg/m2)
    snow_df['ml']=DVsnow[:,1]          # Snow water mass (kg/m2)
    snow_df['u']=DVsnow[:,2]           # Snow temperature (deg C)
    snow_df['qR_snow']=DVsnow[:,3]          # Cumulative rainfall (kg/m2)
    snow_df['qS']=DVsnow[:,4]          # Cumulative snowfall (kg/m2)
    snow_df['qD']=DVsnow[:,5]          # Cumulative snowmelt drainage (kg/m2)
    snow_df['qE_snow']=DVsnow[:,6]          # Cumulative evaporation (kg/m2)
    snow_df['qM']=DVsnow[:,7]          # Cumulative melt (kg/m2)
    snow_df['qR_ground']=DVsnow[:,8]        # Cumulative rain on ground (kg/m2)
    snow_df['jRn_snow']=DVsnow[:,9]         # Cumulative net radiation (J/m2)
    snow_df['jG']=DVsnow[:,10]         # Cumulative heat flux into the soil (J/m2)
    snow_df['jE_snow']=DVsnow[:,11]         # Cumulative latent heat flux (J/m2)
    snow_df['jH_snow']=DVsnow[:,12]         # Cumulative sensible heat flux (J/m2)
    snow_df['jAPS']=DVsnow[:,13]       # Cumulative advection from snowfall (J/m2)
    snow_df['jAPR']=DVsnow[:,14]       # Cumulative advection from rainfall (J/m2)
    snow_df['jRn_ground']=DVsnow[:,15]         # Cumulative advection from drainage (J/m2)
    snow_df['rhoS']=rhoS
    snow_df['d']=df.index
    snow_df['Ta']=df['T']
    
    # Temperature
    mi=DVsnow[:,0]          # Snow ice mass (kg/m2)
    ml=DVsnow[:,1]          # Snow water mass (kg/m2)
    u=DVsnow[:,2]           # Snow temperature (deg C)
    Ts=TFun(mi,ml,u,const)
    Ts[mi<pars['zeroMass']]=np.nan
    snow_df['Ts']=Ts

    return t,snow_df,qST,jST,qSB,jSB,psie,T,thetaL,thetaI
    ###################################





 
    
