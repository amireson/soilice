import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
import time
from IPython.display import clear_output

from scipy.integrate import ode #, solve_ivp
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

# Functions to define u and T
@jit(nopython=True)
def uFun(mi,ml,T,const):
    u=mi*(const['cp_ice']*T-const['lambda_f'])+ml*const['cp_liq']*T
    return u

@jit(nopython=True)
def TFun(mi,ml,u,const):
    numer=(u+mi*const['lambda_f'])
    denom=(mi*const['cp_ice']+ml*const['cp_liq'])
    T=numer/denom
    return T

# Snow thermal conductivity functions:

# Mellor (1997)
@jit(nopython=True)
def Mellor(rhoS):
    ks = 2.576e-6*(rhoS**2)+0.074
    return ks

# Essery
@jit(nopython=True)
def Essery(rhoS, pars,const):
    ks = pars['Ki']*((rhoS/const['rho_ice'])**pars['b'])
    return ks

# Sturm (1997)
@jit(nopython=True)
def Sturm(rhoS):
    rhoS_u=rhoS/1000   # rhoS in g/cm3
    if rhoS >= 156:
        ks = 3.233*(rhoS_u**2)-(1.01*rhoS_u)+0.138
    else:
        ks = 0.234*rhoS_u+0.023
    return ks

@jit(nopython=True)
def snowThermalConductivity(rhos,pars,const):
#     ks=Mellor(rhos)
    ks=Essery(rhos,pars,const)
    # ks=Sturm(rhos)
#     ks=0.
    # ks=0.2
    # print(ks,rhos)
    return ks


# Function to calculate turbulent fluxes 
@jit(nopython=True)
def turbulentFluxes(Ts,Ta,U,BP,SH,const,pars,fs,t):
    
    # After Essery, 2015
    qE=0.  
    jE=0.
    jH=0.
    
    # Transfer coefficient
    z0=pars['z0s']**fs*pars['z0g']**(1-fs)
    zU=1.5
    zTQ=1.5
    z0h=0.1*z0
    k=0.4 
    bh=5.
    Rair=287.
    RiB=9.81*zU**2*(Ta-Ts)/(zTQ*(Ta+273.15)*U**2)
    if RiB>=0:
        fH=(1+3*bh*RiB*(1+bh*RiB)**0.5)**-1
    else:
        c=3*bh**2*k**2*(zU/z0)**0.5*(np.log(zU/z0))**-2
        fH=1-3*bh*RiB*(1+c*(-RiB)**0.5)**-1
    
    CH=fH*k**2*(np.log(zU/z0)*np.log(zTQ/z0h))**-1

    # Air density
    rho_air=BP/(Rair*(Ta+273.15))

    # Sensible heat
    jH=rho_air*const['cp_air']*CH*U*(Ts-Ta) # x
    
    # Saturated specific humidity
    eps=0.622
    e0=611.213

    TsnowSurface=(Ta+Ts)/2.
    if Ts>const['Tf']:
        # Water:
        es = e0*np.exp(17.5043*Ts / (241.3 + Ts))
        L=const['lambda_v']
        # h_per_m=uFun(0,1,Ts,const)
#         h_per_m=const['cp_liq']*Ts+const['']
    else:
        # Ice:
        es = e0*np.exp(22.4422*Ts / (272.186 + Ts))
        L=const['lambda_s']
        # h_per_m=uFun(1,0,Ts,const)
#         h_per_m=const['cp_ice']*Ts-const['lambda_f']
        
    Qsat=eps*es/BP

    # Latent heat:
    qE=rho_air*CH*U*(Qsat-SH)
    # print(jE)
    jE=jE*L

    # # Fugde with this:
    # jE=jH/2.
    # qE=jE/L
    
#     if SH>Qsat: SH=Qsat
#     # Latent heat
# #     moisture_availability_factor = g/(g+CH*U)
#     qE=rho_air*CH*U*(Qsat-SH)
#     jE=qE*(h_per_m+L)
    
# #     if Ts>const['Tf']:
# #         jE=qE*(const['cp_ice']*const['Tf']+const['cp_liq']*(Ts-const['Tf'])+const['lambda_f']+const['lambda_v'])
# #     else:
# #         jE=qE*(const['cp_ice']*Ts+const['lambda_v'])
#     qE=0.
#     jE=0.
#     jH=0.
    
    return qE,jE,jH

# Function to partition precipitation
@jit(nopython=True)
def precip_type(qP, Ta):
    if Ta > 0.:
        qS = 0
        qR = qP
    else:
        qS = qP
        qR = 0
        
    return qS,qR

# Function for net radiation from snowpack
@jit(nopython=True)
def netRadSnow(pars,const,SWin,LWin,T):
    LWout=pars['emissivitySnow']*const['stefanBoltzmann']*(T+273.15)**4
    Rn=SWin*(1-pars['albedoSnow'])+LWin-LWout
    return Rn

# Functions to calculate ground heat fluxes
@jit(nopython=True)
def groundHeatFluxSnow(Ts,TSS,rhoS,zs,pars,const,dz_soil):
    ks=snowThermalConductivity(rhoS,pars,const)
    jG=ks*(Ts-TSS)/((zs/2+dz_soil/2))
    # jG=ks*(Ts-TSS)/(zs/2)
    return jG

@jit(nopython=True)
def groundHeatFluxNoSnow(TSS,pars,const,SWin,LWin,jH,jE,jAPS,jAPR):
    # Net radiation
    LWout=pars['emissivitySoil']*const['stefanBoltzmann']*(TSS+273.15)**4
    jRn=SWin*(1-pars['albedoSoil'])+LWin-LWout
    jG=jRn+jAPS+jAPR-jH-jE
    return jG,jRn

# odefun
def ODEfun(t,DV,pars,const,SWin,LWin,P,Ta,U,BP,SH,TSS,rhoS,zs,dz_soil):
    return ODEfunCall(t,DV,pars,const,SWin,LWin,P,Ta,U,BP,SH,TSS,rhoS,zs,dz_soil)

@jit(nopython=True)
def ODEfunCall(t,DV,pars,const,SWin,LWin,P,Ta,U,BP,SH,TSS,rhoS,zs,dz_soil):
    
    # Ta is air temperature
    # Ts is snow temperature
    # TSS is soil surface temperature
    
    # Unpack
    mi=DV[0]
    ml=DV[1]
    u=DV[2]
    
    # Partition precipitation
    qS,qR=precip_type(P,Ta)

    # Snow drainage:
    qD=pars['drainage_coeff']*ml
    jD=0 #hFun(0,qD,Ts,const)
    
    if (mi>pars['zeroMass']) | (qS>0): # Yes snow on ground (fresh or old)

        qR_snow=qR
        qR_ground=0.
        
        jAPS=uFun(qS,0,Ta,const)
        jAPR=uFun(0,qR,Ta,const)
        
        # Calculate temperature from Enthalpy
        if u<uFun(mi,0,0,const):
            Ts=TFun(mi,ml,u,const)
        else:
            Ts=0.

        # Radiation
        jRn=netRadSnow(pars,const,SWin,LWin,Ts)                
        
        # Turbulent fluxes
        qE,jE,jH=turbulentFluxes(Ts,Ta,U,BP,SH,const,pars,1.,t) 
        qE=qE*mi/(1+mi) # Monod to prevent excessive E
        jE=jE*mi/(1+mi) # Monod to prevent excessive E
        
        # Ground heat flux
        jG=groundHeatFluxSnow(Ts,TSS,rhoS,zs,pars,const,dz_soil)       
        
        # Sum of heat fluxes
        zeta=jRn+jAPS+jAPR-jH-jE-jG-jD                            
        
        # If snow is cold, warm/cool the snowpack, no phase change other than freezing rain:
        if u<uFun(mi,0,0,const):
            qM=-qR  # This says, add rainfall to the ice in in the snow, as negative melt
                    # but if there is no rain, there is no melt.
            
        # Otherwise snow is at zero - heat goes into phase change
        else:
            potentialMelt=zeta/const['lambda_f']+qS-qE
            if (potentialMelt>=0) | (ml>0):
                qM=potentialMelt
            else:
                qM=0.
        dmidt=qS-qM-qE 
        
    else: # No snow on the ground     
        qR_ground=qR
        qR_snow=0.
        qE=0.
        jE=0.
        jH=0.
        jAPS=0.
        jAPR=0.
        qE=0.
        jE=0.
        jRn=0.        
        dum,jRn_ground=groundHeatFluxNoSnow(TSS,pars,const,SWin,LWin,0,0,0,0)
        # Sum of heat fluxes (with no snow on the ground, the only heat flux
        # is associated with draining residual liquid water)
        zeta=0.
        # Snowmelt:
        qM=0.
        
        dmidt=0
        
    # Solve balance equations:
    dmldt=qR_snow+qM-qD
    dudt=zeta
    
    # Pack up
    dDVdt=np.array([dmidt,dmldt,dudt,qR_snow,qS,qD,qE,qM,qR_ground,jRn,jG,jE,jH,jAPS,jAPR,jRn_ground])      
    return dDVdt
    
# solver
def runModel(df,parsIN,constIN,msini,Tsini,TSS,dz_soil):
    # solver(dt,t,parsIN,constIN,msini,Tsini,P,Ta,SWin,LWin,U,BP,SH,TSS):
    pars=MakeDictFloat()
    for k in parsIN: pars[k]=parsIN[k]

    const=MakeDictFloat()
    for k in constIN: const[k]=constIN[k]
    
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

    DV=np.zeros((len(t),16))
#     DV=np.zeros((len(t),3))
    DV[0,0]=msini          # Initial snow ice mass
    DV[0,1]=0.             # Initial snow liquid mass
    DV[0,2]=uFun(msini,0,Tsini,const) # Initial snow enthalpy
    DV[0,3:]=0.            # Cumulative fluxes
        
    rhoS=np.zeros(len(t))
    rhoS[0]=parsIN['rho_newSnow']
    zs=msini/parsIN['rho_newSnow']
    
    potentialInfiltration=np.zeros(len(t))

    r = ode(ODEfun)
    
    rtol = np.concatenate([np.full(3, 1e-6),
                       np.full(13, 1e-1)])
    r.set_integrator('vode',method='BDF',rtol=rtol,atol=rtol,first_step=0.01)
    

    tic=time.time()
    for i in range(1,nt):
        r.set_initial_value(DV[i-1,:], 0)
        params=(pars,const,SWin[i],LWin[i],P[i],Ta[i],U[i],BP[i],SH[i],TSS,rhoS[i-1],zs,dz_soil)
        r.set_f_params(*params)
        r.integrate(dt)
         
        DV[i,:]=r.y
 
        # Calculate new snow density: first average with new snow:
        if DV[i,0]>DV[i-1,0]:             # If the snowpack got deeper, average the fresh snow:
            ms_new=DV[i,0]-DV[i-1,0]
            Delta_zs=ms_new/pars['rho_newSnow']
            zs_old=DV[i-1,0]/rhoS[i-1]
            zs=zs_old+Delta_zs
            rhoS[i]=DV[i,0]/zs
        else:                             # Otherwise keep the snowpack the same
            rhoS[i]=rhoS[i-1]
        
        if DV[i,0]>0.1:                   # If there is snow on the ground, age the snow
            # Calculate new snow density: next age the snowpack:
            rhoS[i]=(rhoS[i]-pars['rho_maxSnow'])*np.exp(-pars['snow_aging_coefficient']*dt)+pars['rho_maxSnow']
        else:                             # Otherwise reset to the density of fresh snow
            rhoS[i]=pars['rho_newSnow']
 
    runtime=time.time()-tic
    print('ode, with jac runtime = %.2f seconds'%(runtime))

    # Convert fluxes from non-cumulative total to timestep totals:
    DV[1:,3:]=np.diff(DV[:,3:],axis=0)
    
    # Pack output:
    snow_df=pd.DataFrame(index=df.index)
    snow_df['mi']=DV[:,0]          # Snow ice mass (kg/m2)
    snow_df['ml']=DV[:,1]          # Snow water mass (kg/m2)
    snow_df['u']=DV[:,2]           # Snow temperature (deg C)
    snow_df['qR_snow']=DV[:,3]          # Cumulative rainfall (kg/m2)
    snow_df['qS']=DV[:,4]          # Cumulative snowfall (kg/m2)
    snow_df['qD']=DV[:,5]          # Cumulative snowmelt drainage (kg/m2)
    snow_df['qE_snow']=DV[:,6]          # Cumulative evaporation (kg/m2)
    snow_df['qM']=DV[:,7]          # Cumulative melt (kg/m2)
    snow_df['qR_ground']=DV[:,8]        # Cumulative rain on ground (kg/m2)
    snow_df['jRn_snow']=DV[:,9]         # Cumulative net radiation (J/m2)
    snow_df['jG']=DV[:,10]         # Cumulative heat flux into the soil (J/m2)
    snow_df['jE_snow']=DV[:,11]         # Cumulative latent heat flux (J/m2)
    snow_df['jH_snow']=DV[:,12]         # Cumulative sensible heat flux (J/m2)
    snow_df['jAPS']=DV[:,13]       # Cumulative advection from snowfall (J/m2)
    snow_df['jAPR']=DV[:,14]       # Cumulative advection from rainfall (J/m2)
    snow_df['jRn_ground']=DV[:,15]         # Cumulative advection from drainage (J/m2)
    snow_df['rhoS']=rhoS
    snow_df['d']=df.index
    snow_df['Ta']=df['T']
    
    # Temperature
    mi=DV[:,0]          # Snow ice mass (kg/m2)
    ml=DV[:,1]          # Snow water mass (kg/m2)
    u=DV[:,2]           # Snow temperature (deg C)
    Ts=TFun(mi,ml,u,const)
    Ts[mi<pars['zeroMass']]=np.nan
    snow_df['Ts']=Ts
    snow_df['zs']=mi/rhoS
    
    return snow_df