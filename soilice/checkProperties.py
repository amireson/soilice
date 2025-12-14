import numpy as np
import matplotlib.pyplot as pl
import pandas as pd
from numba import jit 
from numba import types
from numba.typed import Dict

from .src_soil import GCEfun
from .src_soil import SFCslope
from .src_soil import CBfun
from .src_soil import thermalKfun
from .src_soil import thetaFun
from .src_soil import CFun
from .src_soil import KFun
from .src_soil import MakeDictFloat

def getProperties(psie,T,pars,const):
    
    orig_shape=psie.shape
    
    psie=psie.ravel()
    T=T.ravel()
    
    parsD=MakeDictFloat()
    for key in pars: parsD[key]=pars[key]

    constD=MakeDictFloat()
    for key in const: constD[key]=const[key]
        
    # equivalent and frozen matric potential:
    psif=GCEfun(T,parsD,constD)
    psif=np.minimum(psie,psif)

    # partition theta:
    thetaL=thetaFun(psif,parsD)
    thetaT=thetaFun(psie,parsD)
    thetaI=const['rho_liq']/const['rho_ice']*(thetaT-thetaL)

    # properties:
    dthdT=SFCslope(T,parsD,constD)
    CB=CBfun(psie,psif,parsD,constD)
    kappa=thermalKfun(psie,psif,T,parsD,constD)
    C=CFun(psie,parsD)
    Fdash=SFCslope(T,parsD,constD)
    Kf=KFun(psie,psif,parsD,constD)
    Ke=KFun(psie,psie,parsD,constD)

    psie=psie.reshape(orig_shape)
    T=T.reshape(orig_shape)
    psif=psif.reshape(orig_shape)
    thetaL=thetaL.reshape(orig_shape)
    thetaI=thetaI.reshape(orig_shape)
    thetaT=thetaT.reshape(orig_shape)
    dthdT=dthdT.reshape(orig_shape)
    CB=CB.reshape(orig_shape)
    kappa=kappa.reshape(orig_shape)
    C=C.reshape(orig_shape)
    Fdash=Fdash.reshape(orig_shape)
    Kf=Kf.reshape(orig_shape)
    Ke=Ke.reshape(orig_shape)
    
    return psie,psif,thetaL,thetaI,thetaT,dthdT,CB,kappa,C,Fdash,Kf,Ke

def plotProperties(pars,const):

    pl.figure(figsize=(10,8))
    n=1001
    psi=np.linspace(-150,1,n)
    i=-1
    mycolor = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
    for Ti in [-10,-1,-0.1,-0.01,0]:
        i+=1
        
        T=np.linspace(Ti,Ti,n)
        psie,psif,thetaL,thetaI,thetaT,dthdT,CB,kappa,C,Fdash,Kf,Ke=getProperties(psi,T,pars,const)
        
        pl.subplot(4,2,1)
        pl.semilogx(-psi,thetaL,color=mycolor[i],label=f'T={Ti}')
        pl.ylabel('theta_L'); pl.grid('on')
        pl.subplot(4,2,2)
        pl.semilogx(-psi,thetaI,color=mycolor[i])
        pl.ylabel('theta_I'); pl.grid('on')
        pl.subplot(4,2,3)
        pl.semilogx(-psi,CB/1e6,mycolor[i])
        pl.ylabel('C_B 1e6'); pl.grid('on')
        pl.subplot(4,2,4)
        pl.semilogx(-psi,kappa,mycolor[i])
        pl.ylabel('kappa'); pl.grid('on')
        pl.subplot(4,2,5)
        pl.semilogx(-psi,C,mycolor[i])
        pl.ylabel('C'); pl.grid('on')
        pl.subplot(4,2,6)
        pl.semilogx(-psi,Fdash,mycolor[i])
        pl.ylabel('Fdash'); pl.grid('on')
        pl.subplot(4,2,7)
        pl.loglog(-psi,Kf,mycolor[i])
        pl.xlabel('psi'); pl.ylabel('Kf'); pl.grid('on')
        pl.subplot(4,2,8)
        pl.loglog(-psi,Ke,mycolor[i],label=f'T={Ti}')
        pl.xlabel('psi'); pl.ylabel('Ke'); pl.grid('on')
    
    pl.subplot(4,2,8)
    pl.legend(ncol=10,loc='center left',bbox_to_anchor=(-0.75,-0.5))
