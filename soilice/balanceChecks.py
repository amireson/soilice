import numpy as np
import matplotlib.pyplot as pl
import pandas as pd

def soilBalanceCheck(t,thetaL,thetaI,T,qT,qB,jT,jB,dz,const,pars):
    nt,nz=thetaL.shape
    ml=thetaL*dz*const['rho_liq']
    mi=thetaI*dz*const['rho_ice']
    ms=np.zeros((nt,nz))+((1-pars['thetaS'])*dz*pars['rho_soil'])
    u=(ml*const['cp_liq']+mi*const['cp_ice']+ms*pars['cp_soil'])*T-mi*const['lambda_f']
    ml=np.sum(ml,axis=1)
    mi=np.sum(mi,axis=1)
    ms=np.sum(ms,axis=1)
    u=np.sum(u,axis=1)
    m=ml+mi
    qT=qT.cumsum()
    qB=qB.cumsum()
    jT=jT.cumsum()
    jB=jB.cumsum()
    
    pl.figure()
    pl.subplot(2,1,1)
    pl.plot(t,(qT-qB)*const['rho_liq'],'.',label='net mass in')
    pl.plot(t,m-m[0],label='change in mass')
    pl.grid()
    pl.legend()
    
    pl.subplot(2,1,2)
    pl.plot(t,(jT-jB)/1e6,'.',label='net heat in')
    pl.plot(t,(u-u[0])/1e6,label='change in internal energy')
    pl.grid()
    pl.legend()
    pl.ylabel('Energy (MJ/m2)',fontsize=16)

def snowBalanceCheck(snow_df):
    # Balance plots 
    pl.figure(figsize=(10,5))
    
    pl.subplot(1,2,1)
    pl.title('Snow water balance')
    m=(snow_df['mi']+snow_df['ml']).values
    pl.plot(snow_df.index,m-m[0],'.',label='cum change in mass')
    pl.plot(snow_df.index,(snow_df['qR_snow']+snow_df['qS']-snow_df['qE_snow']-snow_df['qD']).cumsum(),label='net flux')
    pl.legend()
    pl.grid()
    pl.ylabel('mm')
    
    pl.subplot(1,2,2)
    pl.title('Snow energy balance')
    pl.plot(snow_df.index,(snow_df['u']-snow_df['u'].iloc[0])/1e6,'.',label='cum change in internal energy')
    pl.plot(snow_df.index,(snow_df['jRn_snow']+snow_df['jAPS']+snow_df['jAPR']-snow_df['jG']-snow_df['jH_snow']-snow_df['jE_snow']).cumsum()/1e6,label='net flux')
    pl.ylabel('MJ')
    pl.legend()
    pl.grid()
    
    pl.gcf().autofmt_xdate()