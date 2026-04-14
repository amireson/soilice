import dill
import numpy as np
import matplotlib.pyplot as pl
from pathlib import Path
import shutil
import importlib.resources as resources


def save(filename,output):
    from .src_soil import modelInOut
    data=modelInOut()
    for k,v in vars(output).items():
        setattr(data, k, v)
    with open(filename, "wb") as f:
        dill.dump(data, f)

def load(filename):
    with open(filename, "rb") as f:
        return dill.load(f)

def loadModel(filename):
    from .src_soil import model
    sim=model(opts=None)
    for k,v in vars(load(filename)).items():
        setattr(sim, k, v)
    return sim

def writeDefaultPars(daily=True,filename='def'):
    
    # Check if pars file already exists
    dest_path = Path('.') / f'{filename}_pars.txt'
    if dest_path.exists():
        print(f"{filename}_pars.txt already exists in local folder")
        return

    # Save default parameter and constant values to textfiles in the working folder
    
    dailyScaling=86400 if daily else 1
    # hydraulic properties based on VG 1980 Guelph Loam (Drying)
    pars={}
    pars['thetaR']=0.218
    pars['thetaS']=0.520
    pars['alpha']=1.15
    pars['n']=2.03
    pars['eta']=0.5
    pars['Ks']=0.316/86400*dailyScaling # m/s
    pars['neta']=0.5
    pars['Ss']=1e-6
    pars['m']=1-(1/pars['n']) 
    pars['theta_org']=0.0
    pars['theta_mineral'] = (1-pars['thetaS'])
    
    # Thermal Parameters
    # Specific heat capacities
    pars['cp_soil']= 850.                   # J/kg/K
    pars['cp_org']= 580.                    # J/kg/K 
    
    # Thermal conductivity
    pars['kappa_soil']=2.9*dailyScaling     # J/s/m/K
    pars['kappa_org']=0.25*dailyScaling     # J/s/m/K added
    
    # Densities 
    pars['rho_soil']=2600.        # kg/m3 changed 
    pars['rho_org']=1300.         # kg/m3 added

    # Misc
    pars['q']=0.                            # Default flow rate when not solving RE
    pars['impedance']=0.                    # Impedance for frozen soil K (not recommended)

    const = {}
    const['stefanBoltzmann']=5.670374419e-8 # J/m2/s/K
    const['rho_liq']=1000.                  # Density of liquid water kg/m3
    const['rho_ice']=918.                   # Density of ice kg/m3
    const['rho_air']=1.293                  # Default air density kg/m3
    const['cp_liq']=4180.                   # Specific heat capacity of liquid water J/kg/K
    const['cp_ice']=2100.                   # Specific heat capacity of liquid ice J/kg/K
    const['cp_air']=1006.                   # Specific heat capacity of air J/kg/K
    const['kappa_liq']=0.56*dailyScaling    # Thermal conductivity of liquid water J/s/m/K
    const['kappa_ice']=2.2*dailyScaling     # Thermal conductivity of ice J/s/m/K
    const['kappa_air']=0.025*dailyScaling   # Thermal conductivty of air J/s/m/K
    const['lambda_f']=0.334e6               # Latent heat of fusion J/kg
    const['lambda_v']=2.26e6                # Latent heat of vaporization J/kg
    const['lambda_s']=2.835e6               # Latent heat of sublimation J/kg Latent heat of sublimation
    const['g']=9.81                         # Gravity m/s2
    const['T0']=273.15                      # Freezing temperature in K
    const['Tf']=0.                          # Freezing temperature in deg C

    f=open(f'{filename}_pars.txt','w')
    for k in pars: f.write(f'{k:>18}, {pars[k]}\n')
    f.close()

    f=open(f'{filename}_const.txt','w')
    for k in const: f.write(f'{k:>18}, {const[k]}\n')
    f.close()

# class for saving model output
class modelInOut:
    def __init__(self):
        pass

    # Balance closure functions:
    def balanceClosure(self):

        def err(qT,qB,dm):
            rmse=np.sqrt(np.mean((qT-qB-dm)**2))
            bias=(qT[-1]-qT[0])-(qB[-1]-qB[0])-dm[-1]
            return rmse,bias
        
        nt,nz=self.thetaL.shape
        ml=self.thetaL*self.dz*self.const['rho_liq']
        mi=self.thetaI*self.dz*self.const['rho_ice']
        ms=np.zeros((nt,nz))+((1-self.pars['thetaS'])*self.dz*self.pars['rho_soil'])
        u=(ml*self.const['cp_liq']+mi*self.const['cp_ice']+ms*self.pars['cp_soil'])*self.T-mi*self.const['lambda_f']
        ml=np.sum(ml,axis=1)
        mi=np.sum(mi,axis=1)
        u=np.sum(u,axis=1)
        du=u-u[0]
        m=ml+mi
        dm=m-m[0]
        qT=self.qT.cumsum()*self.const['rho_liq']
        qB=self.qB.cumsum()*self.const['rho_liq']
        jT=self.jT.cumsum()
        jB=self.jB.cumsum()
        rmseW,biasW=err(qT,qB,dm)
        rmseH,biasH=err(jT,jB,du)
        
        if self.opts['simulateFlow']:
            print(f'     Mass balance rmse: {rmseW: .2e} kg')
            print(f'                  bias: {biasW: .2e} kg')
        if self.opts['simulateTransport']:
            print(f'   Energy balance rmse: {rmseH: .2e} J')
            print(f'                  bias: {biasH: .2e} J')
        
        self.u=u
        self.m=m


    def plotBalance(self):
        
        t=self.t
        nt,nz=self.thetaL.shape
        ml=self.thetaL*self.dz*self.const['rho_liq']
        mi=self.thetaI*self.dz*self.const['rho_ice']
        ms=np.zeros((nt,nz))+((1-self.pars['thetaS'])*self.dz*self.pars['rho_soil'])
        u=(ml*self.const['cp_liq']+mi*self.const['cp_ice']+ms*self.pars['cp_soil'])*self.T-mi*self.const['lambda_f']
        ml=np.sum(ml,axis=1)
        mi=np.sum(mi,axis=1)
        u=np.sum(u,axis=1)
        du=u-u[0]
        m=ml+mi
        dm=m-m[0]
        qT=self.qT.cumsum()*self.const['rho_liq']
        qB=self.qB.cumsum()*self.const['rho_liq']
        jT=self.jT.cumsum()
        jB=self.jB.cumsum()
        pl.subplot(2,1,1)
        pl.plot(t,qT-qB,'.',label='Net water balance flux')
        pl.plot(t,dm,'-',label='Cumulative change in mass')
        pl.grid(); pl.legend()
        pl.subplot(2,1,2)
        pl.plot(t,jT-jB,'.',label='Net heat balance flux')
        pl.plot(t,du,'-',label='Cumulative change in internal energy')
        pl.grid(); pl.legend()


def copyConstitutiveFuns(dest_dir="."):
    filename = "src_constitutiveFunctions.py"
    dest_path = Path(dest_dir) / filename

    # 1. Check if file already exists
    if dest_path.exists():
        print(f"{filename} already exists at {dest_path}")
        return

    # 2. Access file inside package
    try:
        with resources.files("soilice").joinpath(filename).open("rb") as src:
            with open(dest_path, "wb") as dst:
                shutil.copyfileobj(src, dst)
    except FileNotFoundError:
        raise FileNotFoundError(f"{filename} not found in package")

    print(f"Copied {filename} to {dest_path}")

