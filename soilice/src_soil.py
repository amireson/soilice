#############################################
#
# soilice - source code 
# Andrew Ireson 
# https://github.com/amireson/soilice
#
#############################################

# Import statements
import numpy as np

from numba import njit 
from numba import types
from numba.typed import Dict
from numba.experimental import jitclass

import time
from scipy.integrate import ode

import dill
import ast
from pathlib import Path

# First try to import constitutive functions from the local folder, where
# the user may put custom functions. If the customized file does not exist
# import the default functions
try:
    from src_constitutiveFunctions import ( 
        thetaFun, CFun, KFun, thermalKfun, SFCslope, GCEFun, CBFun )
except ImportError:
    from soilice.src_constitutiveFunctions import ( 
        thetaFun, CFun, KFun, thermalKfun, SFCslope, GCEFun, CBFun )
    
# As above if user specified edits to the conservation equations exist
# import these, otherwise import the defaults
try: 
    from src_conservationFunctions import Richards, heatbalanceFun, ODEfunCall
except ImportError:
    from soilice.src_conservationFunctions import Richards, heatbalanceFun, ODEfunCall

from .utils import modelInOut, writeDefaultPars

#############################################
#
# CONFIGURATION AND DATA HANDLING FUNCTIONS
# Do not edit
#
#############################################

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

# Model function wrapper, called by ODE solver
def ODEfun(t,DV,upperBC,TTop,TBot,TInf,jTopBC,jBotBC,dz,pars,const,opts,nz):
    return ODEfunCall(t,DV,upperBC,TTop,TBot,TInf,jTopBC,jBotBC,dz,pars,const,opts,nz)

# Main class that setups and runs the model
class model:
    def __init__(self,opts=None,rtol=1e-7):
        self.opts=opts
        self.rtol=rtol
            
    def setOpts(self):
        
        def AskQuestion(question,default,trueVal):
            answer=input(question).lower()
            answer=default if answer == '' else answer
            return 1.0 if answer == trueVal else 0.0

        opts={}
        # Simulate flow
        opts['simulateFlow']=AskQuestion(
                'Simulate flow, y/n? (default = y)',
                'y','y')

        if opts['simulateFlow']:
            # Horizontal/vertical
            opts['gravity']=AskQuestion(
                    'Horizontal or vertical, h/v? (default = v)',
                    'v','v')
            
            # Infiltration or fixed psi at upper boundary:
            opts['infiltration']=AskQuestion(
                    'Infiltration at upper boundary (n means fixed psi), y/n? (default = y)',
                    'y', 'y')
            
            # Cryo hydraulic conductivity (yes means K based on liquid water content, no means K based on total water content)
            opts['cryoK']=AskQuestion(
                    'Hydraulic conductivity based on liquid water content, y/n? (default = y)',
                    'y','y')

            # Cryosuction (yes means gradient based on psif, no means gradient based on psie)
            opts['cryoGradient']=AskQuestion(
                    'Use cryosuction gradient, y/n? (default = n)',
                    'n','y')

            # Simulate free Drainage
            opts['freeDrainage']=AskQuestion(
                    'Simulate free drainage on lower boundary (n means no flow), y/n? (default = y)',
                    'y','y')

        else:
            opts['gravity']=1.
            opts['infiltration']=1.
            opts['cryoK']=1.
            opts['cryoGradient']=0.
            opts['freeDrainage']=0.


        # Simulate heat transport
        opts['simulateTransport']=AskQuestion(
                'Simulate heat transport, y/n? (default = y)',
                'y','y')

        if opts['simulateTransport']:

            # Include advection
            opts['withadv']=AskQuestion(
                    'Include advection, y/n? (default = y)',
                    'y', 'y')

            # Conduction at the upper boundary based on TTonp
            opts['conductionTop']=AskQuestion(
                    'Conduction, based on TTop, at upper boundary, y/n? (default = n)',
                    'n','y')

            # Conduction at the lower boundary based on TBot
            opts['conductionBot']=AskQuestion(
                    'Conduction, based on TBop, at lower boundary, y/n? (default = n)',
                    'n','y')

        else:
            opts['withadv']=0.
            opts['conductionTop']=0.
            opts['conductionBot']=0.


        self.opts=opts

    def printOpts(self):
        opts=self.opts
        print('*************************************************************************\n     SUMMARY OF MODEL OPTIONS:')
        
        print(f"     * {'Simulating flow' if opts['simulateFlow']==1  else 'No flow'}")
        if opts['simulateFlow']:
            print(f"     * {'Vertical flow' if opts['gravity']==1  else 'Horizontal flow'}")
            print(f"     * {'Infiltration at ground surface' if opts['infiltration']==1  else 'Fixed psi at ground surface'}")
            print(f"     * {'Hydraulic conductivity based on liquid water content' if opts['cryoK']==1  else 'Hydraulic conductivity based on total water content'}")
            print(f"     * {'Uses cryosuction based gradient' if opts['cryoGradient']==1  else 'No cryosuction'}")
            print(f"     * {'Free draining lower boundary condition' if opts['freeDrainage']==1  else 'No (mass) flow lower boundary condition'}")

        print(f"     * {'Simulating heat transport' if opts['simulateTransport']==1  else 'No heat transport'}")
        if opts['simulateTransport']:
            print(f"     * {'Included heat advection' if opts['withadv']==1  else 'No advection of heat'}")
            print(f"     * {'Conduction on the upper boundary' if opts['conductionTop']==1  else 'No conduction on the upper boundary'}")
            print(f"     * {'Conduction on the lower boundary' if opts['conductionBot']==1  else 'No conduction on the lower boundary'}")
        print('*************************************************************************\n')

    def readPars(self,filename='def'):

        dest_path = Path('.') / f'{filename}_pars.txt'
        if dest_path.exists():
            pass
        else:
            writeDefaultPars(filename=filename)

        # Read parameters and constants from text file.
        filename=filename.replace('_pars','')
        filename=filename.replace('_const','')
        filename=filename.replace('.txt','')
        
        pars={}
        for line in open(f'{filename}_pars.txt','r'): 
            pars[line.split(',',1)[0].strip()]=ast.literal_eval(line.split(',',1)[1].strip())
        self.setPars(pars)

        const={}
        for line in open(f'{filename}_const.txt','r'): 
            const[line.split(',',1)[0].strip()]=ast.literal_eval(line.split(',',1)[1].strip())
        self.const=const
        
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
        
    def setBCs(self,jTopBC=0,jBotBC=0,qI=0,psiT=0,TInf=0,TTop=0,TBot=0):
        
        def _to_timeseries(x):
            if np.ndim(x) == 0:        
                return np.full(self.nt,x)
            else:
                return np.asarray(x)

        self.jTopBC = _to_timeseries(jTopBC)
        self.jBotBC = _to_timeseries(jBotBC)
        self.qI     = _to_timeseries(qI)
        self.psiT   = _to_timeseries(psiT)
        self.TInf   = _to_timeseries(TInf)
        self.TTop   = _to_timeseries(TTop)
        self.TBot   = _to_timeseries(TBot)

    def setICs(self,T0,psi0):

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
            elif len(x)==np.asarray(self.nz):
                # Fully specified 
                return x

        self.T0=_to_array(T0)
        self.psi0=_to_array(psi0)    

#############################################
#
# FUNCTION TO RUN MODEL
# Do not edit
#
#############################################
    def run(self):

        print("""
                                    iiii  lllllll   iiii
                                   i::::i l:::::l  i::::i
                                    iiii  l:::::l   iiii
                                          l:::::l
    ssssssssss      ooooooooooo   iiiiiii  l::::l iiiiiii     cccccccccccccccc    eeeeeeeeeeee
  ss::::::::::s   oo:::::::::::oo i:::::i  l::::l i:::::i   cc:::::::::::::::c  ee::::::::::::ee
ss:::::::::::::s o:::::::::::::::o i::::i  l::::l  i::::i  c:::::::::::::::::c e::::::eeeee:::::ee
s::::::ssss:::::so:::::ooooo:::::o i::::i  l::::l  i::::i c:::::::cccccc:::::ce::::::e     e:::::e
 s:::::s  ssssss o::::o     o::::o i::::i  l::::l  i::::i c::::::c     ccccccce:::::::eeeee::::::e
   s::::::s      o::::o     o::::o i::::i  l::::l  i::::i c:::::c             e:::::::::::::::::e
      s::::::s   o::::o     o::::o i::::i  l::::l  i::::i c:::::c             e::::::eeeeeeeeeee
ssssss   s:::::s o::::o     o::::o i::::i  l::::l  i::::i c::::::c     ccccccce:::::::e
s:::::ssss::::::so:::::ooooo:::::oi::::::il::::::li::::::ic:::::::cccccc:::::ce::::::::e
s::::::::::::::s o:::::::::::::::oi::::::il::::::li::::::i c:::::::::::::::::c e::::::::eeeeeeee
 s:::::::::::ss   oo:::::::::::oo i::::::il::::::li::::::i  cc:::::::::::::::c  ee:::::::::::::e
  sssssssssss       ooooooooooo   iiiiiiiilllllllliiiiiiii    cccccccccccccccc    eeeeeeeeeeeeee

... please wait, model is running...
                """)

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

        if self.opts['infiltration']:
            self.upperBC=self.qI
        else:
            self.upperBC=self.psiT
        
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
            params=(
                    self.upperBC[i],self.TTop[i],self.TBot[i],
                    self.TInf[i],self.jTopBC[i],self.jBotBC[i],
                    self.dz,parsD,constD,optsD,self.nz)
            
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
                  'jTopBC','jBotBC','TInf','TTop','TBot',
                  'T0','psi0','opts','rtol']: 
            setattr(inOut, i, getattr(self, i))

        if self.opts['infiltration']:
            inOut.qI=self.upperBC
        else:
            inOut.psiT=self.upperBC
        
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
        inOut.balanceClosure()
        print('*************************************************************************\n')
        return inOut



       

