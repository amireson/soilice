<img src="logo.png" width="300px">  

# About

$\text{soilice}$ v1.0 is a coupled mass and heat balance solver for frozen soils. The code is written in python and is designed to be concise, readable and easy to customize (try new constituitive relationships, etc) without having to re-compile the code. It will run on any platform. It uses a just-in-time compiler and an ODE solver, so is efficient and has excellent mass/energy conservation. User instructions are provided in this readme file below. 

The technical documentation is provided [here](technicalDocumentation/soiliceDoc.pdf).

Some demonstration simulations are provided [here](notebooks/readme.md).

# Installation

$\text{soilice}$ is installled using the [pip](https://pypi.org/project/pip/) python package installer.

Optionally you might want to create a new python virtual environment before installing this (see [here](https://packaging.python.org/en/latest/guides/installing-using-pip-and-virtual-environments/)).

To install directly from github, enter this command in your terminal/power shell:
```
pip install git+ssh://git@github.com/amireson/soilice.git
```

Alternatively, obtain the source code by cloning this repo, then navigate to the root folder of the repo and install using the command `pip install .`. 

After installation you can import the model to use it from anywhere on your computer. Note, if you are interested in editing the source code, read the section below on customizing $\text{soilice}$.

# Running the testcases

The best way to get started with $\text{soilice}$ is to run the provided test cases. These are simple boundary value problems that are configured in jupyter notebooks, available [here](notebooks/runSoil_Infiltration.ipynb).

# User guide

The $\text{soilice}$ model simulates coupled flow of liquid water and heat transport in a variable saturated and variably frozen soil profile. This can be configured to solve flow only, transport only or coupled flow and transport. There is also flexibility in the boundary conditions that are used, as described below. Below is a step by step guide to setting up and running the model. This is easiest to do within a jupyter notebook, at least for a new model setup and run - here each step of the guide is to carried out within one cell of a jupyter notebook. 

## Step 1: Import and instantiate the model

$\text{soilice}$ is object oriented, so you start by creating a model object - that is you instantiate the model, and all model parameters and variables are saved in that model object. Here we will call our model `sim`. Execute the following commands:

```python
from soilice import model
sim=model()
```

Note, you will also want to import the numpy and pyplot libraries:

```python
import numpy as np
import matplotlib.pyplot as pl
```

## Step 2: Select model options

The model options can be set intuitively using the command:

```python
sim.setOpts()
```

The user is prompted for a number of options. Since you will not want to go through this step every time you run the model, you can output your options to the screen (after the above step), by typing in a cell the command

```python
sim.opts
```

Now you can copy-paste your options and store them in the options variable, using the command

```python
sim.opts={'simulateFlow': 1.0,
 'gravity': 1.0,
 'infiltration': 1.0,
 'cryoK': 1.0,
 'cryoGradient': 0.0,
 'freeDrainage': 1.0,
 'simulateTransport': 1.0,
 'withadv': 1.0,
 'conductionTop': 0.0,
 'conductionBot': 0.0}}
```

The above are the default options. This final command (with options changed to `0.0` or `1.0`) is all you need to define all the model options.

The meaning of each option is described as follows:

option | description
--- | ---
simulateFlow | 1.0/0.0 do/do not solve the mass balance equations
gravity | 0.0/1.0 a horizontal/vertical profile
infiltration | 0.0 use a fixed psi upper boundary condition; 1.0 use a potential infiltration flux upper boundary condition
cryoK | 0.0 hydraulic conductivity is a function of the total water content; 1.0 hydraulic conductivity is a function of the liquid water content
cryoGradient | 1.0/0.0 hydraulic gradient does/does not include cryosuction
freeDrainage | 0.0 zero mass flux flow boundary condition; 1.0 free draining lower boundary condition
simulateTransport | 1.0/0.0 do/do not solve the heat balance equations
withadv | 1.0/0.0 do/do not include the heat flux due to advection
conductionTop | 1.0/0.0 do/do not include a conductive heat flux on the upper boundary (based on the temperature variable TTop)
conductionBot | 1.0/0.0 do/do not include a conductive heat flux on the lower boundary (based on the temperature variable TBot)

Note that turning off flow and transport allows the user to quickly see the equilibrium distribution of liquid water and ice in a soil profile for a given initial (steady-state) condition. 

## Step 3: The model space grid

$\text{soilice}$ uses a 1D finite difference grid, with `z` positive in the downward direction (depth below ground, typically) and flexible space steps - that is, `dz` does not have to be constant. The user defines a numpy array describing the cell boundaries, `bz`, which will have `nz+1` values. The model calculates the midpoints which are used for the state variables, and for which there are `nz` values. The numpy array describing the boundaries is sent into the model with the command `sim.zGrid(bz)`. 

For example, if we have a grid from ground surface to a depth of 1 m below ground, with a 0.01 m space step, we can use this command:

```python
sim.zGrid(np.arange(0,1.01,0.01))
```

After doing this, the space grid is saved in the variable `sim.z` and the number of grid cells is `sim.nz`

## The model time grid

The model uses a regular time grid, with time step `dt`, start time `t0`, and end time `tMax`. Note, that time units must be consistent with the hydraulic and thermal conductivity parameter values. The time grid is assigned to $\text{soilice}$ with the command

```python
sim.tGrid(t0,tMax,dt)
```

Now the time grid is saved in the variable `sim.t` and the number of time steps is `sim.nt`.

## Parameters and constants


$\text{soilice}$ requires the user to define a number of parameters (variables that might change from soil to soil) and constants (variables that are unlikely to change from soil to soil) in python dictionaries. It is essentially that every parameter that is used in the source code (see Customizing a function section below). If using the default constitutive functions, it is recommended that the user imports default parameter values, and edits them as needed. To import the default parameters and constants into a text file in the working folder run the following commands:

```python
from soilice import writeDefaultPars
writeDefaultPars(filename='myPars')
```

The filename here is optional. This will create two text files containing default constants and parameters for a guelph loam defined originally by van Genuchten (1980). These parameters are not yet associated with the model - to do this execute the command:

```python
sim.readPars('myPars')
```

This will read the parameters from the file `myPars_pars.txt` into the library `sim.pars` and the constants from the file `myPars_const.txt` into the library `sim.const`. 

Parameters can also be assigned different values for different layers in the following manner. For a **uniform profile**, every parameter must be assigned a scalar value (e.g. `thetaS=0.4`)

If any single parameter has more than one value assigned, then all parameters are defined for every layer individually, and the model uses a **non-unform** or **layered profile**. In a layered profile the following options are available and each parameter is handled individually:

number of values defined | How this is treated | Example
--- | --- | ---
1 (scalar or array)  | The parameter is duplicated for each layer - i.e. this parameter is uniform | `thetaS=0.4` or `thetaS=[0.4]`
2 (array) | The parameter is scaled linearly from top to bottom | `thetaS=[0.4,0.2]`
3 (array) | The parameter is scaled exponentially, with the 3 values representing the value at `z=0`, `z=1.0`, and `z=infinity` | `thetaS=[0.4,0.3,0.2]`
`n` (array) | The parameter is explicitly defined for each individual layer | `thetaS=[0.4,0.4,0.3,0.3,... ]`

The file def_pars.txt can be edited to the formats shown in the examples above - that is to say the model can parse lists in square brackets. Then, after editing that file, the parameters are assigned in the same way, with the commands above.

## The model boundary conditions

The time varying boundary conditions must always be defined for every calculation step - even if they take a constant value. The user must define the following time varying boundary conditions:

`qI` the infiltration rate (can be zero) <br>
`jTop` the specified heat flux (can be zero or can be set `qI`) <br>
`TTop` the temperature of the upper boundary face (ususally `z=0`), always defined but only used if `conductionTop` is switched on. <br>
`TBot` the temperature of the lower boundary face (`z=zMax`), always defined but only used if `conductionBot` is switched on.

## The model initial conditions

The depth dependent initial conditions are set for `t=t[0]`. Each is an array with `nz` cells. The following variables must be set:

`T0` is the initial soil temperature <br>
`psi0` is the initial matric potential corresponding to the initial total water content.

## Running the model

To run the model, all above steps must be completed, and then enter:

``` python
psie,T,thetaL,thetaI,qT,qB,jT,jB=
run(dt,t,dz,nz,T0,psi0,qI,TTop,TBot,jTop,parsD,const,opts,rtol=1e-8)
```

The `rtol` parameter is optional (default value is `1e-7`) but can improve accuracy.

## Checking mass and energy balance

To plot the mass and energy balance for the soil profile, import and run the following:

``` python
from balanceChecks import soilBalanceCheck


soilBalanceCheck(t,thetaL,thetaI,T,qT,qB,jT,jB,dz,const,pars)
```

Note, this function must be run after running the model.

## Checking the soil properties

To plot all the constitutive relationships that are hard coded in `src_soil.py`, run

``` python
# Import soil parameters, or create manually:
from soilice.pars_loam import pars

# Import constants:
from soilice.constants import const

# Import plotting function:
from soilice.checkProperties import plotProperties

# Plot:
plotProperties(pars,const)
```

Alternatively to return all the variables for customized plots you may run:

``` python
from soilice.checkProperties import getProperties

psie,psif,thetaL,thetaI,thetaT,dthdT,CB,kappa,fdash,gdash,Kf,Ke=getProperties(psi,T,pars,const)
```

Note, this function needs the parameters and constants to be defined and it needs values of `psi` and `T`.

## Customizing a function

$\text{soilice}$ is designed so that you can make a local copy of parts of the source code in a working folder and edit the functions. The code is organized into three python scripts:

`src_soil.py` is used to configure and run the model. It is unlikely users would ever need to edit this script, and there is no simple way to do that.

`src_constitutiveFunctions.py` is used to define all the constitutive relations - that is the soil hydraulic and thermal properties. To edit these functions, the user must copy this script into their working folder. Whenever $\text{soilice}$ is run from within that folder the local copy of this script is used. Run the following commands (e.g. in a jupyter notebook) to make a local copy of the script:

```python
from soilice import copyConstitutiveFuns
copyConstitutiveFuns()
```

Now you have a local copy of `src_constitutiveFunctions.py` that you can easily edit. The only thing you cannot change is the function inputs/outputs. The specific functions, with their inputs and outputs that you can edit are:

```python
theta = thetaFun(psi,pars)
C     = CFun(psi,pars)
K     = KFun(psie,psif,pars,const)
kappa = thermalKfun(psie,psif,T,pars,const)
CB    = CBFun(psie,psif,pars,const)
dthdT = SFCslope(psie,psif,pars,const)
psi   = GCEFun(T,pars,const)
```
`src_conservationFunctions.py` is used to define the mass and energy conservation equations. This includes two main functions - a mass balance function that is essentially solving Richards' Equation; and an energy balance function that is essentially solving the advection-diffusion equation, in both cases with modifications to account for phase change/latent heat, as fully described in the  [technical documentation](technicalDocumentation/soiliceDoc.pdf). These functions can also be copied locally and edited in the same way as the constitutive functions, though it is less likely that a user would want to do that. To do this, run the following script:

```python
from soilice import copyConservationFuns
copyConservationFuns()
```

The specific functions, with their inputs and outputs that you can edit are:

```python
dthetaTdt,dpsiedt,q = Richards(
    t,psif,psie,dz,pars,const,
    opts,nz,upperBC)

dTdt,j = heatbalanceFun(
    t,psie,psif,T,TTop,TBot,jTopAdv,jTopNonAdv,
    dz,pars,const,opts,nz,dthetaTdt,q)
```
