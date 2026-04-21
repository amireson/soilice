<img src="logo.png" width="300px">  

# About

$\text{soilice}$ v1.0 is a coupled mass and heat balance solver for frozen soils. The code is written in python and is designed to be concise, readable and easy to customize (try new constitutive relationships, etc) without having to re-compile the code. It will run on any platform. It uses a just-in-time compiler and an ODE solver, so is efficient and has excellent mass/energy conservation. User instructions are provided in this readme file below. 

The technical documentation is provided [here](technicalDocumentation/soiliceDoc.pdf).

Some demonstration simulations are provided [here](notebooks/readme.md).

# Installation

$\text{soilice}$ is installed using the [pip](https://pypi.org/project/pip/) python package installer.

Optionally you might want to create a new python virtual environment before installing this (see [here](https://packaging.python.org/en/latest/guides/installing-using-pip-and-virtual-environments/)).

To install directly from github, enter this command in your terminal/PowerShell:
```
pip install git+ssh://git@github.com/amireson/soilice.git
```

Alternatively, obtain the source code by cloning this repo, then navigate to the root folder of the repo and install using the command `pip install .`. 

After installation you can import the model to use it from anywhere on your computer. Note, if you are interested in editing the source code, read the section below on customizing $\text{soilice}$.

# Quick start 

The following code can be used for a very simple model setup and run for a soil profile subject to heat conduction from fixed temperature boundaries at the top and bottom.

```python
from soilice import model
import numpy as np
import matplotlib.pyplot as pl

sim = model()

# Options (simulation with no flow, conduction only)
sim.opts = {
 'simulateFlow': 0.0,
 'gravity': 1.0,
 'infiltration': 1.0,
 'cryoK': 1.0,
 'cryoGradient': 0.0,
 'freeDrainage': 1.0,
 'simulateTransport': 1.0,
 'withadv': 0.0,
 'conductionTop': 1.0,
 'conductionBot': 1.0
}

# Grid
sim.zGrid(np.arange(0,1.01,0.01))
sim.tGrid(0,10,0.1)

# Parameters
from soilice import writeDefaultPars
writeDefaultPars()
sim.readPars()

# Initial and boundary conditions
sim.setICs(T0=-1, psi0=0.)
sim.setBCs(TTop=1.,TBot=-1)

# Run
out = sim.run()

pl.plot(out.T.T,out.z,'b')
pl.ylim(1,0); pl.xlabel('Temperature'); pl.ylabel('Depth')
pl.grid(); pl.show()
```

# Running the testcases

For a more complete introduction of how to setup and run the model, and visualize the output, refer to the example notebooks. These are simple boundary value problems that are configured in Jupyter notebooks, available [here](notebooks/runSoil_Infiltration.ipynb).

# User guide

The $\text{soilice}$ model simulates coupled flow of liquid water and heat transport in a variable saturated and variably frozen soil profile. This can be configured to solve flow only, transport only or coupled flow and transport. There is also flexibility in the boundary conditions that are used, as described below. Below is a step by step guide to setting up and running the model. This is easiest to do within a Jupyter notebook, at least for a new model setup and run - here each step of the guide is to be carried out within one cell of a Jupyter notebook. Note, the steps don't necessarily have to be carried out in this order.

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

Optionally you may change the `rtol` parameter of the ODE solver when setting up the model, using the syntax

```python
sim=model(rtol=1e-8)
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

## Step 4: The model time grid

The model uses a regular time grid, with time step `dt`, start time `t0`, and end time `tMax`. Note, that time units must be consistent with the hydraulic and thermal conductivity parameter values. The time grid is assigned to $\text{soilice}$ with the command

```python
sim.tGrid(t0,tMax,dt)
```

Now the time grid is saved in the variable `sim.t` and the number of time steps is `sim.nt`.

## Step 5: Assigning parameters and constants


$\text{soilice}$ requires the user to define a number of parameters (variables that might change from soil to soil) and constants (variables that are unlikely to change from soil to soil) in python dictionaries. It is essential that every parameter that is used in the source code (see Customizing a function section below) is defined. If using the default constitutive functions, it is recommended that the user imports default parameter values, and edits them as needed. To import the default parameters and constants into a text file in the working folder run the following commands:

```python
from soilice import writeDefaultPars
writeDefaultPars(filename='myPars')
```

The filename here is optional. This will create two text files containing default constants and parameters for a guelph loam defined originally by van Genuchten (1980). These parameters are not yet associated with the model - to do this execute the command:

```python
sim.readPars('myPars')
```

This will read the parameters from the file `myPars_pars.txt` into the library `sim.pars` and the constants from the file `myPars_const.txt` into the library `sim.const`. 

Parameters can also be assigned different values for different layers in the following manner. For a **uniform profile**, every parameter must be assigned a scalar value (e.g. `thetaS=0.4`). If any single parameter has more than one value assigned, then all parameters are defined for every layer individually, and the model uses a **non-uniform** or **layered profile**. In a layered profile the following options are available and each parameter is handled individually:

number of values defined | How this is treated | Example
--- | --- | ---
1 (scalar or array)  | The parameter is duplicated for each layer - i.e. this parameter is uniform | `thetaS=0.4` or `thetaS=[0.4]`
2 (array) | The parameter is scaled linearly from top to bottom | `thetaS=[0.4,0.2]`
3 (array) | The parameter is scaled exponentially, with the 3 values representing the value at `z=0`, `z=1.0`, and `z=infinity` | `thetaS=[0.4,0.3,0.2]`
`nz` (array) | The parameter is explicitly defined for each individual layer | `thetaS=[0.4,0.4,0.3,0.3,... ]`

The file `myPars_pars.txt` can be edited to the formats shown in the examples above - that is to say the model can parse lists in square brackets. After editing the parameter file, read the parameters into the model with the command `sim.readPars('myPars_pars.txt')`.

Finally, it is also possible to define the parameters by creating a dictionary with their values (scalars or arrays), say called `pars`, and then assigning this dictionary to the model by simply running `sim.pars=pars`

## Step 6: The model initial conditions

Two initial conditions are always required for the model - the initial temperature, `T0`, and the initial effective matric potential (that is the matric potential associated with the total water content), `psi0`. These are assigned with the command:

```python
sim.setICs(T0=T0,psi0=psi0)
```

Each must be defined for every node. There are three ways to define these:

1. Provide a scalar value, e.g. `T0=0.`, which will assign uniform initial conditions;

2. Provide two values in an array, e.g. `psi0=[-1,0]`, which will assign values scaled linearly with `z`, with the first and second values representing the value at the top and bottom of the domain, respectively;

3. Provide `nz` values in an array, e.g. `T0=np.array([10,9,8,7,...])`, which will assign these values to each node.

As an example, for a 1 m depth profile, a hydrostatic initial condition, with a uniform temperature of -1 deg, would be defined by 

```python
sim.setICs(T0=-1,psi0=[-1,0])
```

## Step 7: The model boundary conditions

Model boundary conditions are assigned by a combination of setting the model options (see above) and assigning time series values. This table describes all the variables and when they are required

Variable | Meaning | Required when
--- | --- | ---
`qI` | Potential infiltration rate | `opts['infiltration']=1.0`
`psiT` | Matric potential on upper boundary | `opts['infiltration']=0.0`
`jTopBC` | Direct heat flux on upper boundary | Can always be used or set to zero.
`jBotBC` | Direct heat flux on lower boundary | Can always be used or set to zero.
`TInf` | Temperature of infiltrating water | `opts['withadv']=1.0`
`TTop` | Temperature of the upper boundary | `opts['conductionTop']=1.0`
`TBot` | Temperature of the lower boundary | `opts['conductionBot']=1.0`

Here we describe all the possible options.

### Upper flow boundary

There are two options for the upper flow boundary, set with the `opts['infiltration']` option, where `1.0` means a potential infiltration flux is applied (that generates infiltration and/or runoff, see the technical documentation for details), or `0.0` means a specified `psi` is applied. For the flux type boundary, the variable `qI` must be defined. For the type I boundary, the variable `psiT` must be defined.

### Lower flow boundary

There are two options for the lower flow boundary, set with the `opts['freeDrainage']` option, where `1.0` means a free draining lower boundary condition and `0.0` means a zero flux lower boundary condition. Neither option requires the user to define any variables.

### Upper heat transport boundary

On the upper boundary, there are three ways to add heat into the model, and they may all be acting simultaneously, with the net heat flux given by the sum of all of them.  The boundary fluxes are as follows:

1. **A directly applied heat flux** that is independent of any movement of water, for example representing the effect of net radiation. To assign this the user defines the variable `jTopBC` (J m$^{-2}$ d$^{-1}$). If `jTopBC` is not defined, it will take on a zero value, which will switch off this particular boundary flux.

2. **An advective heat flux** is defined as a function of the infiltration rate, defined using the upper flow boundary, and the temperature of the infiltrating water, defined with the variable `TInf`. Note that if `TInf` is not defined, it will take on a zero value, which will switch off this heat flux (but note this does not change the mass flux over the boundary).

3. **A conductive heat flux** which depends on the thermal conductivity and the temperature gradient. To turn this on it is necessary to set `opts['conductionTop']=1.0`, and now the gradient is calculated based on the variable `TTop`. If `TTop` is not defined it would take on a zero value, but that would still be used to define the conductive heat flux, likely in error. Therefore, do ensure that `opts['conductionTop']=0.0` if you do not wish for there to be a conductive flux.

### Lower heat transport boundary

The lower heat boundary is essentially the same as the upper boundary - with three possible boundary fluxes. A direct heat flux can be applied, for example to represent the geothermal heat flux, using the `jBotBC` variable. Again, assigning a zero value to `jBotBC` turns that flux off, so no need to set any options. There will be an advective heat flux if there is free drainage occurring, and that will depend on the temperature of the lower grid cell, so the user does not set anything for this. In addition, there can be a conductive heat flux, and this is the same as the upper boundary: to turn it on set `opts['conductionBot']=1.0` and assign the temperature to `TBot`. To turn it off set `opts['conductionBot']=0.0`.

### Assigning the boundary conditions

For all the boundary condition variables, the variable must be defined either as a scalar value (i.e. constant in time) or as a time series that matches the time grid of the model (i.e. has `nt` points). If no value of any given variable is provided, a zero value is assigned. So, to define all boundary conditions you would use

```python
sim.setBCs(jTopBC=jTopBC, qI=qI, psiT=psiT, TInf=TInf, TTop=TTop, TBot=TBot)
```

But in practice one or more are likely to be zeros (e.g. you don't need to set `qI` and `psiT` as only one would be used). 

As an example, consider we are running a simulation with no flow (and hence no advection) and only heat conduction over the boundaries from constant temperature of +1 deg C at the ground surface and -1 deg C at the base. Then we would run

```python 
sim.setBCs(TTop=1.,TBot=-1.)
```

It would also be essential to ensure the options are defined appropriately.

Note that it is always necessary to run the command `sim.setBCs()` even if all boundary condition variables are zeros. 

## Step 8: Running the model

To run the model, all above steps must be completed, and then enter:

``` python
out=sim.run()
```

The variable `out` is an object that stores all the model output, including the following:

Variable | Meaning | Dimensions
--- | --- | ---
`out.t` | Time grid | `(nt)`
`out.z` | Space grid | `(nz)`
`out.thetaL` | Liquid water content | `(nt,nz)`
`out.thetaT` | Total water content | `(nt,nz)`
`out.thetaI` | Ice water content | `(nt,nz)`
`out.psie` | Effective matric potential | `(nt,nz)`
`out.psif` | Matric potential associated with liquid water | `(nt,nz)`
`out.T` | Temperature | `(nt,nz)`
`out.qT` | Actual infiltration flux | (`nt`)
`out.qB` | Actual drainage flux | (`nt`)
`out.jT` | Total upper heat flux | (`nt`)
`out.jB` | Total lower heat flux | (`nt`)

## Step 9: Checking mass and energy balance

After running the model you can report the mass and energy balance closure to the screen with the command:

```python
out.balanceClosure()
```

To plot the mass and energy balance for the soil profile run the command:

``` python
out.plotBalance()
```
## Step 10: Checking the constitutive functions

It is important to understand that the constutive functions in $\text{soilice}$ depend on both `psie` and `T`. As such it is a good idea to always plot your relationships carefully against `psie` and `T` to ensure everything is as you expect. Note that the impedance model for hydraulic conductivity that is fairly widely used is not monotonic under certain parameter combinations. The notebook developed [here](notebooks/00_Properties.ipynb) can be used for this purpose.

## Customizing a function

$\text{soilice}$ is designed so that you can make a local copy of parts of the source code in a working folder and edit the functions. The code is organized into three python scripts:

`src_soil.py` is used to configure and run the model. It is unlikely users would ever need to edit this script, and there is no simple way to do that.

`src_constitutiveFunctions.py` is used to define all the constitutive relations - that is the soil hydraulic and thermal properties. Review the default configuration [here](soilice/src_constitutiveFunctions.py). To edit these functions, the user must copy this script into their working folder. Whenever $\text{soilice}$ is run from within that folder the local copy of this script is used. Run the following commands (e.g. in a Jupyter notebook) to make a local copy of the script:

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
`src_conservationFunctions.py` is used to define the mass and energy conservation equations [see here](soilice/src_conservationFunctions.py). This includes two main functions - a mass balance function that is essentially solving Richards' Equation; and an energy balance function that is essentially solving the advection-diffusion equation, in both cases with modifications to account for phase change/latent heat, as fully described in the  [technical documentation](technicalDocumentation/soiliceDoc.pdf). These functions can also be copied locally and edited in the same way as the constitutive functions, though it is less likely that a user would want to do that. To do this, run the following script:

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
