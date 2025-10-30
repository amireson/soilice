---

<img src=logo.png width=300px>  

---

# Model overview

The $\text{soilice}$ model simulates couples flow of liquid water and heat transport in snow and soil. The model has two particular unique features:

1. The physically based soil mass and heat transport model is fully coupled with a physically based snow model, such that the soil upper boundary condition includes radiation (heat from the atmosphere), conduction (heat from the snowpack) and advection (heat from infiltration, and turbulent heat fluxes)
2. The model is written in a python, such that the code is concise, highly readable and easy to customize (try new constituitive relationships, etc) without having to re-compile the code. It will run on any platform. It uses a just-in-time compiler and an ODE solver, so is efficient and has excellent mass/energy conservation.

The model can be configured in the following basic way:

### Snowpack model
Simulates a standalone snowpack, driven by hydrometric observations. The model simulates the snow mass, depth, heat content, sublimation and snowmelt. To run this configuration you need the import statement `from src_snow import run as run`

### Coupled snow-soil model
Simulates the snowpack as above, but with heat and mass exchange with the underlying soil profile. The soil upper boundary is a mass and heat flux (type II boundary). To run this configuration you need the import statement `from src_coupled import run as run`

### Soil model
Simulates a standalone soil profile. This can be configured to solve flow only, transport only or coupled flow and transport. There is also flexibility in the boundary conditions that are used, as described below. To run this configuration you need the import statement `from src_soil import run as run`

# Configuring the soil model

To configure the soil model you must set the following options, in the dictionary `opts`.

`opts['massflag']` is set to `0.` to make the assumption $\frac{dmc_pT}{dt}=mc_p\frac{dT}{dt}$ 

`opts['massflag']` is set to `1.` to make the assumption $\frac{dmc_pT}{dt}=mc_p\frac{dT}{dt}+Tc_p\frac{dm}{dt}$ (which is more correct). This allows the user to explore the impact of the simplifying assumption.

`opts['gravity']` is set to `0.` for horizontal flow (no gravity component) or `1.` for vertical (positive downwards) flow. Note that in a vertical model configuration (with flow switched on) there is always a free-drainage boundary condition, while in horizontal configuration the lower boundary becomes a no-flow boundary.

`opts['freeDrainage']` is set to `0.` for a no flow lower boundary condition and `1.0` for a free draining lower boundary condition.

`opts['cryoflow']` determines whether the hydraulic gradient in the soil is set based on $\psi_e$ (the matric potential associated with the total water content, using `opts['cryoflow']=0.`) or $\psi_f$ (the matric potential associated with freezing, using `opts['cryoflow']=1.`).

`opts['withadv']` allows the user to turn on (`1.`) or off (`0.`) advection. When off heat transport is by conduction only.

`opts['conductionTop']` is set to 1. to calculate the upper boundary heat flux based on the specified temperature, `TTop`. This conductive heat flux is added to the specified advective heat flux `jTop`. Hence to use a type I upper boundary for heat transport, set `jTop=np.zeros(nt)`, set `opts['conductionTop']=1.` and set `TTop=f(t)`. To use a type II upper boundary for heat transport, set `opts['conduction']=0.` and set `jTop=f(t)` (while the values in `TTop` are not used by the model). A mixed boundary is also possible where both advective and conductive fluxes are set.

`opts['conductionBot']` is the same as the previous option, but for the lower heat transport boundary condition. In this case, the advective flux only comes in if there is free drainage happening, and can only be an outflow.

`opts['simulateFlow']` is set to `True` or `False` to determine whether or not the model simulates flow. If `False` there is no change in total water content, no change in $\psi_e$ and zero mass fluxes $q$, which means the initial condition dictates the distribution of total water content in the profile. 

`opts['simulateTransport']` is set to `True` or `False` to determine whether or not the model simulates flow. If `False` there is change in internal energy or temperature and the heat fluxes $j$ are zero, which means the initial condition dictates the soil temperature and hence partitioning of total water content into ice and liquid. 

Turning off flow and transport allows the user to quickly see the equilibrium distribution of liquid water and ice in a soil profile for a given initial (steady-state) condition. 

# Parameters and constants

In the model `pars` is a dictionary that defines all model parameters - meaning constants that might change under different soil conditions. Each parameter can either by given a scalar value (e.g. `pars['thetaS']=0.4`), for a uniform profile, or it must have a unique value defined for very soil layer (e.g. `pars['thetaS'][:5]=0.4; pars['thetaS'][5:10]=0.3`, where `nz=10`) for a layered profile. Before sending the dictionary into the model it must be converted to a `numba` compatible dictionary with either `MakeDictFloat` (for uniform profile) or `MakeDictArray` (for a layered profile). 

Note that the parameters can be imported from the `HP_xxx.py` and `thermalSoilParameters.py` files.

The dictionary `const` includes all the model parameters that can be considered constants and will not change with depth or between different model runs. Examples include the density of water, `const['rho_liq']`. The constants can be imported from the `constants.py` file.

# The model space grid

The user must specify the following variables to define the model space grid, noting the `z` represents depth below ground in (m):

`nz` an integer representing the number of soil layers <br>
`zMax` the total depth of the soil profile <br>
`dz` an array of dimension `nz` that represents the depth of each soil layer (starting at the top and going down), such that `np.sum(dz)==zMax`. <br>
`z` is not strictly needed by the model, but could be useful for plotting purposes. `z` is an array storing the midpoint depth of each cell. 

# The model time grid

The user must specify the following variables to define the model time grid, noting `t` is time in units that must be consistent with the hydraulic and thermal conductivity parameter values:

`t` an array storing the values of time where calculation outputs are to be stored. <br>
`dt` the time step for each calculation step. The first calculation performed would be stored at time `t[1]=t[0]+dt`. Note `dt` is currently setup to be constant, and this is probably a sensible decision to stick to. <br>
`nt` an integer representing the number of time steps.

# The model boundary conditions

The time varying boundary conditions must always be defined for every calculation step - even if they take a constant value. The user must define the following time varying boundary conditions:

`qI` the infiltration rate (can be zero) <br>
`jTop` the specified heat flux (can be zero or can be set `qI`) <br>
`TTop` the temperature of the upper boundary face (ususally `z=0`), always defined but only used if `conductionTop` is switched on. <br>
`TBot` the temperature of the lower boundary face (`z=zMax`), always defined but only used if `conductionBot` is switched on.

# The model initial conditions

The depth dependent initial conditions are set for `t=t[0]`. Each is an array with `nz` cells. The following variables must be set:

`T0` is the initial soil temperature <br>
`psi0` is the initial matric potential corresponding to the initial total water content.

# Running the model

To run the model, all above steps must be completed, and then enter:

`psie,T,thetaL,thetaI,qT,qB,jT,jB=
run(dt,t,dz,nz,T0,psi0,qI,TTop,TBot,jTop,parsD,const,opts,rtol=1e-8)`

The `rtol` parameter is optional (default value is `1e-7`) but can improve accuracy.

# Checking mass and energy balance

To plot the mass and energy balance for the soil profile, import the function

`from balanceChecks import soilBalanceCheck`

and then run

`soilBalanceCheck(t,thetaL,thetaI,T,qT,qB,jT,jB,dz,const,pars)`

# Checking the soil properties

To plot all the constitutive relationships that are hard coded in `src_soil.py`, import the function 

`from checkProperties import getProperties`

and then run

` psie,psif,thetaL,thetaI,thetaT,dthdT,CB,kappa,fdash,gdash,Kf,Ke=getProperties(psi,T,pars,const)`

Note, this function needs the parameters and constants to be defined and it needs values of `psi` and `T`.

