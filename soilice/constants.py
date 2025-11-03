# Define model constants: (should not be considered free parameters)
const = {}
const['stefanBoltzmann']=5.670374419e-8 # J/m2/s/K
const['rho_liq']=1000.                  # Density of liquid water kg/m3
const['rho_ice']=918.                   # Density of ice kg/m3
const['rho_air']=1.293                  # Default air density kg/m3
const['cp_liq']=4180.                   # Specific heat capacity of liquid water J/kg/K
const['cp_ice']=2100.                   # Specific heat capacity of liquid ice J/kg/K
const['cp_air']=1006.                   # Specific heat capacity of air J/kg/K
const['kappa_liq']=0.56                 # Thermal conductivity of liquid water J/s/m/K
const['kappa_ice']=2.2                  # Thermal conductivity of ice J/s/m/K
const['kappa_air']=0.025                # Thermal conductivty of air J/s/m/K
const['lambda_f']=0.334e6               # Latent heat of fusion J/kg
const['lambda_v']=2.26e6                # Latent heat of vaporization J/kg
const['lambda_s']=2.835e6               # Latent heat of sublimation J/kg Latent heat of sublimation
const['g']=9.81                         # Gravity m/s2
const['T0']=273.15                      # Freezing temperature in K
const['Tf']=0.                          # Freezing temperature in deg C