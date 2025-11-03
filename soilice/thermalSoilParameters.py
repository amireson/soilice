# Thermal Parameters
pars={}
# Specific heat capacities
pars['cp_soil']= 850.         # J/kg/K
pars['cp_org']= 580.          # J/kg/K 

# Thermal conductivities
pars['kappa_soil']=2.9       # J/s/m/K
pars['kappa_org']=0.25    # J/s/m/K added

# Densities 
pars['rho_soil']=2600.       # kg/m3 changed 
pars['rho_org']=1300.     # kg/m3 added

# Misc
pars['a']=-0.5            # exponent in the thermal conductivity model (between -1 and 1 and not 0)