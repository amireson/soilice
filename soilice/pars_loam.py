# Loam properties:

# hydraulic properties based on VG 1980 Guelph Loam (Drying)
pars={}
pars['thetaR']=0.218
pars['thetaS']=0.520
pars['alpha']=1.15
pars['n']=2.03
pars['eta']=0.5
pars['Ks']=0.316/86400 # m/s
pars['neta']=0.5
pars['Ss']=1e-6
pars['m']=1-(1/pars['n']) 
pars['theta_org']=0.0
pars['theta_mineral'] = (1-pars['thetaS'])

# Thermal Parameters
# Specific heat capacities
pars['cp_soil']= 850.         # J/kg/K
pars['cp_org']= 580.          # J/kg/K 

# Thermal conductivity
pars['kappa_soil']=2.9        # J/s/m/K
pars['kappa_org']=0.25        # J/s/m/K added

# Densities 
pars['rho_soil']=2600.        # kg/m3 changed 
pars['rho_org']=1300.         # kg/m3 added

# Misc
pars['a']=-0.5                # exponent in the thermal conductivity model (between -1 and 1 and not 0)
pars['q']=0.                  # default value for flow rate when not solving RE
