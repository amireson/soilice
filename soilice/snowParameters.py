# Snow parameters (can change):
pars={}
pars['Tm'] = 0.
pars['rho_newSnow']=100.              # kg/m3
pars['rho_maxSnow']=300.              # kg/m3
pars['snow_aging_coefficient']=0.01/3600 # *86400
pars['albedoSnow']=0.8
pars['albedoSoil']=0.3
pars['emissivitySnow']=1.00
pars['emissivitySoil']=1.00
pars['z0s']=0.01        # m
pars['z0g']=0.1         # m

pars['drainage_coeff']=0.00002
pars['zeroMass']=0.001  # mass (kg) considered negligible - snow less than this infiltrates into the ground.

pars['Kl'] = 0.56   # J/s/m/K water
pars['Ki'] = 2.24   # J/s/m/K ice
pars['Ka'] = 0.025  # J/s/m/K air
pars['b']  = 2  # Essery const
