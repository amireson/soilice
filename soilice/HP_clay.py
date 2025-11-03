# Clay properties
pars={}
pars['thetaR']=0.10
pars['thetaS']=0.5 #0.6
pars['alpha']=0.152  # van Genuchten (1980) m-1
pars['n']=1.15 #1.17       # van Genuchten (1980)  
pars['eta']=0.5
pars['Ks']=9.4e-16 #-9   # m/s (from van Genuchten 0.082 cm/day)
pars['neta']=0.5
pars['Ss']=1e-6
pars['m']=1-(1/pars['n']) 
pars['theta_org']=0.0
pars['theta_mineral'] = (1-pars['thetaS'])