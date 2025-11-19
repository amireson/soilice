# Geostudio function:
@jit(nopython=True)
def thermalKfun(psie,psif,T,pars,const):
    thetaL=thetaFun(psif,pars)
    thetaT=thetaFun(psie,pars)
    thetaI=const['rho_liq']/const['rho_ice']*(thetaT-thetaL)
    thetaG=pars['thetaS']-thetaT
    n = pars['thetaS']
    kunfroz=const['kappa_liq']**thetaT*pars['kappa_soil']**(1-n)*const['kappa_air']**thetaG
    kfroz=const['kappa_ice']**thetaT*pars['kappa_soil']**(1-n)*const['kappa_air']**thetaG
    kappa=thetaL/thetaT*kunfroz+thetaI/thetaT*kfroz
    return kappa

