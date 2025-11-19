def KFun(psie,psif,pars,const):
    Se=(1+(psie*-pars['alpha'])**pars['n'])**(-pars['m'])
    Se[psie>0.]=1.0
    Ke=pars['Ks']*Se**pars['neta']*(1-(1-Se**(1/pars['m']))**pars['m'])**2
    
    Se=(1+(psif*-pars['alpha'])**pars['n'])**(-pars['m'])
    Se[psif>0.]=1.0
    Kf=pars['Ks']*Se**pars['neta']*(1-(1-Se**(1/pars['m']))**pars['m'])**2

    K=(Ke+Kf)/2.
    return K


