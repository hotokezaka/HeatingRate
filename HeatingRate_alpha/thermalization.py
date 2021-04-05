import numpy as np
from numpy import pi
from astropy import constants as const

c = const.c.cgs.value
eV = 1.6021766e-12
MeV = 1.0e6*eV
me = const.m_e.cgs.value
e = const.e.esu.value
mu = const.u.cgs.value
E = const.c*const.c*const.m_e
me_MeV = E.to('MeV').value
ma = 4.001506179127*mu
kappa_beta = 1.0 #MeV cm^2/g
kappa_alpha = -0.15
#alpha_max = 4.
#alpha_min = 1.0
#n = 4.
gam_t = 0.15#0.15
gam_t_alpha = 0.15#0.15

#E is kinetic energy in MeV
def calc_ad(E):
    p = np.sqrt(((E/me_MeV)+1.)*((E/me_MeV)+1.)-1.)
    return 1.+p*p/3./np.sqrt(p*p+1.)/(np.sqrt(p*p+1.)-1.)

#fitting function of stopping power time beta for Xe 0.01-10MeV
a_kappa = 2.1
def calc_kappa_beta(Etmp):
    E = Etmp*MeV
    Z = 54.
    A = 130.
    Imean =  482.*eV
    gamma = 1.0 + E/(me*c*c)
    tmp3 = 1.0 - 1.0/gamma/gamma
    beta = np.sqrt(tmp3)
    tmp1 = gamma*gamma*me*beta*beta*c*c*E/Imean/Imean/2.
    tmp2 = 1./gamma/gamma
    tmp = np.log(tmp1) - (2./gamma-1.+beta*beta)*np.log(2) + 1. - beta*beta + (1.-1./gamma)*(1.-1./gamma)/8.
    omegap = np.sqrt(4.*np.pi*1.0e5*e*e/me)
    tmp_f = np.log(1.123*me*c*c*c*beta*beta*beta/(e*e*omegap))
    kappabeta = beta*2.*np.pi*Z*np.power(e,4.)*tmp/(me*c*c*beta*beta)/MeV/(A*mu)\
                + beta*4.*np.pi*np.power(e,4.)*tmp_f/(me*c*c*beta*beta)/MeV/(A*mu)
    
    return kappabeta 

def calc_kappa_beta_alpha(Etmp):
#    x = E/0.8
#    tmp = (np.power(x,-0.215*a_kappa)+np.power(x,0.4*a_kappa))
#    ai_kappa = 1./a_kappa
    #fitting
    kappa_ion = 25./(1./np.power(Etmp/0.7,1.15) + 1./np.power(Etmp/0.7,-0.17))
    E = Etmp*MeV
    Z = 54.
    A = 130.
    chi = 2.
    Imean =  482.*eV
    gamma = 1.0 + E/(ma*c*c)
    tmp3 = 1.0 - 1.0/gamma/gamma
    beta = np.sqrt(tmp3)
    tmp1 = gamma*gamma*me*beta*beta*c*c*E/Imean/Imean/2.
    tmp2 = 1./gamma/gamma 
    #print tmp1
    omegap = np.sqrt(4.*np.pi*1.0e5*e*e/me)
    tmp = np.log(1.123*me*c*c*c*beta*beta*beta/(e*e*omegap))
    kappa_free =  chi*beta*4.*np.pi*np.power(e,4.)*tmp/(me*c*c*beta*beta)/MeV/(A*mu)
    return kappa_ion + kappa_free#0.79*np.power(tmp,ai_kappa)

#alpha = 0.0522589652371*1.0#n=5, vmax=0.3c, v0 = 0.1c #Sigma = alpha*M/r0^2
#alpha = 0.154157897023##n=9, vmax=0.4c, v0 = 0.1c 
alpha = 0.0895299932779 #n=5, vmax=0.4c, v0 = 0.2c

def calc_thermalization_time(E0, Mej, vej, Aave, alpha_max, alpha_min, n):
    vmax = alpha_max*vej
    vmin = alpha_min*vej
    kappa_beta = calc_kappa_beta(E0)
    rho_inv = (1.-np.power(alpha_max,-n+3.))/(n-3.)+0.5*(1.-alpha_min*alpha_min)
    #print 1./(4.*np.pi*rho_inv) 
    tmp = ((1.0-np.power(alpha_max,-2.*n+3.))/(2.*n-3.) + 1.-alpha_min)/(rho_inv*rho_inv)
    #####tmp2 is rho=rho0(v/vmin), v>vmin####
    tmp2 = (n-3.)*(n-3.)*(1.0-np.power(alpha_max,-2.*n+3.))/((1.-np.power(alpha_max,-n+3.))*(1.-np.power(alpha_max,-n+3.)))/(2.*n-3.)
    #print tmp2/(4.*np.pi), tmp/(4.*np.pi), tmp2/tmp
    gamma_ad = calc_ad(E0)
    te2 = tmp*kappa_beta*c*Mej/(4.*np.pi*E0*vej*vej*vej*3.*(gamma_ad-1.))

    return np.sqrt(te2)



def calc_thermalization_time_alpha(E0, Mej, vej, Aave, alpha_max, alpha_min, n):
    vmax = alpha_max*vej
    vmin = alpha_min*vej
    kappa_beta_alpha = calc_kappa_beta_alpha(E0)
    rho_inv = (1.-np.power(alpha_max,-n+3.))/(n-3.)+0.5*(1.-alpha_min*alpha_min)
    tmp = ((1.0-np.power(alpha_max,-2.*n+3.))/(2.*n-3.) + 1.-alpha_min)/(rho_inv*rho_inv)
    tmp2 = (n-3.)*(n-3.)*(1.0-np.power(alpha_max,-2.*n+3.))/((1.-np.power(alpha_max,-n+3.))*(1.-np.power(alpha_max,-n+3.)))/(2.*n-3.)
    gamma_ad = 5./3.
    te2 = tmp*kappa_beta_alpha*c*Mej/(4.*np.pi*E0*vej*vej*vej*3.*(gamma_ad-1.))
    return np.sqrt(te2)


def calc_density(Mej,vej,t,alpha_max,alpha_min,n):
    r0 = vej*t
    rho_inv = (1.-np.power(alpha_max,-n+3.))/(n-3.)+0.5*(1.-alpha_min*alpha_min)
    tmp = ((1.0-np.power(alpha_max,-2.*n+3.))/(2.*n-3.) + 1.-alpha_min)/(rho_inv*rho_inv)
    
    return tmp*Mej/(4.0*np.pi*np.power(r0,3.))

def calc_zero_energy(tau0, tau1, E0):
    gamma_ad = calc_ad(E0)
    x_ad = 3.*(gamma_ad-1.)
    tmp = (1.0+gam_t)/(2.-(gam_t+1.)*x_ad)
    tmpp = -2.0 + (1.+gam_t)*x_ad
    x = tau1/tau0
    return 1.0 + tmp*(np.power(x,tmpp)-1.0)/(tau0*tau0)

def calc_zero_energy_alpha(tau0, tau1):
    x_ad = 2.
    tmp = (1.0+gam_t_alpha)/(2.-(gam_t_alpha+1.)*x_ad)
    tmpp = -2.0 + (1.+gam_t_alpha)*x_ad
    x = tau1/tau0
    return 1.0 + tmp*(np.power(x,tmpp)-1.0)/(tau0*tau0)   

#atau1 = x_ad
#ntau = -gam_t
#atau1 = x_ad

def epsilon_tau(tau0,tau,E0):
    x = tau/tau0
    gamma_ad = calc_ad(E0)
    atau1 = 3.*(gamma_ad - 1.)
    ntau = -gam_t
    pp = -2. + atau1*(1.-ntau)
    ppp = 1./(1.-ntau)
    tmp = np.power(tau0,-atau1*(1.-ntau))*(1.-ntau)/(atau1*(-ntau+1.)-2.)*(np.power(tau,pp)-np.power(tau0,pp))
    tmpp =  np.power(x,-atau1*(1.-ntau))*(1.-tmp)

    if(tmpp > 0.):
        e_tmp =  np.power(tmpp,ppp)
        if(e_tmp > 0.01):
            return np.power(e_tmp,ntau) #note gam_t is chosen to be constant
        else:
            return 0.
    else:
        return 0.

def epsilon_tau_alpha(tau0,tau,E0):
    x = tau/tau0
    gamma_ad = 5./3.
    atau1 = 3.*(gamma_ad-1.)
    ntau = -gam_t_alpha
    pp = -2.+atau1*(1.-ntau)
    ppp = 1./(1.-ntau)
    tmp = np.power(tau0,-atau1*(1.-ntau))*(1.-ntau)/(atau1*(-ntau+1.)-2.)*(np.power(tau,pp)-np.power(tau0,pp))
    tmpp =  np.power(x,-atau1*(1.-ntau))*(1. - tmp)

    if(tmpp > 0.):
        e_tmp =  np.power(tmpp,ppp)
        if(e_tmp > 0.01):
            return np.power(e_tmp,ntau) #note gam_t is chosen to be constant
        else:
            return 0.
    else:
        return 0.


def calc_gamma_deposition(kappa_eff, t, Mej, vej, alpha_min, alpha_max, n):
    w = alpha_min/alpha_max
    k = n - 3.0
    alpha_gam = 0.1*w + 0.003*k/w
    t0 = np.sqrt(alpha_gam*kappa_eff*Mej/(vej*vej))
   # print t0/day
    return 1.0 - np.exp(-t0*t0/(t*t))
