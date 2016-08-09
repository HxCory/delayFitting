import os
import cmath
import scipy.special as spcl
import scipy.optimize as sopt
import numpy as np
import matplotlib.pylab as plb
import matplotlib.pyplot as plt
import pylab

omega = np.arange(0.5, 1.2, 0.005)
m = range(1, 3)
E_0 = -1.1591
E_m = [-0.8502, -0.3278, -0.1665]
alpha = 0.006

delta_m = np.zeros(np.size(E_m))
for j in m:
    delta_m[j] = E_m[j] - E_0

T = 3/0.02419

def fit_dipoles(omegaVars, m1, b1):
    #dipole_fit = betaA*(np.sin(betaB*omegaVars))\
    #*(np.sin(betaB*omegaVars))*np.exp(-betaC*omegaVars)
    #dipole_fit = betaA*((1/(omegaVars-betaB))+\
    #(np.exp(-betaC*((omegaVars-betaB)*(omegaVars-betaB)))))
    dipole_fit = m1*omegaVars + b1
    return dipole_fit

def fit_dipoles3(omegaVars, m3, b3):
    #dipole_fit = betaA*(np.sin(betaB*omegaVars))\
    #*(np.sin(betaB*omegaVars))*np.exp(-betaC*omegaVars)
    #dipole_fit = (betaA/(omegaVars-betaB))+\
    #(betaC*np.exp(-betaD*((omegaVars-betaE)*(omegaVars-betaE))))\
    #+ (betaF*np.exp(-betaG*((omegaVars-betaH)*(omegaVars-betaH)))) + betaI
    #dipole_fit = dipole_fit = betaA*((1/(omegaVars-betaB))+\
    #(np.exp(-betaC*((omegaVars-betaB)*(omegaVars-betaB)))))
    dipole_fit = m3*omegaVars + b3
    
    return dipole_fit

def fit_alphas_dynamic(omegaV, m_1, b_1, m_3, b_3, alpha_free):
                       
    dawsarg_1st = T*(delta_m[0] - omegaV)
    dawsarg_3rd = T*(delta_m[1] - omegaV)
    #dawsarg_5th = T*(delta_m[2] - omegaV)
    alpha1f = fit_dipoles(omegaV, m_1, b_1)
    alpha3f = fit_dipoles3(omegaV, m_3, b_3)
    real_factor = (alpha1f*np.exp(-alpha_free*T*T*(delta_m[0]-omegaV)\
                    *(delta_m[0]-omegaV)))\
                     + (alpha3f* np.exp(-alpha_free*T*T*\
                    (delta_m[1]-omegaV)*(delta_m[1]-omegaV)))
    imag_factor = (alpha1f*(-2*cmath.sqrt(-1))/\
		(np.sqrt(np.pi))*spcl.dawsn(np.sqrt(alpha_free)*dawsarg_1st)) + \
			(alpha3f*((-2*cmath.sqrt(-1))/\
			(np.sqrt(np.pi)))*spcl.dawsn(np.sqrt(alpha_free)*dawsarg_3rd))
    domegaV = np.gradient(omegaV)
    phase_fit = np.arctan(np.imag(imag_factor)/np.real(real_factor))
    
    dphi_fit = np.gradient(phase_fit, domegaV)
    
    #return phase_fit
    return dphi_fit  