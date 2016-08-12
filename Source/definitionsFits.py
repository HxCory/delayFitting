import os
import cmath
import scipy.special as spcl
import scipy.optimize as sopt
import numpy as np
import matplotlib.pylab as plb
import matplotlib.pyplot as plt

omega = np.arange(0.5, 1.2, 0.005)
m = range(1, 3)
E_0 = -1.1591
E_m = [-0.8502, -0.3278, -0.1665]
alpha = 0.006
T = 3/0.02419

Path = '/Users/cgoldsmith/Desktop/text_files_data'
    # os.chdir(Path)
omega_dip=np.zeros(2000)
i=0
with open(Path+'/omega4dip.txt') as infile:
    for line in infile:
        omega_dip[i]=(line.split()[0])
        i=i+1

dip_1=np.zeros(np.size(omega_dip))
i=0
with open(Path + '/alpha1_3fs.txt') as infile:
    for line in infile:
        dip_1[i]=(line.split()[0])
        i=i+1    

dip_3=np.zeros(np.size(omega_dip))
i=0
with open(Path + '/alpha3_3fs.txt') as infile:
    for line in infile:
        dip_3[i]=(line.split()[0])
        i=i+1       

delta_m = np.zeros(np.size(E_m))
for j in m:
    delta_m[j] = E_m[j] - E_0

def fit_dipoles(omegaVars, m1, b1):
    dipole_fit = m1*omegaVars + b1
    return dipole_fit

def fit_dipoles3(omegaVars, m3, b3):
    dipole_fit = m3*omegaVars + b3 
    return dipole_fit

def fit_alphas_dynamic(omegaV, m_1, b_1, m_3, b_3, alpha_free):
                       
    dawsarg_1st = T*(delta_m[0] - omegaV)
    dawsarg_3rd = T*(delta_m[1] - omegaV)
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
    return dphi_fit  

# def fit_alphas_dynamic(omegaF, c_1, c_3, alpha_free):
#     ideal_size = 2000
#     opt_size = np.size(omegaF)
#     domega_dip = np.gradient(omega_dip)
#     alpha1f = np.zeros(opt_size)
#     alpha3f = np.zeros(opt_size)
#     omegaV = np.zeros(opt_size)
#     real_factor = np.zeros(opt_size)
#     imag_factor = np.zeros(opt_size)
    
#     if opt_size == ideal_size:
#         j = 0
#         alpha1f = dip_1 + c_1
#         alpha3f = dip_3 + c_3
#         omegaV = omega_dip
#     else:
#         fac = int(ideal_size/(opt_size-1))
#         looper = np.arange(0, (ideal_size), fac)
#         i = 0
#     #return(loopertest)
#         for x in looper:
#             omegaV[i] = omega_dip[x]
#             alpha1f[i] = dip_1[x] + c_1
#             alpha3f[i] = dip_3[x] + c_3
#             i = i + 1
#     dawsarg_1st = T*(delta_m[0] - omegaV)
#     dawsarg_3rd = T*(delta_m[1] - omegaV)
#     #dawsarg_5th = T*(delta_m[2] - omegaV)
            
#     real_factor = (alpha1f*np.exp(-alpha_free*alpha_free*T*T*(delta_m[0]-omegaV)\
#                     *(delta_m[0]-omegaV)))\
#                      + (alpha3f* np.exp(-alpha_free*T*T*\
#                     (delta_m[1]-omegaV)*(delta_m[1]-omegaV)))
#     imag_factor = (alpha1f*(-2*cmath.sqrt(-1))/\
#     (np.sqrt(np.pi))*spcl.dawsn(np.sqrt(alpha_free)*dawsarg_1st)) + \
#     (alpha3f*((-2*cmath.sqrt(-1))/\
#     (np.sqrt(np.pi)))*spcl.dawsn(alpha_free*dawsarg_3rd))
#     domegaV = np.gradient(omegaV)
#     phase_fit = np.arctan(np.imag(imag_factor)/np.real(real_factor))
    
#     dphi_fit = np.gradient(phase_fit, domegaV)
    
#     #return omegaV
#     return dphi_fit

