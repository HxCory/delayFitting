import numpy as np
import definitionsFits as defs
import matplotlib.pyplot as plt
import cmath
import os
import scipy.special as spcl
import pylab


Path = '/Users/cgoldsmith/Desktop/text_files_data'
os.chdir(Path)
data = np.loadtxt('TDSE_3fs.txt')
x = data[:, 0]
y = data[:, 4]

def alphaOne(omega, bA, bB, bC, bD, bE, bF):
    termOne = bA / (omega - bB)
    termTwo = bC * np.exp(-bD * np.power(omega - bE, 2.0))
    return termOne + termTwo + bF

def alphaThree(omega, bA, bB, bC, bD, bE, bF, bG, bH, bI):
    termOne = bA / (omega - bB)
    termTwo = bC * np.exp(-bD * np.power(omega - bE, 2.0))
    termThree = bF * np.exp(-bG * np.power(omega - bH, 2.0))
    return termOne + termTwo + termThree + bI

def presetAlphas(omegaV, alpha_free):                       
    dawsarg_1st = defs.T*(defs.delta_m[0] - omegaV)
    dawsarg_3rd = defs.T*(defs.delta_m[1] - omegaV)
    alpha1f = alphaOne(omegaV, -9.81535556e+00,   3.13975075e-01,   4.94111384e+02,
          4.29958163e-01,  -2.09667678e+00,   5.85035181e+00)
    alpha3f = alphaThree(omegaV,\
          3.26519717e-01,  -3.86939149e-01,  -6.38185123e+00,
          9.69371148e+00,   6.33599497e-02,  -1.60485987e-01,
         -1.36571022e-01,   1.15095946e-02,  -7.15781666e-02)
          
    real_factor = (alpha1f*np.exp(-alpha_free*defs.T*defs.T*(defs.delta_m[0]-omegaV)\
                    *(defs.delta_m[0]-omegaV)))\
                     + (alpha3f * np.exp(-alpha_free* defs.T * defs.T *\
                    (defs.delta_m[1]-omegaV)*(defs.delta_m[1]-omegaV)))
    imag_factor = (alpha1f*(-2*cmath.sqrt(-1))/\
        (np.sqrt(np.pi))*spcl.dawsn(np.sqrt(alpha_free)*dawsarg_1st)) + \
            (alpha3f*((-2*cmath.sqrt(-1))/\
                (np.sqrt(np.pi)))*spcl.dawsn(np.sqrt(alpha_free)*dawsarg_3rd))
    domegaV = np.gradient(omegaV)
    phase_fit = np.arctan(np.imag(imag_factor)/np.real(real_factor))
    
    dphi_fit = np.gradient(phase_fit, domegaV)
    return dphi_fit  

def fitAlphas(omegaV, b1A, b1B, b1C, b1D, b1E, b1F, \
		b3A, b3B, b3C, b3D, b3E, b3F, b3G, b3H, b3I, alpha_free):                       
    dawsarg_1st = defs.T*(defs.delta_m[0] - omegaV)
    dawsarg_3rd = defs.T*(defs.delta_m[1] - omegaV)
    alpha1f = alphaOne(omegaV, b1A, b1B, b1C, b1D, b1E, b1F)
    alpha3f = alphaThree(omegaV, b3A, b3B, b3C, b3D, b3E, b3F, b3G, b3H, b3I)          
    real_factor = (alpha1f*np.exp(-alpha_free*defs.T*defs.T*(defs.delta_m[0]-omegaV)\
                    *(defs.delta_m[0]-omegaV)))\
                     + (alpha3f * np.exp(-alpha_free* defs.T * defs.T *\
                    (defs.delta_m[1]-omegaV)*(defs.delta_m[1]-omegaV)))
    imag_factor = (alpha1f*(-2*cmath.sqrt(-1))/\
        (np.sqrt(np.pi))*spcl.dawsn(np.sqrt(alpha_free)*dawsarg_1st)) + \
            (alpha3f*((-2*cmath.sqrt(-1))/\
                (np.sqrt(np.pi)))*spcl.dawsn(np.sqrt(alpha_free)*dawsarg_3rd))
    domegaV = np.gradient(omegaV)
    phase_fit = np.arctan(np.imag(imag_factor)/np.real(real_factor))
    
    dphi_fit = np.gradient(phase_fit, domegaV)
    return dphi_fit 

sampleOne = presetAlphas(defs.omega_dip, 5.62392912e-03)
sampleTwo = fitAlphas(x, -9.81535556e+00,   3.13975075e-01,   4.94111384e+02,
          4.29958163e-01,  -2.09667678e+00,   5.85035181e+00, 3.26519717e-01,  -3.86939149e-01,  -6.38185123e+00,
          9.69371148e+00,   6.33599497e-02,  -1.60485987e-01,
         -1.36571022e-01,   1.15095946e-02,  -7.15781666e-02, 5.62392912e-03)
plt.plot(defs.omega_dip, sampleOne, 'r-', x, y, 'bo', x, sampleTwo)
pylab.ylim([-10, 10])
plt.show()