import numpy as np
import definitionsFits as defs
import matplotlib.pyplot as plt
import cmath
import os
import scipy.special as spcl
import pylab as plb
import twoPhotonDipole as dipole
import defsDipole
import scipy.optimize as sopt

def delay(omegaV, free):
	E_0 = -1.1591
	E_m = [-0.8502, -0.3278, -0.1665]
	T = 3/0.02419
	delta_m = np.zeros(np.size(E_m))
	for j in (0, 1, 2):
	    delta_m[j] = E_m[j] - E_0
	optSize = np.size(omegaV)
	alpha1f = np.zeros(optSize, 'complex')
	alpha3f = np.zeros(optSize, 'complex')
	realFactor = np.zeros(optSize, 'complex')
	imagFactor = np.zeros(optSize, 'complex')
	T = 3 / 0.02419

	for i in range(0, np.size(defs.omega_dip)-1):
	    alpha1f[i] = dipole.constructDipole(omegaV[i], defsDipole.wf1)
	    alpha3f[i] = dipole.constructDipole(omegaV[i], defsDipole.wf3)

	dawsArg1st = T*(delta_m[0] - omegaV)
	dawsArg3rd = T*(delta_m[1] - omegaV)
	realFactor = np.real((alpha1f * np.exp(-free*np.power(dawsArg1st,2.0)))\
            + (alpha3f * np.exp(-np.sqrt(free)*np.power(dawsArg3rd, 2.0))))

	imagFactor = (alpha1f * ((-2 * cmath.sqrt(-1))\
            / (np.sqrt(np.pi))) * spcl.dawsn(free*dawsArg1st))\
                + (alpha3f * ((-2 * cmath.sqrt(-1))\
                    / (np.sqrt(np.pi))) * spcl.dawsn(np.sqrt(free)*dawsArg3rd))

	dE = np.gradient(omegaV)
	dIm = np.gradient(imagFactor, dE)
	dRe = np.gradient(realFactor, dE)
	derivFactor = (realFactor * dIm) - (np.imag(imagFactor * dRe))
	zSquared = realFactor*realFactor + np.imag(imagFactor)*np.imag(imagFactor)
	ans =  (realFactor * np.imag(dIm) - np.imag(imagFactor * dRe)) / zSquared
	return ans

print delay(defs.omega_dip, 0.002)