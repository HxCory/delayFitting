import os
import cmath
import numpy as np
import defsDipole as defs
import scipy.special as spcl


def constructBound(wf):
	i = 0
	dip = 0
	for m in defs.z:
		elmt = -m * defs.wf0[i] * defs.conjugate(wf[i])
		dip = dip + elmt
		i = i + 1
	return dip

def constructPlaneWave(omegaV):
    i = 0
    planeWave = np.zeros(defs.wfsize, 'complex')
    while(i < defs.wfsize):
        for m in defs.z:
            planeWave[i] = np.exp(-1 * cmath.sqrt(-1) * m * defs.getk(omegaV, defs.E_0))
            i = i + 1
    return planeWave

def constructContinuum(omegar, wf):
    realDipper = np.zeros(defs.wfsize)
    imagDipper = np.zeros(defs.wfsize)
    dummyVec = constructPlaneWave(omegar)

    realDipper = -np.real(dummyVec * wf) * defs.z 
    imagDipper = -np.imag(dummyVec * wf) * defs.z

    summer = sum(realDipper) + (cmath.sqrt(-1) * sum(imagDipper))    
    return summer


def constructDipole(photonE, wfInt):
	return constructContinuum(photonE, wfInt) * constructBound(wfInt)

def fitDipoleFactors(omegaF, c_1, c_3):
	E_0 = -1.1591
	E_m = [-0.8502, -0.3278, -0.1665]
	T = 3/0.02419
	delta_m = np.zeros(np.size(E_m))
	for j in (0, 1, 2):
	    delta_m[j] = E_m[j] - E_0
	optSize = np.size(omegaF)
	alpha1f = np.zeros(optSize, 'complex')
	alpha3f = np.zeros(optSize, 'complex')
	omegaV = np.zeros(optSize)
	realFactor = np.zeros(optSize, 'complex')
	imagFactor = np.zeros(optSize, 'complex')
	i = 0
	T = 3 / 0.02419

	for x in omegaF:
		omegaV[i] = x
		alpha1f[i] = constructDipole(x, defs.wf1) + c_1
		alpha3f[i] = constructDipole(x, defs.wf3) + c_3
		i = i + 1

	dawsArg1st = T*(delta_m[0] - omegaV)
	dawsArg3rd = T*(delta_m[1] - omegaV)

	realFactor = (alpha1f * np.exp(-np.power(dawsArg1st,2.0)))\
					+ (alpha3f * np.exp(-np.power(dawsArg3rd, 2.0)))

	imagFactor = (alpha1f * ((-2 * cmath.sqrt(-1))\
					/ (np.sqrt(np.pi))) * spcl.dawsn(dawsArg1st))\
						+ (alpha3f * ((-2 * cmath.sqrt(-1))\
							/ (np.sqrt(np.pi))) * spcl.dawsn(dawsArg3rd))

	domegaV = np.gradient(omegaV)
	phase_fit = np.arctan(np.imag(imagFactor)/np.real(realFactor))

	dphi_fit = np.gradient(phase_fit, domegaV)

	#return omegaV
	return dphi_fit
