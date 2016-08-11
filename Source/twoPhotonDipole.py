import os
import cmath
import numpy as np
import defsDipole as defs

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
