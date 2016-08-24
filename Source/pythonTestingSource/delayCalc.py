import os
import cmath
import scipy.special as spcl
import scipy.optimize as sopt
import numpy as np
import pylab
import matplotlib.pylab as plb
import matplotlib.pyplot as plt
import defsDipole as defs
import twoPhotonDipole as func

omega = np.arange(0.5, 1.2, 0.0005)
m = range(1, 3)
mu1 = np.zeros(np.size(omega), 'complex')
mu3 = np.zeros(np.size(omega), 'complex')
l = 0
for i in omega:
	mu1[l] = func.constructDipole(i, defs.wf1)
	mu3[l] = func.constructDipole(i, defs.wf3)
	l = l + 1
E_0 = -1.1591
E_m = [-0.8502, -0.3278, -0.1665]
T = 3/0.02419
alpha = 1.0

delta_m = np.zeros(np.size(E_m))
for j in m:
    delta_m[j] = E_m[j] - E_0

dawsarg_1st = T*(delta_m[0] - omega)
dawsarg_3rd = T*(delta_m[1] - omega)
dawsarg_5th = T*(delta_m[2] - omega)

imagCf1 = (-2*cmath.sqrt(-1))/(np.sqrt(np.pi))\
				*spcl.dawsn(np.sqrt(alpha)*dawsarg_1st)
imagCf3 = (-2*cmath.sqrt(-1))/(np.sqrt(np.pi))\
				*spcl.dawsn(np.sqrt(alpha)*dawsarg_3rd)
imagCf5 = (-2*cmath.sqrt(-1))/(np.sqrt(np.pi))\
				*spcl.dawsn(np.sqrt(alpha)*dawsarg_5th)

realCf1 = np.exp(-alpha*T*T*(delta_m[0]-omega)\
                    *(delta_m[0]-omega))
realCf3 = np.exp(-alpha*T*T*(delta_m[1]-omega)\
                    *(delta_m[1]-omega))
realCf5 = np.exp(-alpha*T*T*(delta_m[2]-omega)\
                    *(delta_m[2]-omega))

cf1 = mu1 * (realCf1 + imagCf1)
cf3 = mu3 * (realCf3 + imagCf3)

# plt.plot(omega, np.imag(cf1 + cf3), omega, np.real(cf1 + cf3))
# plb.legend(['imag', 'real'])

quotientSum = np.imag(cf1 + cf3)/np.real(cf1 + cf3)
plt.plot(omega, quotientSum)

# phaseSum = np.arctan(quotientSum)
# dOmega = np.gradient(omega)
# dPhiSum = np.gradient(phaseSum,dOmega)

# plt.plot(omega, phaseSum)
# plt.plot(omega, dPhiSum)
plb.show()