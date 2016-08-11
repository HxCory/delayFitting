import os
import cmath
import scipy.special as spcl
import scipy.optimize as sopt
import numpy as np
import pylab
import matplotlib.pylab as plb
import matplotlib.pyplot as plt

omega = np.arange(0.5, 1.2, 0.005)
m = range(1, 3)

c01 = 17.2795915535
c03 = 0.342313948196
E_0 = -1.1591
E_m = [-0.8502, -0.3278, -0.1665]
T = 3/0.02419
alpha = 0.006

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

imagCf1p3 = (c01 * np.imag(imagCf1)) + (c03 * np.imag(imagCf3))
realCf1p3 = (c01 * realCf1) + (c03 * realCf3)

# plt.plot(omega, imagCf1p3, omega, realCf1p3)
plt.plot(omega, c01 * np.imag(imagCf1), omega, c03 * np.imag(imagCf3)) 
plb.legend(['1st', '3rd'])
plb.show()