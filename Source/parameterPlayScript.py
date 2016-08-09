import os
import cmath
import scipy.special as spcl
import scipy.optimize as sopt
import numpy as np
import matplotlib.pylab as plb
import matplotlib.pyplot as plt
import pylab
import definitionsFits as fit
 
Path = '/Users/cgoldsmith/repos/delayFitting/Data'
os.chdir(Path)

omegaD=np.zeros(32)
i=0
with open('TDSE_3fs.txt') as infile:
    for line in infile:
        omegaD[i]=(line.split()[0])
        i=i+1

delays=np.zeros(32)
i=0
with open('TDSE_3fs.txt') as infile:
    for line in infile:
        delays[i]=(line.split()[4])
        i=i+1

omega_dip=np.zeros(2000)
i=0
with open('omega4dip.txt') as infile:
    for line in infile:
        omega_dip[i]=(line.split()[0])
        i=i+1

dip_1=np.zeros(np.size(omega_dip))
i=0
with open('alpha1_3fs.txt') as infile:
    for line in infile:
        dip_1[i]=(line.split()[0])
        i=i+1

dip_3=np.zeros(np.size(omega_dip))
i=0
with open('alpha3_3fs.txt') as infile:
    for line in infile:
        dip_3[i]=(line.split()[0])
        i=i+1

plt.plot(omegaD, delays, omega_dip, fit.fit_alphas_dynamic\
         (omega_dip,\
			-4.08147964, -16.96877514,   0.3398741 ,\
	          -0.21666972,   0.03361551))
pylab.ylim([-5, 5])
pylab.legend(['TDSE data', 'fit'])
pylab.xlabel('omega (a.u.)')
pylab.ylabel('time delay (a.u.)')
pylab.title('first order/linear approx. for dipoles')
pylab.show()


