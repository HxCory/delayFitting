import os
import cmath
import numpy as np
import matplotlib.pylab as plb
import matplotlib.pyplot as plt

def getk(omega, IP):
	k =  np.sqrt(2 * ((2 * omega) + IP))
	return k

def conjugate(z):
	return np.real(z) - (cmath.sqrt(-1) * np.imag(z))

Path = '/Users/cgoldsmith/repos/delayFitting/Data/eigenstates/wavefxCalc__2'
os.chdir(Path)

m = range(1, 3)
E_0 = -1.1591
E_m = [-0.8502, -0.3278, -0.1665]

omegaAU = 0.67

wfsize = 0


with open('wf1D0.txt') as infile:
	for line in infile:
		row = (line.split())[0]
		wfsize = wfsize + 1

z = np.zeros(wfsize)
wf0Re = np.zeros(wfsize)
wf0Im = np.zeros(wfsize)
wf0 = np.zeros(wfsize, 'complex')
wf1Re = np.zeros(wfsize)
wf1Im = np.zeros(wfsize)
wf1 = np.zeros(wfsize, 'complex')
wf2Re = np.zeros(wfsize)
wf2Im = np.zeros(wfsize)
wf2 = np.zeros(wfsize, 'complex')
wf3 = np.zeros(wfsize, 'complex')
wf3Re = np.zeros(wfsize)
wf3Im = np.zeros(wfsize)
planeWaveTest = np.zeros(wfsize, 'complex')
zRealPlaneWave = np.zeros(wfsize)
zImagPlaneWave = np.zeros(wfsize)

i = 0
with open('wf1D0.txt') as infile:
	for line in infile:
		z[i] = (line.split())[0]
		wf0Re[i] = (line.split())[1]
		wf0Im[i] = (line.split())[2]
		wf0[i] = wf0Re[i] + cmath.sqrt(-1) * wf0Im[i]
		i = i + 1

i = 0
with open('wf1D1.txt') as infile:
	for line in infile:
		wf1Re[i] = (line.split())[1]
		wf1Im[i] = (line.split())[2]
		wf1[i] = wf2Re[i] + cmath.sqrt(-1) * wf2Im[i]
		i = i + 1

i = 0
with open('wf1D2.txt') as infile:
	for line in infile:
		wf2Re[i] = (line.split())[1]
		wf2Im[i] = (line.split())[2]
		wf2[i] = wf2Re[i] + cmath.sqrt(-1) * wf2Im[i]
		i = i + 1

i = 0
with open('wf1D3.txt') as infile:
	for line in infile:
		wf3Re[i] = (line.split())[1]
		wf3Im[i] = (line.split())[2]
		wf3[i] = wf3Re[i] + cmath.sqrt(-1) * wf3Im[i]
		i = i + 1

i = 0
while(i < wfsize):
	for m in z:
		planeWaveTest[i] = np.exp(-1 * cmath.sqrt(-1) * m * getk(omegaAU, E_0))
		i = i + 1
	

m = 0
for line in z:
	zRealPlaneWave[m] = np.real(planeWaveTest[m]) * z[m]
	zImagPlaneWave[m] = np.imag(planeWaveTest[m]) * z[m]
	m = m +1

plt.plot(z, wf2Re, z, wf2Im)
plb.xlim([-10, 10])
plt.show()

