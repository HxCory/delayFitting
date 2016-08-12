import os
import cmath
import numpy as np
import matplotlib.pylab as plb
import matplotlib.pyplot as plt

def getk(omega, IP):
	k =  np.sqrt(2 * ((2 * omega) - IP))
	return k

def conjugate(z):
	return np.real(z) - (cmath.sqrt(-1) * np.imag(z))

Path = '/Users/cgoldsmith/repos/delayFitting/Data/eigenstates/wavefxCalc__2'
os.chdir(Path)

m = range(1, 3)
E_0 = -1.1591
E_m = [-0.8502, -0.3278, -0.1665]

omegaAU = 1.15

wfsize = 0


with open('wf1D0.txt') as infile:
	for line in infile:
		row = (line.split())[0]
		wfsize = wfsize + 1

z = np.zeros(wfsize)
wf0 = np.zeros(wfsize)
wf1 = np.zeros(wfsize)
wf2 = np.zeros(wfsize)
wf3 = np.zeros(wfsize)
planeWaveTest = np.zeros(wfsize, 'complex')
zRealPlaneWave = np.zeros(wfsize)
zImagPlaneWave = np.zeros(wfsize)

i = 0
with open('wf1D0.txt') as infile:
	for line in infile:
		z[i] = (line.split())[0]
		wf0[i] = (line.split())[1]
		i = i + 1

i = 0
with open('wf1D1.txt') as infile:
	for line in infile:
		wf1[i] = (line.split())[1]
		i = i + 1

i = 0
with open('wf1D2.txt') as infile:
	for line in infile:
		wf2[i] = (line.split())[1]
		i = i + 1

i = 0
with open('wf1D3.txt') as infile:
	for line in infile:
		wf3[i] = (line.split())[1]
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


