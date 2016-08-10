import os
import cmath
import numpy as np
import matplotlib.pylab as plb
import matplotlib.pyplot as plt

Path = '/Users/cgoldsmith/repos/delayFitting/Data/eigenstates/wavefxCalc__2'
os.chdir(Path)

wfsize = 0

with open('wf1D0.txt') as infile:
	for line in infile:
		row = (line.split())[0]
		wfsize = wfsize + 1
print wfsize

z = np.zeros(wfsize)
wf0 = np.zeros(wfsize)
wf1 = np.zeros(wfsize)
wf2 = np.zeros(wfsize)
wf3 = np.zeros(wfsize)

i = 0
# state = 0

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
		wf0[i] = (line.split())[1]
		i = i + 1



plt.plot(z, wf0, z, wf3)
plb.legend(['0', '1', '2', '3'])
plb.xlim([-40, 40])
plb.xlabel('z (a.u.)')
plb.ylabel('Re (wf0)')
plb.show()
