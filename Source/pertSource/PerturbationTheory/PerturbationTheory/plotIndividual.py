import numpy as np
from scipy.interpolate import UnivariateSpline
import matplotlib.pyplot as plt
import pylab as plb
import cmath
import os

Path = '/Users/cgoldsmith/repos/delayFitting/Data/pertOutput/'
os.chdir(Path)

Path2 = '/Users/cgoldsmith/Desktop/text_files_data'
os.chdir(Path2)
data = np.loadtxt('TDSE_3fs.txt')
x = data[:, 0]
y = data[:, 4]

sizeO = 0
with open(Path + 'omega.txt') as infile:
    for line in infile:
        sizeO = sizeO + 1

omega = np.zeros(sizeO)
alphaOne = np.zeros(sizeO)
alphaThree = np.zeros(sizeO)
recf = np.zeros(sizeO)
imcf = np.zeros(sizeO)

i = 0
with open(Path + 'omega.txt') as infile:
    for line in infile:
        omega[i] = (line.split()[0])
        i = i + 1
        
i = 0
with open(Path + 'alphaOne.txt') as infile:
    for line in infile:
        alphaOne[i] = (line.split()[0])
        i = i + 1
        
i = 0
with open(Path + 'alphaThree.txt') as infile:
    for line in infile:
        alphaThree[i] = (line.split()[0])
        i = i + 1
        
i = 0
with open(Path + 'recfOne.txt') as infile:
    for line in infile:
        recf[i] = (line.split()[0])
        i = i + 1
        
i = 0
with open(Path + 'imcfOne.txt') as infile:
    for line in infile:
        imcf[i] = (line.split()[0])
        i = i + 1


domega = np.gradient(omega)
drecf = np.gradient(recf, domega)
dimcf = np.gradient(imcf, domega)

cfSquare = recf**2 + imcf**2
delay = (recf*dimcf - imcf*drecf)/cfSquare

# plt.plot(omega, cfSquare)

plt.plot(omega, delay, 'r-')
# plb.ylim([-10, 10])
# plb.ylabel('streaking delay (a.u.)')
# plb.xlabel('central frequency (a.u.)')
# plb.legend(['model', 'TDSE'])
plt.show()
