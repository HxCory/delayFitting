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
with open(Path + 'recf.txt') as infile:
    for line in infile:
        recf[i] = (line.split()[0])
        i = i + 1
        
i = 0
with open(Path + 'imcf.txt') as infile:
    for line in infile:
        imcf[i] = (line.split()[0])
        i = i + 1


domega = np.gradient(omega)
drecf = np.gradient(recf, domega)
dimcf = np.gradient(imcf, domega)

cfSquare = (recf + cmath.sqrt(-1)*imcf) * (recf - cmath.sqrt(-1)*imcf)
delay = (recf*dimcf - imcf*drecf)

plt.plot(omega, -delay/cfSquare, 'r-', x, y, 'bo')
plb.ylim([-20, 20])
plt.show()
