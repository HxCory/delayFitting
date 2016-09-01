import numpy as np
import matplotlib.pyplot as plt
import pylab as plb
import cmath
import os

Path2 = '/Users/cgoldsmith/Desktop/text_files_data'
os.chdir(Path2)
data = np.loadtxt('TDSE_3fs.txt')
x = data[:, 0]
y = data[:, 4]

# PathTest = '/Users/cgoldsmith/repos/delayFitting/Data/pertOutput/Tests/'
PathTest = '/Users/cgoldsmith/repos/delayFitting/Data/pertOutput/'

def plotDelays(Path):
    os.chdir(Path)

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

    cfSquare = recf**2 + imcf**2
    delay = (recf*dimcf - imcf*drecf)/cfSquare

    plt.plot(omega, delay, 'r-', x, y, 'bo')
    plb.ylim([-10, 10])
    plb.ylabel('streaking delay (a.u.)')
    # plt.plot(omega, np.arctan(imcf/recf))
    # plb.ylabel('phase, first and third')
    # plb.xlabel('central frequency (a.u.)')
    # plb.legend(['real', 'imag'])
    # plt.plot(omega, alphaOne, omega, alphaThree)
    plt.show()

plotDelays(PathTest)
