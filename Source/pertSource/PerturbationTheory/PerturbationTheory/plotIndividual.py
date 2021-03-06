import numpy as np
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
reAlphaTwo = np.zeros(sizeO)
imAlphaTwo = np.zeros(sizeO)
alphaThree = np.zeros(sizeO)
recf = np.zeros(sizeO)
imcf = np.zeros(sizeO)
recfOne = np.zeros(sizeO)
imcfOne = np.zeros(sizeO)
recfTwo = np.zeros(sizeO)
imcfTwo = np.zeros(sizeO)
recfThree = np.zeros(sizeO)
imcfThree = np.zeros(sizeO)

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
with open(Path + 'alphaTwo.txt') as infile:
    for line in infile:
        reAlphaTwo[i] = (line.split()[0])
        imAlphaTwo[i] = (line.split()[0])
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
i = 0
with open(Path + 'recfOne.txt') as infile:
    for line in infile:
        recfOne[i] = (line.split()[0])
        i = i + 1
        
i = 0
with open(Path + 'imcfOne.txt') as infile:
    for line in infile:
        imcfOne[i] = (line.split()[0])
        i = i + 1

i = 0
with open(Path + 'recfTwo.txt') as infile:
    for line in infile:
        recfTwo[i] = (line.split()[0])
        i = i + 1
        
i = 0
with open(Path + 'imcfTwo.txt') as infile:
    for line in infile:
        imcfTwo[i] = (line.split()[0])
        i = i + 1

i = 0
with open(Path + 'recfThree.txt') as infile:
    for line in infile:
        recfThree[i] = (line.split()[0])
        i = i + 1
        
i = 0
with open(Path + 'imcfThree.txt') as infile:
    for line in infile:
        imcfThree[i] = (line.split()[0])
        i = i + 1

domega = np.gradient(omega)

drecfOne = np.gradient(recfOne, domega)
dimcfOne = np.gradient(imcfOne, domega)
cfSquareOne = recfOne**2 + imcfOne**2
delayOne = (recfOne*dimcfOne - imcfOne*drecfOne)/cfSquareOne

drecfTwo = np.gradient(recfTwo, domega)
dimcfTwo = np.gradient(imcfTwo, domega)
cfSquareTwo = recfTwo**2 + imcfTwo**2
delayTwo = (recfTwo*dimcfTwo - imcfThree*drecfTwo)/cfSquareTwo

drecfThree = np.gradient(recfThree, domega)
dimcfThree = np.gradient(imcfThree, domega)
cfSquareThree = recfThree**2 + imcfThree**2
delayThree = (recfThree*dimcfThree - imcfThree*drecfThree)/cfSquareThree

# plt.plot(omega, cfSquareThree)
# plt.plot(omega, reAlphaTwo, omega, imAlphaTwo)
# plt.plot(omega, recfTwo, omega, imcfTwo)
# plt.plot(omega, alphaOne, 'r-', omega, alphaThree, 'b-')
# plb.ylim([-10, 10])
# plb.ylabel('alphas')
plt.plot(omega, np.gradient(np.arctan((imcfThree + imcfTwo)/(recfTwo + recfThree)), domega))
# plt.plot(omega, recf, omega, imcf )
# plb.ylim([-10, 10])
plb.ylabel('real and imag cf, just two')
plb.xlabel('central frequency (a.u.)')
plb.legend(['real', 'imag'])
plt.show()
