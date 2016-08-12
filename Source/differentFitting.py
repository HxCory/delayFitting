from numpy import sqrt, pi, exp, linspace, loadtxt
import numpy as np
from lmfit import  Model
import os
import matplotlib.pyplot as plt
import definitionsFits as defs
import twoPhotonDipole as func

Path = '/Users/cgoldsmith/Desktop/text_files_data'
os.chdir(Path)
data = loadtxt('TDSE_3fs.txt')
x = data[:, 0]
y = data[:, 4]

gmod = Model(func.fitDipoleFactors)
result = gmod.fit(y, omegaF = x, m1 = 1.0, m3 = 1.0,\
	c_1 = -4.0088389208,  c_3 = 1.0094309140)

print(result.fit_report())

plt.plot(x, y,'bo')
plt.plot(x, result.init_fit, 'k--')
plt.plot(x, result.best_fit, 'r-')
plt.show()
