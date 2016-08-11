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

def gaussian(x, amp, cen, wid):
    "1-d gaussian: gaussian(x, amp, cen, wid)"
    return (amp/(sqrt(2*pi)*wid)) * exp(-(x-cen)**2 /(2*wid**2))


gmod = Model(func.fitDipoleFactors)
result = gmod.fit(y, omegaF = x,\
 c_1 = -9.88389208,  c_3 = 4.94309140)

print(result.fit_report())

plt.plot(x, y,'bo')
# plt.plot(x, result.init_fit, 'k--')
plt.plot(x, result.best_fit, 'r-')
# plt.plot(defs.omega_dip, func.fitDipoleFactors\
# 	(defs.omega_dip, -9.88230074, 4.94203098, 0.00187428 ))
plt.show()