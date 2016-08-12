from numpy import sqrt, pi, exp, linspace, loadtxt
import numpy as np
from lmfit import  Model
import os
import matplotlib.pyplot as plt
import definitionsFits as defs
import twoPhotonDipole as func
import tooManyParameters as too

Path = '/Users/cgoldsmith/Desktop/text_files_data'
os.chdir(Path)
data = loadtxt('TDSE_3fs.txt')
x = data[:, 0]
y = data[:, 4]

gmod = Model(too.fitAlphas)
result = gmod.fit(y, omegaV = x,
	b1A = -9.81535556e+00,
	b1B = 3.13975075e-01,
	b1C = 4.94111384e+02,
	b1D = 4.29958163e-01,
	b1E = -2.09667678e+00,
	b1F = 5.85035181e+00, 
	b3A = 3.26519717e-01,
	b3B = -3.86939149e-01,
	b3C = -6.38185123e+00,
	b3D = 9.69371148e+00,
	b3E = 6.33599497e-02,
	b3F = -1.60485987e-01,
	b3G = -1.36571022e-01,
	b3H = 1.15095946e-02,
	b3I = -7.15781666e-02,
	alpha_free = 5.62392912e-03)

# gmod = Model(too.presetAlphas)
# result = gmod.fit(y, omegaV = x, alpha_free = 0.00562)

print(result.fit_report())

plt.plot(x, y,'bo')
plt.plot(x, result.init_fit, 'k--')
plt.plot(x, result.best_fit, 'r-')
plt.show()
