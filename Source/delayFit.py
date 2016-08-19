import numpy as np
import matplotlib.pyplot as plt
import pylab as plb
import defDelay as defD
import scipy.optimize as sopt
from lmfit import Model
import os

Path = '/Users/cgoldsmith/Desktop/text_files_data'
os.chdir(Path)
data = np.loadtxt('TDSE_3fs.txt')
x = data[:, 0]
y = data[:, 4]

gmod = Model(defD.delay)
result = gmod.fit(y, omegaV = x, free = 0.002)

print(result.fit_report())
plt.plot(x, y, 'bo')
plt.plot(x, result.init_fit, 'k--')
plt.plot(x, result.bes_fit, 'r-')
plt.show()