import os
import cmath
import numpy as np
import matplotlib.pylab as plb
import matplotlib.pyplot as plt
import defsDipole as dip

dip01 = 0
dip03 = 0
i = 0
for m in dip.z:
	elmt01 = -m * dip.wf0[i] * dip.conjugate(dip.wf1[i])
	elmt03 = -m * dip.wf0[i] * dip.conjugate(dip.wf3[i])
	dip01 = dip01 + elmt01 
	dip03 = dip03 + elmt03
	i = i + 1
print dip01
print dip03

plt.plot(dip.z, dip.zRealPlaneWave, dip.z, dip.zImagPlaneWave)
plb.xlim([-40, 40])
plb.ylim([-40, 40])
plb.show()

