import os
import cmath
import numpy as np
import matplotlib.pylab as plb
import matplotlib.pyplot as plt
import defsDipole as dip

# plt.plot(dip.z, dip.wf0, dip.z, dip.wf3)
# plb.legend(['0', '1', '2', '3'])
# plb.xlim([-40, 40])
# plb.xlabel('z (a.u.)')
# plb.ylabel('Re (wf0)')
# plb.show()
plt.plot(dip.z, dip.zRealPlaneWave, dip.z, dip.zImagPlaneWave)
plt.plot(dip.z, np.real(dip.planeWaveTest))
plb.xlim([-40, 40])
plb.show()
