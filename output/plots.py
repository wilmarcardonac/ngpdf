import numpy as np
import pylab as py
from scipy import stats

Ak = np.loadtxt('Ak-distribution-simulations.dat')

P1 = np.loadtxt('one-point-distribution-functions.dat')

py.plot(Ak[:,0],Ak[:,1],label='simulations')

#for index in range(len(hist)):
#    print float(hist[index])/float(hist.sum())

#print hist.sum()

py.plot(P1[:,0],P1[:,1],label='Gaussian')

py.plot(P1[:,0],P1[:,2],label='Non-Gaussian')

py.legend()

py.show()

exit()
