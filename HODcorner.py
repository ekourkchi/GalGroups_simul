from numpy import *
import numpy as np
from math import sqrt, pi, sin, cos, log as ln, e, log10, exp
import scipy.optimize as op
import emcee
import sys
import corner
import matplotlib as plt


filename = 'results_seed73_nchain3000_3pars.txt'
dat = loadtxt(filename)

naccept = dat[:,0]

fsat = dat[:,1]

Mm = dat[:,2]
delm = dat[:,3]
chi2 = dat[:,4]

Mm= np.asarray(Mm)

samples2=np.zeros((len(fsat),3))
samples2[:,0]= fsat
samples2[:,1]= delm
samples2[:,2]= Mm

fsat_true, mass_true , delm_true = np.mean(samples2[:,0]),np.mean(samples2[:,2]),np.mean(samples2[:,1])

print fsat_true, mass_true , delm_true

fig = corner.corner(samples2, bins=15, labels=[r"$f_{sat}$", r"$\Delta_{m}$", r"$M_{m}$"],truths=[fsat_true,delm_true, mass_true],color='green',fontsize=25,truth_color='red',scale_hist=False,space=0)

print fsat_true, mass_true, delm_true



fig.show()
fig.savefig("triangle2.eps")
