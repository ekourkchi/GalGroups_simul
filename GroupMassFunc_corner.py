from numpy import *
import numpy as np
from math import sqrt, pi, sin, cos, log as ln, e, log10, exp
import scipy.optimize as op
#import emcee
import sys
import corner
import matplotlib as plt


######################################################################
def median_error(X):
  
  size = len(X)
  X    = np.sort(X) 
  mean = np.median(X)
  
  u_err = X[int(round(0.84*size))] - mean
  l_err = mean - X[int(round(0.16*size))]
  
  return mean, u_err, l_err

######################################################################




#filename = 'v3seed200.csv'
filename = 'set4.csv'
table = np.genfromtxt(filename , delimiter=',', filling_values=0, names=True, dtype=None)

alfa = table['alfa']
beta = table['beta']
gama = table['gama']

alfa = alfa[100:]
beta = beta[100:]
gama = gama[100:]


alfa_med = np.median(alfa)
beta_med = np.median(beta)
gama_med = np.median(gama)
                     
alfa_std = np.std(alfa)
beta_std = np.std(beta)
gama_std = np.std(gama)

N = len(alfa)
f1 = np.zeros(N)
f2 = np.zeros(N)
f3 = np.zeros(N)
f4 = np.zeros(N)
f5 = np.zeros(N)
f6 = np.zeros(N)

c = 3.0

f1[np.where(alfa<alfa_med+c*alfa_std)] = 1
f2[np.where(alfa>alfa_med-c*alfa_std)] = 1
f3[np.where(beta<beta_med+c*beta_std)] = 1
f4[np.where(beta>beta_med-c*beta_std)] = 1
f5[np.where(gama<gama_med+c*gama_std)] = 1
f6[np.where(gama>gama_med-c*gama_std)] = 1

f = f1+f2+f3+f4+f5+f6
alfa = alfa[np.where(f==6)]
beta = beta[np.where(f==6)]
gama = gama[np.where(f==6)]

samples=np.zeros((len(alfa),3))
samples[:,0]= alfa
samples[:,1]= beta
samples[:,2]= gama

alf, bet, gam = np.median(samples[:,0]),np.median(samples[:,1]),np.median(samples[:,2])

al, al_uerr, al_lerr = median_error(alfa)
bt, bt_uerr, bt_lerr = median_error(beta)
ga, ga_user, ga_lerr = median_error(gama)





fig = corner.corner(samples, bins=25, labels=[r"$\alpha$", r"$\beta$", r"$\gamma$"],fontsize=25,truth_color='green',scale_hist=False,space=0, quantiles=[0.16, 0.84], levels=(1-np.exp(-0.5),1-np.exp(-2.0),), truths=[alf,bet,gam], show_titles=True, title_kwargs={"fontsize": 12})




#fig.show()
fig.savefig("M2L_fit_CornerPlot.eps")
