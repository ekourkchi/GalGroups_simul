## Importing Important Python Libraries
import sys
import os
import subprocess
from time import time
from math import *
import random

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import ascii
from astropy.table import Table, Column 
from scipy  import interpolate

# **************************************
# Global Variables
# Physical Constants
# **************************************
H0 = 75.           # hubble constant
sglV = 102.8806    # M87 center - super galactic longitude
sgbV = -2.3479     # M87 center - super galactic latitude
G = 6.67384E-11   # Gravitational constant
M_sun = 1.9891E30  # Solar maas [Kg]
Mpc_km = 3.08567758E19
#t0_age = Mpc_km / H0   # Universe Age [sec]
Yr= 365.25 * 24 * 3600    # [sec]
t0_age = 13.8E9 * Yr # Universe Age [sec] - 13.8 Gyr
h = 0.75  # hubble constant
# **************************************


# This part is responsible for the curved part 
#################################################################
# Originally made by M_L_ratio_curve_v2.py
table = np.genfromtxt('M_L_curve_v2.csv' , delimiter=',', filling_values=0, names=True, dtype=None)
Lumin_in   = 10**table['log10_L']
Mass_out   = 10**table['log10_M']
LtoM_func = interpolate.interp1d(Lumin_in, Mass_out)
#################################################################


# Returns mass in solar unit, given the absolute K-band luminosity
def Mass(L_k):
  
  if L_k==0:
    return 0
  
  L = L_k / 1.E10
  
  # ORIGINAL
  if L <= 0.0927:
    MtoL = 32.0*(L**-0.5)     # Original
  elif L > 0.0927 and L < 4.423:
    #MtoL = 58.0*(L**-0.25)
    return LtoM_func(L_k)
  elif L >= 4.423:
    MtoL = 32*(L**0.15)
  
  Mass_out = L_k * MtoL
  
  return Mass_out


#################################################################

if __name__ == '__main__':
  

  fig = plt.figure(figsize=(45/7., 6), dpi=100)
  ax = fig.add_axes([0.15, 0.13, 0.80,  0.80])  
  
  
  L = np.logspace(7,14,100)
  M = np.zeros(len(L))
  for i in range(len(L)):
      M[i] = Mass(L[i])
  
  ax.plot(L, M/L, 'b.')
  
          
  
  ax.set_ylabel('M'+r'$_v$'+r'$^{exp}$'+'/L'+r'$_{K_s}$'+' ['+r'$M_\odot/L_\odot$'+']', fontsize=16)
  ax.set_xlabel('K'+r'$_s$'+'-band Luminosity ['+r'$L_\odot$'+']', fontsize=16)
  plt.tick_params(which='major', length=7, width=1.5)
  plt.tick_params(which='minor', length=4, color='#000033', width=1.0)     
  plt.yticks(fontsize=16)
  plt.xticks(fontsize=16)
  ax.annotate(r'$L^{-0.5}$', (3.E7, 141), rotation=0, color='black', size=18)
  ax.annotate(r'$L^{0.15}$', (4.5E11, 35), rotation=0, color='black', size=18)
  
  plt.xscale('log')
  plt.yscale('log')
  plt.xlim(1.E7,1.E13)
  plt.ylim(1,1.E4)

  plt.show()  
  
  



