import sys
import os
import random
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import numpy as np
from math import *
from time import time
import random
from astropy.io import ascii
from astropy.table import Table, Column 
import pyfits
import pylab as py

from astropy import coordinates as coord
from astropy import units as unit
import subprocess
from scipy  import interpolate

#################################################################
def median(target, array, mi, mf):
    
    array = np.asarray(array)
    N = len(array)
    
    f1 = np.zeros(N)
    f2 = np.zeros(N)
    
    f1[np.where(array>=mi)] = 1
    f2[np.where(array<mf)] = 1
    f = f1 + f2
    
    ML = target[np.where(f==2)]
    media = np.median(ML)
    stdev = np.std(ML)
    
    N = len(ML)
    f1 = np.zeros(N)
    f2 = np.zeros(N)
    f1[np.where(ML<=media+2*stdev)] = 1
    f2[np.where(ML>=media-2*stdev)] = 1
    f = f1 + f2    
    ML = ML[np.where(f==2)]
    media = np.median(ML)
    stdev = np.std(ML)        

    N = len(ML)
    f1 = np.zeros(N)
    f2 = np.zeros(N)
    f1[np.where(ML<=media+3*stdev)] = 1
    f2[np.where(ML>=media-3*stdev)] = 1
    f = f1 + f2    
    ML = ML[np.where(f==2)]
    media = np.median(ML)
    stdev = np.std(ML)   

    

    return media, stdev
#################################################################

def myM2L_(L_k, alfa):
      
      
       L = L_k / 1.E10
       if L>=1: 
           return 32.*(L**0.15)
       else:
           return 32.*(L**alfa)
    
def myM2L__(alfa):
      
       L_lst = [1.E13, 4.423E10, 9E8, 1.E7]
       M2L_lst = [myM2L_(1.E13, alfa), myM2L_(4.423E10, alfa), myM2L_(9E8, alfa), myM2L_(1.E7, alfa)]
       
       L_lst = np.asarray(L_lst)
       M2L_lst = np.asarray(M2L_lst)
      
       myM2L = interpolate.interp1d(np.log10(L_lst), np.log10(M2L_lst))
       
       return myM2L 



def L_exp(M, alfa, beta, gama):
    
    M12 = M/1.E12
    
    log10L = log10(alfa) + 10. + beta*log10(M12) + (gama/M12)*log10(exp(1.))
    
    return 10**log10L


def M_exp_f(alfa, beta, gama):
    
    Lbin = np.logspace(5,13,10000)
    Mbin = np.logspace(0.01,50,10000) 
    
    for i in range(len(Mbin)): Lbin[i] = L_exp(Mbin[i], alfa, beta, gama)
    
    M_f = interpolate.interp1d(np.log10(Lbin), np.log10(Mbin))
    
    return M_f
#################################################################

if __name__ == '__main__':
  
  
  fig = plt.figure(figsize=(45/7., 6), dpi=100)
  #ax = fig.add_axes([0.13, 0.1, 0.83,  0.85])  # m-L relation
  ax = fig.add_axes([0.15, 0.13, 0.80,  0.80])  
  
  table = np.genfromtxt('halo.0.100.v01.group' , delimiter='|', filling_values=0, names=True, dtype=None)
  flag    = table['flag']
  grpID   = table['grpID']
  logK    = table['logK']
  Mv_lum  = table['Mv_lum']
  dist    = table['dist']
  Vls     = table['Vls']
  N = len(flag)
    
  table = np.genfromtxt('halos_simple_mock-14.0.csv' , delimiter=',', filling_values=0, names=True, dtype=None)
  grpID2   =  table['grpID']
  HaloMass =  table['HaloMass']
  
  M_halo =  np.zeros(N)
  L_halo =  np.zeros(N)
  for i in range(N):
        
        if flag[i]==2:
            grpID[i] = grpID[i+1]
        
        if flag[i]!=1:  # !=1
            index=np.where(grpID2==grpID[i])[0]
            
            #print i, index, grpID2[index], grpID[i]
            M_halo[i]=HaloMass[index]
            #M_halo[i]=Mv_lum[i]
            L_halo[i]=10**logK[i]
          
  ax.plot(L_halo,M_halo/L_halo, '.', markersize=1, color='black', alpha=0.1)
  
  
  M2L_halo = M_halo/L_halo
  

  alfa =  2
  M1   =  0.5E7
  
  m0 = M1
    
  x = []
  y = []
  y_err = []
    
  while m0<1.E14:
        
        m1 = m0*alfa 
        med, stdev = median(M2L_halo, L_halo, m0, m1)
        y.append(med)
        y_err.append(stdev)
        x.append(sqrt(m0*m1))
        m0 = m1
        if m0>1.E11: alfa=2.5
  
  x = np.asarray(x)
  y = np.asarray(y)
  y_err = np.asarray(y_err)
  #ax.plot(x,y*3, 'g*')
  plt.errorbar(x,y*3.5, yerr=y_err, fmt='*', color='brown')  
  
  
  
  

  
  
  

  Lbin = np.logspace(7,13,200)
  Mbin = np.logspace(10,15,200)
  
  f = myM2L__(-0.9393)
  for i in range(len(Lbin)): Mbin[i] = 10**f(np.log10(Lbin[i]))*Lbin[i]
  ax.plot(Lbin,Mbin/Lbin, 'b.') 


  f = myM2L__(-0.5)
  for i in range(len(Lbin)): Mbin[i] = 10**f(np.log10(Lbin[i]))*Lbin[i]
  ax.plot(Lbin,Mbin/Lbin, 'r-')   
  
  
  for i in range(len(Lbin)): Mbin[i] = myM2L_(Lbin[i], -0.7)*Lbin[i]
  ax.plot(Lbin,Mbin/Lbin, '--', color='orange') 

  
  #alfa = 3.7043   
  #beta = 0.8827   
  #gama = -0.4347
  
  alfa = 4.0 
  beta = 0.82   
  gama = -0.42 
  
  M_f = M_exp_f(alfa, beta, gama)
  for i in range(len(Lbin)): Mbin[i] = 10**M_f(log10(Lbin[i]))
  ax.plot(Lbin,Mbin/Lbin, 'g-')  
  
  
  alfa = 3.25   
  beta = 0.59
  gama = -0.6
  M_f = M_exp_f(alfa, beta, gama)
  for i in range(len(Lbin)): Mbin[i] = 10**M_f(log10(Lbin[i]))
  ax.plot(Lbin,Mbin/Lbin, '--', color='black')  
     
  
  #######################################
  

  
  ax.set_ylabel('M'+r'$_v$'+r'$^{exp}$'+'/L'+r'$_{K_s}$'+' ['+r'$M_\odot/L_\odot$'+']', fontsize=16)
  ax.set_xlabel('K'+r'$_s$'+'-band Luminosity ['+r'$L_\odot$'+']', fontsize=16)
 
  #plt.minorticks_on()
  plt.tick_params(which='major', length=7, width=1.5)
  plt.tick_params(which='minor', length=4, color='#000033', width=1.0)     

  plt.yticks(fontsize=16)
  plt.xticks(fontsize=16)

  #ax.annotate(r'$L^{-0.5}$', (3.E7, 141), rotation=0, color='black', size=18)
  ax.annotate(r'$L^{0.15}$', (4.5E11, 90), rotation=0, color='black', size=18)
  
  ax.annotate(r'$L^{-0.94}$', (2.E8, 1300), rotation=0, color='blue', size=18)
  
  plt.xscale('log')
  plt.yscale('log')
  plt.xlim(1.E7,1.E13)
  plt.ylim(1,1.E4)

  
  
  
  plt.show()
  

