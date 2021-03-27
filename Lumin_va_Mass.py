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
from scipy  import interpolate




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
    
    Lbin = np.logspace(7,13,1000)
    Mbin = np.logspace(10,20,1000) 
    
    for i in range(len(Mbin)): Lbin[i] = L_exp(Mbin[i], alfa, beta, gama)
    
    M_f = interpolate.interp1d(np.log10(Lbin), np.log10(Mbin))
    
    return M_f


#################################################################

if __name__ == '__main__':
  
  fig = plt.figure(figsize=(45/7., 6), dpi=100)
  #ax = fig.add_axes([0.13, 0.1, 0.83,  0.85])  # m-L relation
  ax = fig.add_axes([0.15, 0.13, 0.80,  0.80])  
  
  
  Lbin = np.logspace(7,13,50)
  Mbin = np.logspace(10,15,50)
  
  f = myM2L__(-0.9393)
  for i in range(len(Lbin)): Mbin[i] = 10**f(np.log10(Lbin[i]))*Lbin[i]
  ax.plot(Mbin,Lbin, 'b.') 
  
  
  f = myM2L__(-0.5)
  for i in range(len(Lbin)): Mbin[i] = 10**f(np.log10(Lbin[i]))*Lbin[i]
  ax.plot(Mbin,Lbin, 'r-')   
  
  
  #alfa = 3.7043   
  #beta = 0.8827   
  #gama = -0.4347
  alfa = 4.0 
  beta = 0.82   
  gama = -0.42 
  
  M_f = M_exp_f(alfa, beta, gama)
  for i in range(len(Lbin)): Mbin[i] = 10**M_f(log10(Lbin[i]))
  ax.plot(Mbin,Lbin, 'g-')  
  
  
  alfa = 3.25   
  beta = 0.59
  gama = -0.6
  M_f = M_exp_f(alfa, beta, gama)
  for i in range(len(Lbin)): Mbin[i] = 10**M_f(log10(Lbin[i]))
  ax.plot(Mbin,Lbin, '--', color='black')  
   
  
  #alfa = 4.0 
  #beta = 0.82   
  #gama = -0.42  
  #M_f = M_exp_f(alfa, beta, gama)
  #for i in range(len(Lbin)): Mbin[i] = 10**M_f(log10(Lbin[i]))
  #ax.plot(Mbin,Lbin, 'g--')   
  
  

  ax.set_xlabel('M'+r'$_v$'+' ['+r'$M_\odot$'+']', fontsize=16)
  ax.set_ylabel('K'+r'$_s$'+'-band Luminosity ['+r'$L_\odot$'+']', fontsize=16)
 
  #plt.minorticks_on()
  plt.tick_params(which='major', length=7, width=1.5)
  plt.tick_params(which='minor', length=4, color='#000033', width=1.0)     

  plt.yticks(fontsize=16)
  plt.xticks(fontsize=16)

  ##ax.annotate(r'$L^{-0.5}$', (3.E7, 141), rotation=0, color='black', size=18)
  #ax.annotate(r'$L^{0.15}$', (4.5E11, 35), rotation=0, color='black', size=18)
  
  #ax.annotate(r'$L^{-0.7}$', (2.E8, 800), rotation=0, color='blue', size=18)

  plt.xscale('log')
  plt.yscale('log')
  plt.xlim(1.E10,2.E15)
  plt.ylim(5.E7,1.E13)

  
  
  
  plt.show()
