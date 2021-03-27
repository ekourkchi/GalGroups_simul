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


table = np.genfromtxt('M_L_curve_v2.csv' , delimiter=',', filling_values=0, names=True, dtype=None)
Lumin_in   = table['log10_L']
Mass_out   = table['log10_M']
MtoL_mod = interpolate.interp1d(Lumin_in, Mass_out)
#################################################################

# Returns mass in solar unit, given the absolute K-band luminosity
def M2L_paper(L_k):
  
  global MtoL_mod

  if L_k==0:
    return 0
  
  L = L_k / 1.E10
  
  if L <= 0.0927:
    MtoL = 32.0*(L**-0.5)
  elif L > 0.0927 and L < 4.423:
    return 10**MtoL_mod(np.log10(L_k))/L_k
  elif L >= 4.423:
    MtoL = 32*(L**0.15)

  
  return MtoL
#################################################################

######################################################################
def median_error(X):
  
  size = len(X)
  X    = np.sort(X) 
  mean = np.median(X)
  
  u_err = X[int(round(0.84*size))] - mean
  l_err = mean - X[int(round(0.16*size))]
  
  return mean, u_err, l_err

######################################################################

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
    f1[np.where(ML<=media+2.5*stdev)] = 1
    f2[np.where(ML>=media-2.5*stdev)] = 1
    f = f1 + f2    
    ML = ML[np.where(f==2)]
    media = np.median(ML)
    stdev = np.std(ML)        

    N = len(ML)
    f1 = np.zeros(N)
    f2 = np.zeros(N)
    f1[np.where(ML<=media+4*stdev)] = 1
    f2[np.where(ML>=media-4*stdev)] = 1
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
            M_halo[i]=HaloMass[index]
            L_halo[i]=10**logK[i]
          
  #ax.plot(L_halo,M_halo/L_halo, '.', markersize=1, color='brown', alpha=0.05)
  
  
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
  
  #######################################
  ### Plotting the Mock catalog M/L ratio
  
  x = np.asarray(x)  # L
  y = np.asarray(y)  # M/L
  y_err = np.asarray(y_err) 
  
  y = y*3.2
  
  z = x*y   # M
  z_err = x*y_err
  
  plt.fill_betweenx(x, z - z_err, z + z_err, color='brown', alpha=0.2)
  line2, = plt.plot(z,x, '-', color='brown', label='mock catalog') 
  
  #######################################

  
  
  filename = 'set4.all.csv'
  table = np.genfromtxt(filename , delimiter=',', filling_values=0, names=True, dtype=None)
  
  alfa = table['alfa']
  beta = table['beta']
  gama = table['gama']
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
  
  Lbin = np.logspace(7,13,50)
  M_lst = np.logspace(7,13,50)
  M_err_l = np.logspace(7,13,50)
  M_err_u = np.logspace(7,13,50)

  for i in range(len(Lbin)):
      
      Mbin = table['Mbin'+str(i)]
      Mbin = Mbin[np.where(f==6)]
      
      mean, u_err, l_err = median_error(Mbin)
      
      M_lst[i]   = mean
      M_err_l[i] = l_err
      M_err_u[i] = u_err
      

  plt.fill_betweenx(Lbin, M_lst - M_err_l, M_lst + M_err_u, color='green', alpha=0.2)
  line1, = plt.plot(M_lst, Lbin, 'g-', label='best fit')  
  
  
     

  #######################################
  
  LLL = np.logspace(7,13,100)
  MMM = np.logspace(7,13,100)
  for i in range(len(LLL)):
      MMM[i] = LLL[i] * M2L_paper(LLL[i])
      
  line3, = plt.plot(MMM, LLL, '-', color='blue', label='group paper (KT17)') 
  
  
  #######################################  
  plt.legend(handles=[line1, line2, line3], loc=2)
  #plt.legend(handles=[line1, line2], loc=2)
  #######################################  
  

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

  # additional Y-axis (on the right)
  y_ax = ax.twinx()
  y_ax.set_ylim(5.E7,1.E13)
  y_ax.set_yscale('log')
  y_ax.set_yticklabels([])
  y_ax.minorticks_on()
  y_ax.tick_params(which='major', length=5, width=1.5, direction='in')
  y_ax.tick_params(which='minor', length=2, color='#000033', width=1., direction='in')

  # additional X-axis (on the right)
  x_ax = ax.twiny()
  x_ax.set_xlim(1.E10,2.E15)
  x_ax.set_xscale('log')
  x_ax.set_xticklabels([])
  x_ax.minorticks_on()
  x_ax.tick_params(which='major', length=5, width=1.5, direction='in')
  x_ax.tick_params(which='minor', length=2, color='#000033', width=1., direction='in')   
  
  
  plt.show()
