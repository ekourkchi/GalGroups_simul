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

#################################################################
def stat(array, mi, mf):
    
    array = np.asarray(array)
    N = len(array)
    
    f1 = np.zeros(N)
    f2 = np.zeros(N)
    
    f1[np.where(array>=mi)] = 1
    f2[np.where(array<mf)] = 1
    f = f1 + f2

    return len(array[np.where(f==2)])


#########
#################################################################
def luminFunc(l):
    
    l1, n1 = 3.E12, 1.
    l2, n2 = 4.09E10, 95.
    l3, n3 = 1.2E7, 400.
    
    # right side
    if l>l2:
        alfa12 = log10(n2/n1)/log10(l2/l1)
        beta12 = log10(n2)-alfa12*log10(l2)
        return 10**(alfa12*log10(l)+beta12)
    
    # left side 
    else:
        alfa23 = log10(n3/n2)/log10(l3/l2)
        beta23 = log10(n3)-alfa23*log10(l3)
        return 10**(alfa23*log10(l)+beta23)   



def luminFunc2(l):
    
    l1, n1 = 3.E12, 1.
    l2, n2 = 1.156E11, 28.
    l3, n3 = 3.398E9, 150.
    l4, n4 = 1.2E7, 400.
    
    # right side
    if l>l2:
        alfa12 = log10(n2/n1)/log10(l2/l1)
        beta12 = log10(n2)-alfa12*log10(l2)
        return 10**(alfa12*log10(l)+beta12)
    
    # middle side 
    elif l>l3 and l<=l2:
        alfa23 = log10(n3/n2)/log10(l3/l2)
        beta23 = log10(n3)-alfa23*log10(l3)
        return 10**(alfa23*log10(l)+beta23)  
    
    # left side
    else:
        alfa34 = log10(n4/n3)/log10(l4/l3)
        beta34 = log10(n4)-alfa34*log10(l4)
        return 10**(alfa34*log10(l)+beta34)          
        
    
#################################################################

if __name__ == '__main__':
    
    fig = plt.figure(figsize=(5,6), dpi=100)
    ax = fig.add_subplot()
    
    # N total:  1594
    random.seed(1)
    
    
    
    alfa = 3 # 2 # 3
    M1 =  1.E7 # 1.E10 # 1.E7
    #plt.xlabel('Mass [M'+r'$_\odot$'+']')
    plt.xlabel('Luminosity [L'+r'$_{k\odot}$]')
    
    
      
    

    N = 200000
    L = np.zeros(N)
    M = np.zeros(N)
    i = 0
    while i<N:
        
        power = random.uniform(6,12)
        coeff = random.uniform(0,10)
        l = coeff * (10**power)
        
        Number = random.uniform(0,500)
        if Number < luminFunc2(l):
            L[i] = l
            M[i] = Mass(l)
            i+=1
        


    m0 = M1
    n_lst = []
    m_lst = []
    while m0<1.E14:
        
        m1 = m0*alfa
        n_lst.append(stat(L, m0, m1))
        #n_lst.append(stat(M, m0, m1))
        m_lst.append(sqrt(m0*m1))
        m0 = m1
    
    m_lst = np.asarray(m_lst)
    n_lst = np.asarray(n_lst)
    plt.plot(m_lst,1.*n_lst/N, 'g-')

    



    N = 200000
    L = np.zeros(N)
    M = np.zeros(N)
    i = 0
    while i<N:
        
        power = random.uniform(6,12)
        coeff = random.uniform(0,10)
        l = coeff * (10**power)
        
        Number = random.uniform(0,500)
        if Number < luminFunc(l):
            L[i] = l
            M[i] = Mass(l)
            i+=1
        


    
    
    m0 = M1
    n_lst = []
    m_lst = []
    while m0<1.E14:
        
        m1 = m0*alfa
        n_lst.append(stat(L, m0, m1))
        #n_lst.append(stat(M, m0, m1))
        m_lst.append(sqrt(m0*m1))
        m0 = m1
    
    m_lst = np.asarray(m_lst)
    n_lst = np.asarray(n_lst)
    plt.plot(m_lst,1.*n_lst/N, 'r.')
    
    
    plt.ylim(0.0001,0.3)
    plt.ylabel('Number')
    plt.xscale('log')
    plt.yscale('log')
    plt.show()    
    
    
    
    
    
    
    
    
    
    
    
