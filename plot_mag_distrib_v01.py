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


#################################################################

if __name__ == '__main__':
    
    fig = plt.figure(figsize=(5,6), dpi=100)
    ax = fig.add_subplot()
    
    
    #simul.3D.iter.0.v01.group
    #all.iter.2.v43.group
    table = np.genfromtxt('all.iter.2.v44.group' , delimiter='|', filling_values=0, names=True, dtype=None)
    flag   =  table['flag']
    dist =  table['dcf2']
    Vls =  table['Vls']
    Ks =  table['Ks']
    
    N = len(flag)
    
    for i in range(N):
       if dist[i]==0:
           dist[i] = Vls[i]/75.
    
    MagK = Ks - 5*np.log10(dist) - 30 + 5

    indices = np.where(flag!=2)
    MagK = MagK[indices]
    dist   = dist[indices]
    
    indices = np.where(dist>=5)
    MagK = MagK[indices]
    dist   = dist[indices]    
    
    indices = np.where(dist<25)
    MagK = MagK[indices]
    dist   = dist[indices]      
    
    
    print min(MagK), max(MagK)

    
    
    
    alfa = 1
    m0 = -27
    
    n_lst = []
    m_lst = []
    
    while m0<-13:
        
        m1 = m0+alfa
        n_lst.append(stat(MagK, m0, m1))
        m_lst.append(0.5*(m0+m1))
        m0 = m1
    
    m_lst = np.asarray(m_lst)
    n_lst = np.asarray(n_lst)
    plt.plot(m_lst,n_lst, 'r-')
    
    
    

    table = np.genfromtxt('halo.0.100.v01.group' , delimiter='|', filling_values=0, names=True, dtype=None)
    flag   =  table['flag']
    MagK =  table['MagK']
    dist =  table['dist']
    Vls =  table['Vls']
    
    N = len(flag)

    indices = np.where(flag!=2)
    MagK = MagK[indices]
    dist   = dist[indices]
    
    indices = np.where(dist>=5)
    MagK = MagK[indices]
    dist   = dist[indices]    
    
    indices = np.where(dist<25)
    MagK = MagK[indices]
    dist   = dist[indices]      
    

    
    alfa = 1
    m0 = -27
    
    n_lst = []
    m_lst = []
    
    while m0<-13:
        
        m1 = m0+alfa
        n_lst.append(stat(MagK, m0, m1))
        m_lst.append(0.5*(m0+m1))
        m0 = m1
    
    
    m_lst = np.asarray(m_lst)
    n_lst = np.asarray(n_lst)
    plt.plot(m_lst,n_lst, 'b-')    
    
    
   
    
    
    plt.title('5-25 Mpc')
    plt.xlabel('Mag K [mag]')
    plt.ylabel('Number')
    #plt.xscale('log')
    plt.yscale('log')
    plt.show()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
