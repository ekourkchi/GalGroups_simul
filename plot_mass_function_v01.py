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
    
    
    #simul.3D.iter.0.v01.group
    #all.iter.2.v43.group
    table = np.genfromtxt('all.iter.2.v44.group' , delimiter='|', filling_values=0, names=True, dtype=None)
    flag   =  table['flag']
    Mv_lum =  table['Mv_lum']
    dist =  table['dcf2']
    Vls =  table['Vls']
    logK = table['logK']
    
    N = len(flag)
    
    for i in range(N):
       if dist[i]==0:
           dist[i] = Vls[i]/75.
    
    

    indices = np.where(flag!=1)
    Mv_lum = Mv_lum[indices]
    dist   = dist[indices]
    logK = logK[indices]
    
    indices = np.where(dist>=5)
    Mv_lum = Mv_lum[indices]
    dist   = dist[indices]    
    logK = logK[indices]
    
    indices = np.where(dist<20)
    Mv_lum = Mv_lum[indices]
    dist   = dist[indices]      
    logK = logK[indices]
    
    
    
    print min(dist), max(dist)
    print min(Mv_lum), max(Mv_lum)
    
    
    
    alfa = 3 # 2 # 3
    M1 = 1.E7# 1.E10 # 1.E7
    #plt.xlabel('Mass [M'+r'$_\odot$'+']')
    plt.xlabel('Luminosity [L'+r'$_{k\odot}$]')
    
    m0 = M1
    
    n_lst = []
    m_lst = []
    
    while m0<1.E14:
        
        m1 = m0*alfa
        #n_lst.append(stat(Mv_lum, m0, m1))
        n_lst.append(stat(10**logK, m0, m1))
        m_lst.append(sqrt(m0*m1))
        m0 = m1
    
    m_lst = np.asarray(m_lst)
    n_lst = np.asarray(n_lst)
    #plt.plot(m_lst,1.*n_lst, 'ro')
    
    
    

    table = np.genfromtxt('halo.0.100.v01.group' , delimiter='|', filling_values=0, names=True, dtype=None)
    flag   =  table['flag']
    grpID   =  table['grpID']
    Mv_lum =  table['Mv_lum']
    logK = table['logK']
    dist =  table['dist']
    Vls =  table['Vls']
    N = len(flag)
    
    table = np.genfromtxt('halos_simple_mock-14.0.csv' , delimiter=',', filling_values=0, names=True, dtype=None)
    grpID2   =  table['grpID']
    HaloMass =  table['HaloMass']
    
    M_halo =  np.zeros(N)
    for i in range(N):
        
        if flag[i]==2:
            grpID[i] = grpID[i+1]
        
        if flag[i]!=1:
            index=np.where(grpID2==grpID[i])[0]
            
            #print i, index, grpID2[index], grpID[i]
            M_halo[i]=HaloMass[index]
        
    


    indices = np.where(flag!=1)
    Mv_lum = Mv_lum[indices]
    dist   = dist[indices]
    M_halo = M_halo[indices]
    logK = logK[indices]
    
    
    indices = np.where(dist>=2)
    Mv_lum = Mv_lum[indices]
    dist   = dist[indices]  
    M_halo = M_halo[indices]
    logK = logK[indices]
    
    indices = np.where(dist<15)
    Mv_lum = Mv_lum[indices]
    dist   = dist[indices]
    M_halo = M_halo[indices]
    logK = logK[indices]
    

    
    m0 = M1
    
    n_lst = []
    m_lst = []
    
    while m0<1.E14:
        
        m1 = m0*alfa
        #n_lst.append(stat(Mv_lum, m0, m1))
        #n_lst.append(stat(M_halo, m0, m1))
        n_lst.append(stat(10**logK, m0, m1))
        m_lst.append(sqrt(m0*m1))
        m0 = m1
    
    
    m_lst = np.asarray(m_lst)
    n_lst = np.asarray(n_lst)
    plt.plot(m_lst,1.*n_lst/sum(n_lst), 'bo') 
    print 'N total: ', sum(n_lst) 
    
    
    
    L = np.logspace(7,13,20)
    N1 = np.zeros(len(L))
    N2 = np.zeros(len(L))
    for i in range(len(L)):
        N1[i] = luminFunc(L[i])
        N2[i] = luminFunc2(L[i])
    
    plt.plot(L,N1/sum(n_lst), 'r.')
    plt.plot(L,N2/sum(n_lst), 'g-')
    
    
    
    
    
    
    #table = np.genfromtxt('simul.25.45.2D.iter.0.v01.group' , delimiter='|', filling_values=0, names=True, dtype=None)
    #flag   =  table['flag']
    #Mv_lum =  table['Mv_lum']
    #dist =  table['dist']
    #Vls =  table['Vls']
    
    #N = len(flag)

    #indices = np.where(flag!=1)
    #Mv_lum = Mv_lum[indices]
    #dist   = dist[indices]
    
    #indices = np.where(dist>=20)
    #Mv_lum = Mv_lum[indices]
    #dist   = dist[indices]    
    
    #indices = np.where(dist<30)
    #Mv_lum = Mv_lum[indices]
    #dist   = dist[indices]      
    

    
    #alfa = 1.5
    #m0 = 1.E10
    
    #n_lst = []
    #m_lst = []
    
    #while m0<1.E14:
        
        #m1 = m0*alfa
        #n_lst.append(stat(Mv_lum, m0, m1))
        #m_lst.append(sqrt(m0*m1))
        #m0 = m1
    
    
    #m_lst = np.asarray(m_lst)
    #n_lst = np.asarray(n_lst)
    ##plt.plot(m_lst,n_lst, 'g.')       
    
    
    plt.title('Groups ... 5-20 Mpc')
    
    plt.ylim(0.0001,0.3)
    plt.ylabel('Number')
    plt.xscale('log')
    plt.yscale('log')
    plt.show()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
