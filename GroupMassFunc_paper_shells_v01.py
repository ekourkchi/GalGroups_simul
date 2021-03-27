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


#################################################################
#################################################################

# Returns mass in solar unit, given the absolute K-band luminosity
def Mass2(L_k):

  L = L_k / 1.E10
  

  if L <1:
    MtoL = 32.0*(L**-0.60)
  elif L >= 1:
    MtoL = 32*(L**0.15)
  

  
  
  Mass_out = L_k * MtoL
  
  return Mass_out

#################################################################

# Returns mass in solar unit, given the absolute K-band luminosity
def Mass_brent(L_k):
  
  if L_k==0:
    return 0
  
  L = L_k / 1.E10
  
  # ORIGINAL
  if L <= 1:
    MtoL = 32.0*(L**-0.5)     # Original
  else:
    MtoL = 32*(L**0.15)
  
  Mass_out = L_k * MtoL
  
  return Mass_out


#################################################################

def stat(array, mi, mf):
    
    #array = np.asarray(array)
    #N = len(array)
    
    #f1 = np.zeros(N)
    #f2 = np.zeros(N)
    
    #f1[np.where(array>=mi)] = 1
    #f2[np.where(array<mf)] = 1
    #f = f1 + f2
    
    #N = len(f[np.where(f==2)])
    
    N=0
    for a in array:
        if a>=mi and a<mf:
            N+=1
    
    

    return N
#################################################################
def plot_shell(Mv_lum, dist, M_bins, d_min=None, d_max=None, color=None, plot=False):
    

    
    if d_min!=None:
        indices = np.where(dist>d_min)
        Mv_lum_b = Mv_lum[indices]
        dist_b   = dist[indices]  
    else:
        dist_b = dist
        Mv_lum_b = Mv_lum
        
    
    if d_max!=None:
        indices = np.where(dist_b<d_max)
        Mv_lum_b = Mv_lum_b[indices]
        dist_b   = dist_b[indices]      
    
    if d_min==None and d_max==None:
        Mv_lum_b = Mv_lum
        dist_b   = dist
    
    if color==None:
        color='black'
    
    M = M_bins

    n_lst = []
    m_lst = []
    
    for i in range(len(M)-1):
        
        n_lst.append(stat(Mv_lum_b, M[i], M[i+1]))
        m_lst.append(sqrt(M[i]*M[i+1]))
    
    m_lst = np.asarray(m_lst)
    n_lst = np.asarray(n_lst)
    #plt.ploterr(m_lst,1.*n_lst, '.', color=color)
    
    if plot: 
      plt.errorbar(m_lst,1.*n_lst, yerr=np.sqrt(n_lst), fmt='.', color=color)    
    
    return m_lst, n_lst, np.sqrt(n_lst)
#################################################################
def xi2(M_func, m0,n00, n00_err):
    
    xi2_min = 1.E12
    alf = 0
    for alfa in np.arange(0,30,0.1):
        xi2 = 0
        for i in range(len(m0)):
           if m0[i] > 4.5E12 and m0[i] <1.E14:
               xi2 += ((alfa*n00[i]-M_func(m0[i]))/(alfa*n00_err[i]))**2
        if xi2<xi2_min:
            alf = alfa
            xi2_min = xi2
        #print  alfa, xi2
    return alf    # alf = 16. for min xi2


def xi2_tail(M_func, m0,n00, n00_err):
    
    xi2 = 0
    for i in range(len(m0)):
       if m0[i] < 1.0E14 and m0[i] >= 1.0E11:
          xi2 += ((n00[i]-M_func(m0[i]))/(n00_err[i]))**2

    return xi2
##################################################################
def calcMassFunc(Mv_lum, dist, plot=False):
    
    M_bins = np.logspace(10,14,15)
    M_bins = np.append(M_bins, [1.E15])
    Nbins = len(M_bins)

    
    
    m0,n0, nerr0 = plot_shell(Mv_lum, dist, M_bins, d_min=20, d_max=40)
    for i in [0,1,2]: n0[i]=0;nerr0[i]=0
    if plot: plt.errorbar(m0,n0,yerr=nerr0, fmt='.', color='red')
    
    m1,n1, nerr1 = plot_shell(Mv_lum, dist, M_bins, d_max=6)
    if plot: plt.errorbar(m1,n1,yerr=nerr1, fmt='.', color='black')
    
    m2,n2, nerr2 = plot_shell(Mv_lum, dist, M_bins, d_min=6, d_max=10)
    for i in [0,1]: n2[i]=0;nerr2[i]=0
    if plot: plt.errorbar(m2,n2,yerr=nerr2, fmt='.', color='blue')
    
    m3,n3, nerr3 = plot_shell(Mv_lum, dist, M_bins, d_min=10, d_max=20)
    for i in [0,1]: n3[i]=0;nerr3[i]=0
    if plot: plt.errorbar(m3,n3,yerr=nerr3, fmt='.', color='green')
     
    
    N1=0; N0=0
    for i in range(len(n1)):
        if n1[i]!=0 and n0[i]!=0:
            N1+=n1[i]
            N0+=n0[i]
    n11     = np.zeros(len(n1))
    n11_err = np.zeros(len(n1))
    for i in range(len(n1)):
        n11[i] = n1[i]*N0/N1
        n11_err[i] = nerr1[i]*N0/N1
    #plt.errorbar(m1,n11,yerr=n11_err, fmt='.', color='black')

 
    N2=0; N0=0
    for i in range(len(n2)):
        if n2[i]!=0 and n0[i]!=0:
            N2+=n2[i]
            N0+=n0[i]
    n22     = np.zeros(len(n2))
    n22_err = np.zeros(len(n2))
    for i in range(len(n2)):
        n22[i] = n2[i]*N0/N2
        n22_err[i] = nerr2[i]*N0/N2
    #plt.errorbar(m2,n22,yerr=n22_err, fmt='.', color='blue')
 
 
    N3=0; N0=0
    for i in range(len(n3)):
        if n3[i]!=0 and n0[i]!=0:
            N3+=n3[i]
            N0+=n0[i]
    n33     = np.zeros(len(n3))
    n33_err = np.zeros(len(n3))
    for i in range(len(n3)):
        n33[i] = n3[i]*N0/N3
        n33_err[i] = nerr3[i]*N0/N3
    #plt.errorbar(m3,n33,yerr=n33_err, fmt='.', color='green')    
    
    
    n00     = np.zeros(len(n0))
    n00_err = np.zeros(len(n0))
    for i in range(len(n0)):
        sum = 0
        semerr = 0
        if n33[i]!=0:
            sum += n33[i]/n33_err[i]
            semerr += 1./n33_err[i]
        if n22[i]!=0:
            sum += n22[i]/n22_err[i]
            semerr += 1./n22_err[i]
        if n11[i]!=0:
            sum += n11[i]/n11_err[i]
            semerr += 1./n11_err[i]
        if n0[i]!=0:
            sum += n0[i]/nerr0[i]
            semerr += 1./nerr0[i]
        if semerr!=0: 
            n00[i] = sum/semerr
            n00_err[i] = sqrt(2.)/semerr
        else:
            n00[i] = 0    
    
    #plt.errorbar(m0,n00*4,yerr=n00_err*4 ,fmt='.', color='magenta')  
    
    alf = xi2(M_func, m0,n00, n00_err)
    #alf = 16.0
    #print "results: ", alf
    if plot: plt.errorbar(m0,n00*alf,yerr=n00_err*alf , fmt='.', color='magenta')  
    
    return m0, n00*alf, n00_err*alf
      
################################################################
def M2L_init():
    
      L = 1E13
      L_lst = []
      M2L_lst = []
      while L>= 9.E6:
          M = Mass_brent(L)
          M2L_lst.append(M/L)
          L_lst.append(L)
          L = L/1.77827941
      
      return L_lst, M2L_lst
#################################################################
def M2L_mutate(L_lst, M2L_lst):
      
      M2L_lst_new = []
      j = 0
      while(L_lst[j]>=1.E11): j+=1
      
      p = random.randint(j,len(M2L_lst))
      for i in range(len(M2L_lst)):
          alfa = 1.
          if i==p: alfa = random.normalvariate(0, 1.5)
          
          if alfa>=0: M2L_lst_new.append(M2L_lst[i]*alfa)
          if alfa<0: M2L_lst_new.append(M2L_lst[i]/abs(alfa))
      
      return M2L_lst_new  
  
  
#################################################################
def L_015(L_k):
      L = L_k / 1.E10
      return 32*(L**0.15)
  
def MtoLrand():
      
      
      slopes = np.arange(0.15,-1.01,-0.01)
      L = 1E13
      L_lst = []
      M2L_lst = []
      while L>= 9.E6:
          
          if L>= 9.E10:
              M2L_ = L_015(L)
          else:
              p = random.randint(0,116)
              alfa = slopes[p]
              M2L_ = (M2L_/L_**alfa)*L**alfa

          M2L_lst.append(M2L_)
          L_lst.append(L)
          L_ = L
          L = L/sqrt(10)
      
      myM2L = interpolate.interp1d(L_lst, M2L_lst)
      
      return  myM2L
            
#################################################################

if __name__ == '__main__':
    
    ##################################################
    table = np.genfromtxt('MW_17ST_H100.csv' , delimiter=',', filling_values=0, names=True, dtype=None)
    logMass = table['logMass']
    nl20 = table['nl20']
    
    
    ##################################################
    fig = plt.figure(figsize=(5,5), dpi=100)
    ax = fig.add_subplot()
    
    table = np.genfromtxt('groupmassfunction3_10a.csv' , delimiter=',', filling_values=0, names=True, dtype=None)
    logM = table['logM']
    logN = table['logN']
    logNerr = table['Nerr']
    
    
    table = np.genfromtxt('groupmassfunction3_6.5b.csv' , delimiter=',', filling_values=0, names=True, dtype=None)
    logMb = table['logM']
    logNb = table['logN']
    logNerrb = table['Nerr']   
    
    plt.errorbar(10**logM, 10**logN, yerr=10**logNerr , fmt='.', color='black')
    plt.errorbar(10**logMb, 5*10**logNb, yerr=5*10**logNerrb ,  fmt='.', color='black')
    
    #plt.plot(10**logMass, 10**nl20, 'r--')
    plt.plot(10**logMass, 10**(nl20+0.662), '--', color='black')
    
    M_func = interpolate.interp1d(10**logMass, 10**(nl20+0.662))
   
    
    ##################################################
    table = np.genfromtxt('all.iter.2.v44.group' , delimiter='|', filling_values=0, names=True, dtype=None)
    flag   =  table['flag']
    Mv_lum =  table['Mv_lum']
    dcf2 =  table['dcf2']
    mDist =  table['mDist']
    
    Vls =  table['Vls']
    logK = table['logK']
    ######Mv_lum = table['Mv_dyn']
    
    N = len(flag)
    dist = np.zeros(N)

    for i in range(N):
        if flag[i] ==2:
            dist[i] = mDist[i]
        else:
            dist[i] = dcf2[i]

    
    for i in range(N):
       if dist[i]==0:
           dist[i] = Vls[i]/75.
    
    
    indices = np.where(flag!=1)
    Mv_lum = Mv_lum[indices]
    logK = logK[indices]
    dist   = dist[indices] 

    indices = np.where(logK>7)
    Mv_lum = Mv_lum[indices]
    logK = logK[indices]
    dist   = dist[indices] 
    
    N = len(logK)
    
    #for i in range(N):
       #Mv_lum[i] = Mass_brent(10**logK[i])

    
    
    #myM2L = MtoLrand()
    
    myTable = Table()
    empty = []
    myTable.add_column(Column(data=empty,name='n', dtype=np.dtype(int)))
    myTable.add_column(Column(data=empty,name='chi2', format='%2f'))
    
    
    L_lst, M2L_lst = M2L_init()
    myM2L = interpolate.interp1d(L_lst, M2L_lst)
    for i in range(N): Mv_lum[i] = myM2L(10**logK[i])*(10**logK[i])
    m0,n00, n00_err = calcMassFunc(Mv_lum, dist, plot=False)
    chi2 =  xi2_tail(M_func, m0, n00, n00_err) 
    
    for p in range(len(L_lst)):
        myTable.add_column(Column(data=empty,name='M_L'+str(p), format='%0.4f'))
    
    row = [-1, chi2]
    for l in range(len(L_lst)):
             row.append(M2L_lst[l])
    myTable.add_row(row)
    
    
    
    random.seed(40)
    iiter=0  
    while (iiter<500000): 
        
        M2L_lst_new = M2L_mutate(L_lst, M2L_lst)       
        myM2L = interpolate.interp1d(L_lst, M2L_lst_new)
        for i in range(N):
           Mv_lum[i] = myM2L(10**logK[i])*(10**logK[i])
        
                        
        m0,n00, n00_err = calcMassFunc(Mv_lum, dist, plot=False)
        chi2_new =  xi2_tail(M_func, m0, n00, n00_err) 
        
        delta_chi2 = chi2-chi2_new
        if delta_chi2>0: 
          ratio = 1
        else: 
          ratio = exp(0.5*delta_chi2)        
        
        #print chi2, chi2_new
        if random.uniform(0, 1.) < ratio:
            M2L_lst = M2L_lst_new
            chi2 = chi2_new
            
            row = [iiter, chi2]
            for l in range(len(L_lst)):
                row.append(M2L_lst[l])
            myTable.add_row(row)
            
            if iiter%5==0:
                myTable.write('seed40_normal1.5.txt', format='ascii.fixed_width',delimiter=' ', bookend=False)
                
                
            print iiter
            iiter+=1   # loop index
            









    
    ################################################## PLOT
  
    #plt.title('all.iter.2.v44.group + new curved M/L ratio')
    #plt.xlabel('Group Mass [M'+r'$_\odot$'+']')
    #plt.ylim(0.15,1.E6)
    #plt.xlim(1.E9,1.E16)
    #plt.ylabel('Number')
    #plt.xscale('log')
    #plt.yscale('log')
    #plt.show()
    ##################################################


    