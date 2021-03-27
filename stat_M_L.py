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
    
  Lbin = np.logspace(7,13,50)
  Mbin = np.logspace(10,15,50)
  
  
  
  filename = 'set4.csv'
  table = np.genfromtxt(filename , delimiter=',', filling_values=0, names=True, dtype=None)
  alfa = table['alfa']
  beta = table['beta']
  gama = table['gama']
  
  
  myTable = Table()
  empty = []
  myTable.add_column(Column(data=empty,name='n', dtype=np.dtype(int)))
  myTable.add_column(Column(data=empty,name='alfa', format='%0.4f'))
  myTable.add_column(Column(data=empty,name='beta', format='%0.4f'))
  myTable.add_column(Column(data=empty,name='gama', format='%0.4f'))  
  for i in range(len(Lbin)): 
      myTable.add_column(Column(data=empty,name='Mbin'+str(i))) 
  
  
  for i in range(len(alfa)):
      M_f = M_exp_f(alfa[i], beta[i], gama[i])
      for l in range(len(Lbin)): Mbin[l] = 10**M_f(log10(Lbin[l]))
      row = [i, alfa[i], beta[i], gama[i]]
      for l in range(len(Lbin)): row.append(Mbin[l])
      myTable.add_row(row)
      print i
             
             
  myTable.write('set4.all.csv', format='ascii.fixed_width',delimiter=',', bookend=False)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
