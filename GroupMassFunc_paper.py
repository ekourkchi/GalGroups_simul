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

if __name__ == '__main__':
    
    fig = plt.figure(figsize=(5,6), dpi=100)
    ax = fig.add_subplot()
    
    table = np.genfromtxt('groupmassfunction3_10a.csv' , delimiter=',', filling_values=0, names=True, dtype=None)
    logM = table['logM']
    logN = table['logN']
    Nerr = table['Nerr']
    
    
    table = np.genfromtxt('groupmassfunction3_6.5b.csv' , delimiter=',', filling_values=0, names=True, dtype=None)
    logMb = table['logM']
    logNb = table['logN']
    Nerrb = table['Nerr']   
    
    plt.plot(10**logM, 10**logN, 'b.')
    plt.plot(10**logMb, 5*10**logNb, 'r.')
    
    #plt.title('Groups ... 5-20 Mpc')
    plt.xlabel('Group Mass [M'+r'$_\odot$'+']')
    plt.ylim(0.15,1.E6)
    plt.xlim(1.E9,1.E16)
    plt.ylabel('Number')
    plt.xscale('log')
    plt.yscale('log')
    plt.show()
     


    
