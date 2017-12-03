#!/usr/bin/python
__author__ = "Ehsan Kourkchi"
__copyright__ = "Copyright 2017"
__credits__ = ["Ehsan Kourkchi"]
__version__ = "1.0"
__maintainer__ = "Ehsan Kourkchi"
__email__ = "ehsan@ifa.hawaii.edu"
__status__ = "Production"
# As of Nov, 21, 2017

## Importing Important Python Libraries
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

def XYZ(ra, dec, dist):
    
   cl1 = cos(ra*pi/180.)
   sl1 = sin(ra*pi/180.)
   cb1 = cos(dec*pi/180.)
   sb1 = sin(dec*pi/180.)   
   
   X = dist * cl1 * cb1
   Y = dist * sl1 * cb1
   Z = dist * sb1 

   return X, Y, Z

#################################################################

if __name__ == '__main__':
    
    inFile = 'simple_mock-14.0.csv'
    table = np.genfromtxt(inFile , delimiter=',', filling_values=0, names=True, dtype=None)
    
    dist   = table['d_kpc']/1000.
    ra     =  table['ra']
    dec    = table['dec']
    
    for i in range(len(dist)):
        X, Y, Z = XYZ(ra[i], dec[i], dist[i])
        print X, Y, Z
    
    
    
    
