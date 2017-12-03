## Importing Important Python Libraries

from math import *
from matplotlib.widgets import Cursor

import sys
import os
import subprocess
import numpy as np
from datetime import *
import time
import matplotlib

from matplotlib.figure import Figure
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
from pylab import *
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table, Column 
import random
from kapteyn import wcs
from optparse import OptionParser

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
#################################################
def extractPGC_lst(id, grp=False, supergrp=False):
  
  idl = []
  for i in id: idl.append(extractPGC(i, grp=grp, supergrp=supergrp))
  
  return np.asarray(idl)
      
      
  
def extractPGC(id, grp=False, supergrp=False):
  
  if not grp and not supergrp:
    return id
  
  
  if grp:
    pgc = int(id)%100000000
  
  
  if supergrp:
    grp = int(id)%10000000000
    pgc = int(grp)%100000000
  
  return pgc  
########################################################################## 
def create_tables():
  
  #filee_super = 'simul.MagK16.d5.d100.2D.v01.supergroup'
  filee_super = 'halo.0.100.v01.supergroup'
  
  try:
    intable_super = np.genfromtxt(filee_super , delimiter='|', filling_values="-100000", names=True, dtype=None )
  except:
    print "\n[Error] The catalog \""+filee_super+"\" is not available, or has a wrong format ..."
    print >> sys.stderr, "You can use \"python "+sys.argv[0]+" -h\" for help ... (or 'skyplot -h') \n"
    exit(1)  
  

  pgc_super          = intable_super['ID']
  flag_super         = intable_super['flag']
  ra_super           = intable_super['ra']
  dec_super          = intable_super['dec']
  Ks_super           = intable_super['Ks']
  logK_super         = intable_super['logK']
  Vls_super          = intable_super['Vls']
  mDist_super        = intable_super['mDist']
  ID1_super         = intable_super['nest']
  r1t_lum_super      = intable_super['r1t_lum']
  No_Galaxies_super  = intable_super['No_Galaxies']
  Mv_lum_super       = intable_super['Mv_lum']
  
  
  pgc_p_super = [] 
  mass_p_super = [] 
  pgc_tmp_super =  pgc_super[0]
  mass_tmp_super =  Mv_lum_super[0]
  
  for i in range(len(pgc_super)):
      if flag_super[i] == 5 or flag_super[i] == 0 or flag_super[i] == 2:
         pgc_tmp_super  = pgc_super[i]
         mass_tmp_super = Mv_lum_super[i]
         
      if flag_super[i] == 5:
         ID1_super[i] = ID1_super[i+1]
         
      
      pgc_p_super.append(extractPGC(pgc_tmp_super, grp=False, supergrp=True))
      mass_p_super.append(mass_tmp_super)
  
  pgc_p_super = np.asarray(pgc_p_super)    # PGC+   supergroup ID
  mass_p_super = np.asarray(mass_p_super)
  
  
  for i in range(len(Ks_super)): 
      if Ks_super[i] == 0:
          Ks_super[i] = None
          

  for i in range(len(logK_super)): 
      if logK_super[i] == 0:
          logK_super[i] = None
                                  
  for i in range(len(mDist_super)): 
      if mDist_super[i] == 0:
          mDist_super[i] = None
          mDistErr_super[i] = None  
  
    
  indices            = np.where(flag_super==5)
  pgc_p_super_t      = pgc_p_super[indices]
  ID1_super_t        = ID1_super[indices]
  ra_super_t         = ra_super[indices]
  dec_super_t        = dec_super[indices]
  logK_super_t       = logK_super[indices]
  Vls_super_t        = Vls_super[indices]
  mDist_super_t      = mDist_super[indices]
  r1t_lum_super_t    = r1t_lum_super[indices]
  Mv_lum_super_t     = Mv_lum_super[indices]
  
  
  indices            = np.argsort(Mv_lum_super_t)
  indices = indices[::-1]
  pgc_p_super_t      = pgc_p_super_t[indices]
  ID1_super_t        = ID1_super_t[indices]
  ra_super_t         = ra_super_t[indices]
  dec_super_t        = dec_super_t[indices]
  logK_super_t       = logK_super_t[indices]
  Vls_super_t        = Vls_super_t[indices]
  mDist_super_t      = mDist_super_t[indices]
  r1t_lum_super_t    = r1t_lum_super_t[indices]
  Mv_lum_super_t     = Mv_lum_super_t[indices]  
  

  
  myTable = Table()
  myTable.add_column(Column(data=pgc_p_super_t, name='ID1+'))
  myTable.add_column(Column(data=ra_super_t, name='RA', format='%0.4f'))
  myTable.add_column(Column(data=dec_super_t, name='DEC', format='%0.4f'))
  myTable.add_column(Column(data=logK_super_t, name='logK', format='%0.2f'))
  myTable.add_column(Column(data=Vls_super_t, name='Vls', format='%0.0f'))
  myTable.add_column(Column(data=mDist_super_t, name='D', format='%0.2f'))
  myTable.add_column(Column(data=r1t_lum_super_t, name='r1t', format='%0.3f'))
  myTable.add_column(Column(data=np.log10(Mv_lum_super_t), name='Mass_lum', format='%0.3f'))


  #myTable.write('SGroup_table_halo.0.100.v01.csv', format='ascii.fixed_width',delimiter=',', bookend=False)   
  
  
########################################################################## 
########################################################################## 
     
  #filee = 'simul.MagK16.d5.d100.2D.v01.group'
  filee = 'halo.0.100.v01.group'
  
  try:
    intable = np.genfromtxt(filee , delimiter='|', filling_values="-100000", names=True, dtype=None )
  except:
    print "\n[Error] The catalog \""+filee+"\" is not available, or has a wrong format ..."
    print >> sys.stderr, "You can use \"python "+sys.argv[0]+" -h\" for help ... (or 'skyplot -h') \n"
    exit(1)  
  

  pgc          = intable['ID']
  grpID        = intable['grpID']
  flag         = intable['flag'] 
  ra           = intable['ra']
  dec          = intable['dec']
  MagK         = intable['MagK']
  Ks           = intable['Ks']
  logK         = intable['logK']
  Vls          = intable['Vls']
  dist         = intable['dist']
  mDist        = intable['mDist']
  R_2t2        = intable['R2t_lum']
  R_2t2_dyn    = intable['R2t_dyn']
  Rg_dyn       = intable['Rg_dyn']
  ID1          = intable['nest']
  sigmaP_lum   = intable['sigmaP_lum']    
  sigmaP_dyn   = intable['sigmaP_dyn']    
  Mv_lum       = intable['Mv_lum']
  Mv_dyn       = intable['Mv_dyn']
  No_Galaxies  = intable['No_Galaxies']
  
  D_counter = np.zeros((len(pgc),), dtype=np.int)
  i = 0 
  while i < len(pgc):
      
      if flag[i] <= 0:
          if dist[i] != 0:
             D_counter[i]+=1
          i+=1
      
      elif flag[i] == 2:
          j = i + 1
          while flag[j] == 1:
              if dist[j] != 0 : 
                  D_counter[i]+=1
                  D_counter[j]+=1
              j+=1
          i = j    
  

  indices = np.where(flag!=1)
  pgc_grp          = pgc[indices]
  flag_grp         = flag[indices]
  ra_grp           = ra[indices]
  dec_grp          = dec[indices]
  Ks_grp           = Ks[indices]
  logK_grp         = logK[indices]
  Vls_grp          = Vls[indices]
  D_counter_grp    = D_counter[indices]
  dist_grp         = dist[indices]
  mDist_grp        = mDist[indices]
  R_2t2_grp        = R_2t2[indices]
  R_2t2_dyn_grp    = R_2t2_dyn[indices]
  Rg_dyn_grp       = Rg_dyn[indices]
  ID1_grp          = ID1[indices]
  sigmaP_lum_grp   = sigmaP_lum[indices]   
  sigmaP_dyn_grp   = sigmaP_dyn[indices]  
  Mv_lum_grp       = Mv_lum[indices]
  Mv_dyn_grp       = Mv_dyn[indices]
  No_Galaxies_grp  = No_Galaxies[indices] 
  
  indices = np.argsort(logK_grp)
  #indices = indices[::-1]
  pgc_grp          = pgc_grp[indices]
  flag_grp         = flag_grp[indices]
  ra_grp           = ra_grp[indices]
  dec_grp          = dec_grp[indices]
  Ks_grp           = Ks_grp[indices]
  logK_grp         = logK_grp[indices]
  Vls_grp          = Vls_grp[indices]
  D_counter_grp    = D_counter_grp[indices]
  dist_grp         = dist_grp[indices]
  mDist_grp        = mDist_grp[indices]
  R_2t2_grp        = R_2t2_grp[indices]
  R_2t2_dyn_grp    = R_2t2_dyn_grp[indices]
  Rg_dyn_grp       = Rg_dyn_grp[indices]
  ID1_grp          = ID1_grp[indices]
  sigmaP_lum_grp   = sigmaP_lum_grp[indices]   
  sigmaP_dyn_grp   = sigmaP_dyn_grp[indices]  
  Mv_lum_grp       = Mv_lum_grp[indices]
  Mv_dyn_grp       = Mv_dyn_grp[indices]
  No_Galaxies_grp  = No_Galaxies_grp[indices]   
  
  for i in range(len(Ks_grp)): 
      if Ks_grp[i] == 0:
          Ks_grp[i] = None
          

  for i in range(len(logK_grp)): 
      if logK_grp[i] == 0:
          logK_grp[i] = None
                                  
  for i in range(len(mDist_grp)): 
      if mDist_grp[i] == 0:
          mDist_grp[i] = None
  
  


  for i in range(len(Mv_dyn_grp)): 
      if Mv_dyn_grp[i] == 0:
          Mv_dyn_grp[i] = None

  for i in range(len(R_2t2_dyn_grp)): 
      if R_2t2_dyn_grp[i] == 0:
          R_2t2_dyn_grp[i] = None

  for i in range(len(Rg_dyn_grp)): 
      if Rg_dyn_grp[i] == 0:
          Rg_dyn_grp[i] = None


  ID_p  = []
  Mass_p = []
  for pgc1 in ID1_grp:
      i = 0 
      try:
          while pgc1 != ID1_super[i]: i+=1
          ID_p.append(pgc_p_super[i])
          Mass_p.append(mass_p_super[i])
      except:
          print 'Warning ...'
     
        
  ID_p = np.asarray(ID_p)
  Mass_p = np.asarray(Mass_p)
  
  
  indices = np.argsort(Mass_p, kind='mergesort')
  indices = indices[::-1]
  Mass_p           = Mass_p[indices]
  ID_p             = ID_p[indices]
  pgc_grp          = pgc_grp[indices]
  flag_grp         = flag_grp[indices]
  ra_grp           = ra_grp[indices]
  dec_grp          = dec_grp[indices]
  Ks_grp           = Ks_grp[indices]
  logK_grp         = logK_grp[indices]
  Vls_grp          = Vls_grp[indices]
  D_counter_grp    = D_counter_grp[indices]
  dist_grp         = dist_grp[indices]
  mDist_grp        = mDist_grp[indices]
  R_2t2_grp        = R_2t2_grp[indices]
  R_2t2_dyn_grp    = R_2t2_dyn_grp[indices]
  Rg_dyn_grp       = Rg_dyn_grp[indices]
  ID1_grp          = ID1_grp[indices]
  sigmaP_lum_grp   = sigmaP_lum_grp[indices]   
  sigmaP_dyn_grp   = sigmaP_dyn_grp[indices]  
  Mv_lum_grp       = Mv_lum_grp[indices]
  Mv_dyn_grp       = Mv_dyn_grp[indices]
  No_Galaxies_grp  = No_Galaxies_grp[indices]    
  
      
  myTable = Table()
  myTable.add_column(Column(data=ID1_grp, name='ID1'))
  myTable.add_column(Column(data=ID_p, name='ID1+'))
  myTable.add_column(Column(data=No_Galaxies_grp, name='Mem'))
  myTable.add_column(Column(data=ra_grp, name='Gra', format='%0.4f'))
  myTable.add_column(Column(data=dec_grp, name='Gdec', format='%0.4f'))
  myTable.add_column(Column(data=Ks_grp, name='Ks', format='%0.2f'))
  myTable.add_column(Column(data=logK_grp, name='logK', format='%0.2f'))
  myTable.add_column(Column(data=Vls_grp, name='Vls', format='%0.0f'))
  myTable.add_column(Column(data=D_counter_grp, name='Mem_D'))
  myTable.add_column(Column(data=mDist_grp, name='D', format='%0.2f'))
  myTable.add_column(Column(data=sigmaP_lum_grp, name='Sigma_L', format='%0.0f'))
  myTable.add_column(Column(data=sigmaP_dyn_grp, name='Sigma_V', format='%0.0f'))
  myTable.add_column(Column(data=R_2t2_grp, name='R2t_lum', format='%0.3f'))
  myTable.add_column(Column(data=Rg_dyn_grp, name='Rg_dyn', format='%0.3f'))
  myTable.add_column(Column(data=np.log10(Mv_lum_grp), name='Mass_lum', format='%0.3f'))
  myTable.add_column(Column(data=np.log10(Mv_dyn_grp), name='Mass_dyn', format='%0.3f'))

  #myTable.write('Group_table_halo.0.100.v01.csv', format='ascii.fixed_width',delimiter=',', bookend=False)
  
    #########################################################################

  
  indices          = np.where(flag<2)
  pgc_gal          = pgc[indices]
  grpID_gal        = grpID[indices]
  flag_gal         = flag[indices]
  ra_gal           = ra[indices]
  dec_gal          = dec[indices]
  Ks_gal           = Ks[indices]
  MagK_gal         = MagK[indices]
  logK_gal         = logK[indices]
  Vls_gal          = Vls[indices]
  dist_gal         = dist[indices]
  mDist_gal        = mDist[indices]
  R_2t2_gal        = R_2t2[indices]
  ID1_gal         = ID1[indices]
  sigmaP_lum_gal   = sigmaP_lum[indices]   
  sigmaP_dyn_gal   = sigmaP_dyn[indices]  
  Mv_lum_gal       = Mv_lum[indices]
  No_Galaxies_gal  = No_Galaxies[indices] 

  
  for i in range(len(MagK_gal)): 
      if MagK_gal[i] == -100000.0:
          MagK_gal[i] = None          
   
  for i in range(len(Ks_gal)): 
      if Ks_gal[i] == 0:
          Ks_gal[i] = None
          

  for i in range(len(logK_gal)): 
      if logK_gal[i] == 0:
          logK_gal[i] = None
                                  
  for i in range(len(dist_gal)): 
      if dist_gal[i] == 0:
          dist_gal[i] = None
  

  #########################################################################
  
  
  myTable = Table()
  
  myTable.add_column(Column(data=pgc_gal, name='ID'))
  myTable.add_column(Column(data=grpID_gal, name='grpID'))
  myTable.add_column(Column(data=ra_gal, name='RA', format='%0.4f'))
  myTable.add_column(Column(data=dec_gal, name='DE', format='%0.4f'))
  myTable.add_column(Column(data=MagK_gal, name='MagK', format='%0.2f'))
  myTable.add_column(Column(data=Ks_gal, name='K_t', format='%0.2f'))
  myTable.add_column(Column(data=logK_gal, name='lgL_K', format='%0.2f'))
  myTable.add_column(Column(data=Vls_gal, name='Vls', format='%0.0f'))
  myTable.add_column(Column(data=dist_gal, name='D_i', format='%0.2f'))
  myTable.add_column(Column(data=ID1_gal, name='ID1'))

  #myTable.write('Gal_table_halo.0.100.v01.csv', format='ascii.fixed_width',delimiter=',', bookend=False)


  EQX = np.zeros(len(pgc_gal))
  EQY = np.zeros(len(pgc_gal))
  EQZ = np.zeros(len(pgc_gal))
  for i in range(len(pgc_gal)):
      EQX[i], EQY[i], EQZ[i] = XYZ(ra_gal[i], dec_gal[i], dist_gal[i])

  
  ID1_galgrp = []
  ID_p_gal = []
  No_Galaxies_galgrp = []
  ra_galgrp = []
  dec_galgrp = []
  Ks_galgrp = []
  logK_galgrp = []
  Vls_galgrp = []
  D_counter_galgrp = []
  mDist_galgrp = []
  sigmaP_lum_galgrp = []
  sigmaP_dyn_galgrp = []
  R_2t2_galgrp = []
  Rg_dyn_galgrp = []
  Mv_lum_galgrp = []
  Mv_dyn_galgrp = []
  
  for pgc1 in ID1_gal:
      i = 0 
      try:
          while ID1_grp[i] != pgc1: i+=1
          ID1_galgrp.append(ID1_grp[i])
          ID_p_gal.append(ID_p[i])
          No_Galaxies_galgrp.append(No_Galaxies_grp[i])
          ra_galgrp.append(ra_grp[i])
          dec_galgrp.append(dec_grp[i])
          Ks_galgrp.append(Ks_grp[i])
          logK_galgrp.append(logK_grp[i])
          Vls_galgrp.append(Vls_grp[i])
          D_counter_galgrp.append(D_counter_grp[i])
          mDist_galgrp.append(mDist_grp[i])
          sigmaP_lum_galgrp.append(sigmaP_lum_grp[i])
          sigmaP_dyn_galgrp.append(sigmaP_dyn_grp[i])
          R_2t2_galgrp.append(R_2t2_grp[i])
          Rg_dyn_galgrp.append(Rg_dyn_grp[i])
          Mv_lum_galgrp.append(Mv_lum_grp[i])
          Mv_dyn_galgrp.append(Mv_dyn_grp[i])      
      except:
          print 'Warning ...'

  ID1_galgrp = np.asarray(ID1_galgrp)
  ID_p_gal = np.asarray(ID_p_gal)
  No_Galaxies_galgrp = np.asarray(No_Galaxies_galgrp)
  ra_galgrp = np.asarray(ra_galgrp)
  dec_galgrp = np.asarray(dec_galgrp)
  Ks_galgrp = np.asarray(Ks_galgrp)
  logK_galgrp = np.asarray(logK_galgrp)
  Vls_galgrp = np.asarray(Vls_galgrp)
  D_counter_galgrp = np.asarray(D_counter_galgrp)
  mDist_galgrp = np.asarray(mDist_galgrp)
  sigmaP_lum_galgrp = np.asarray(sigmaP_lum_galgrp)
  sigmaP_dyn_galgrp = np.asarray(sigmaP_dyn_galgrp)
  R_2t2_galgrp = np.asarray(R_2t2_galgrp)
  Rg_dyn_galgrp = np.asarray(Rg_dyn_galgrp)
  Mv_lum_galgrp = np.asarray(Mv_lum_galgrp)
  Mv_dyn_galgrp = np.asarray(Mv_dyn_galgrp)

  
  #myTable.add_column(Column(data=ID1_galgrp, name='g_ID1'))
  myTable.add_column(Column(data=ID_p_gal, name='ID1+'))
  myTable.add_column(Column(data=No_Galaxies_galgrp, name='Ng'))
  myTable.add_column(Column(data=ra_galgrp, name='gRA', format='%0.4f'))
  myTable.add_column(Column(data=dec_galgrp, name='gDE', format='%0.4f')) 
  myTable.add_column(Column(data=Ks_galgrp, name='gK_t', format='%0.2f'))
  myTable.add_column(Column(data=logK_galgrp, name='glgL_K', format='%0.2f'))
  myTable.add_column(Column(data=Vls_galgrp, name='gVls', format='%0.0f'))
  myTable.add_column(Column(data=D_counter_galgrp, name='gND'))
  myTable.add_column(Column(data=mDist_galgrp, name='gD', format='%0.2f'))
  myTable.add_column(Column(data=sigmaP_lum_galgrp, name='gsigL', format='%0.0f'))
  myTable.add_column(Column(data=sigmaP_dyn_galgrp, name='gsigV', format='%0.0f'))
  myTable.add_column(Column(data=R_2t2_galgrp, name='gR2t', format='%0.3f'))
  myTable.add_column(Column(data=Rg_dyn_galgrp, name='gRdyn', format='%0.3f'))
  myTable.add_column(Column(data=np.log10(Mv_lum_galgrp), name='gMassL', format='%0.3f'))
  myTable.add_column(Column(data=np.log10(Mv_dyn_galgrp), name='gMassdyn', format='%0.3f'))

  myTable.add_column(Column(data=EQX, name='EQX', format='%0.2f'))
  myTable.add_column(Column(data=EQY, name='EQY', format='%0.2f'))
  myTable.add_column(Column(data=EQZ, name='EQZ', format='%0.2f'))

  #myTable.write('EDD_table_simul.MagK16.d5.d100.2D.v01.csv', format='ascii.fixed_width',delimiter='|', bookend=False)    
  
  myTable.write('EDD_table_halo.0.100.v01.csv', format='ascii.fixed_width',delimiter='|', bookend=False)
#########################################################################

########################################################################## 
  
if __name__ == '__main__':


  create_tables()
  
  
  
  
