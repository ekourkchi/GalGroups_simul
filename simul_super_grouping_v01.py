#!/usr/bin/python
__author__ = "Ehsan Kourkchi"
__copyright__ = "Copyright 2016"
__credits__ = ["Ehsan Kourkchi"]
__version__ = "1.0"
__maintainer__ = "Ehsan Kourkchi"
__email__ = "ehsan@ifa.hawaii.edu"
__status__ = "Production"
# As of April, 8, 2016
###
# Written by Ehsan Kourkchi (September 2015)
# email: ehsan@ifa.hawaii.edu
# This code, identifies groups of galaxies, given a 
# a cataloge of galaxies. The radial velocity of galaxies would 
# be used to find galaxies with relatively the same radial velocities
# and almost the same position on the sky
###

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



# **************************************
# Global Variables
# Physical Constants
# **************************************
H0 = 75.           # hubble constant
raV = 102.8806    # M87 center - super galactic longitude
decV = -2.3479     # M87 center - super galactic latitude
G = 6.67384E-11   # Gravitational constant
M_sun = 1.9891E30  # Solar maas [Kg]
Mpc_km = 3.08567758E19
#t0_age = Mpc_km / H0   # Universe Age [sec]
Yr= 365.25 * 24 * 3600    # [sec]
t0_age = 13.8E9 * Yr # Universe Age [sec] - 13.8 Gyr
h = 0.75  # hubble constant
# **************************************
LtoM_func = None

# Returns mass in solar unit, given the absolute K-band luminosity
def Mass(L_k):

  L = L_k / 1.E10
  
  if L <= 0.0927:
    MtoL = 32.0*(L**-0.5)
  elif L > 0.0927 and L < 4.423:
    #MtoL = 58.0*(L**-0.25)
    return LtoM_func(L_k)
  elif L >= 4.423:
    MtoL = 32*(L**0.15)
  
  Mass_out = L_k * MtoL
  
  return Mass_out



# **************************************
def extractPGC(id, grp=False, supergrp=False):
  
  if not grp and not supergrp:
    return id
  
  
  if grp:
    pgc = int(id)%100000000
  
  
  if supergrp:
    grp = int(id)%10000000000
    pgc = int(grp)%100000000
  
  return pgc

# **************************************
# returns angular separation of 
# two vectors in radian
# **************************************
def angle(l1, b1, l2, b2):
  
   cl1 = cos(l1*pi/180.)
   sl1 = sin(l1*pi/180.)
   cb1 = cos(b1*pi/180.)
   sb1 = sin(b1*pi/180.)
   
   x1 = cl1 * cb1
   y1 = sl1 * cb1
   z1 = sb1
   
   cl2 = cos(l2*pi/180.)
   sl2 = sin(l2*pi/180.)
   cb2 = cos(b2*pi/180.)
   sb2 = sin(b2*pi/180.)
   
   x2 = cl2 * cb2
   y2 = sl2 * cb2
   z2 = sb2   
   
   XdotY = x1*x2 + y1*y2 + z1*z2
   #X2 = sqrt(x1**2 + y1**2 + z1**2)
   #Y2 = sqrt(x2**2 + y2**2 + z2**2)
   
   if XdotY > 1 :
     theta12 = 0.
   elif XdotY < -1 :
     theta12 = -1.*pi
   else:
     theta12 = acos(XdotY)  
   return theta12   # radian
# **************************************
# returns sign of a number 
def sign(x):
  if x<0 : return -1
  if x>=0 : return 1
# **************************************
# L is luminosity
# l,b are galaxtic coordinates, and d is distance
# returns the barycentric coordiantes of a pair 
def barycenter(L1, l1, b1, d1, L2, l2, b2, d2):
  
   
   cl1 = cos(l1*pi/180.)
   sl1 = sin(l1*pi/180.)
   cb1 = cos(b1*pi/180.)
   sb1 = sin(b1*pi/180.)
   
   x1 = d1 * cl1 * cb1
   y1 = d1 * sl1 * cb1
   z1 = d1 * sb1
   
   cl2 = cos(l2*pi/180.)
   sl2 = sin(l2*pi/180.)
   cb2 = cos(b2*pi/180.)
   sb2 = sin(b2*pi/180.)
   
   x2 = d2 * cl2 * cb2
   y2 = d2 * sl2 * cb2
   z2 = d2 * sb2   


   L_tot = L1 + L2
   x = (L1*x1+L2*x2)/L_tot
   y = (L1*y1+L2*y2)/L_tot
   z = (L1*z1+L2*z2)/L_tot


   r = sqrt(x**2+y**2)
   b = atan(z/r)*180/pi
   l = atan(y/x)*180/pi
   

   if sign(x) < 0 : l+=180
   if sign(x) > 0 and sign(y) < 0: l+=360   
   
   if d1==0 or d2==0:
      d = 0.
   else:
      d = sqrt(x**2+y**2+z**2)


   
   return l, b, d

# **************************************
def touch(gr1, gr2, coeff=1.0, error='None', projected=0):
   
   l1  = gr1.ra
   b1  = gr1.dec
   d1    = gr1.dist
   r1t_1 = gr1.r1t
   
   l2  = gr2.ra
   b2  = gr2.dec
   d2    = gr2.dist
   r1t_2 = gr2.r1t
   
   if angle(l1, b1, l2, b2)*180./pi > 30: return False
   
   
   if gr1.mDist <0.001 or gr2.mDist <0.001:
     v1 = gr1.Vls
     v2 = gr2.Vls
     if v1<0: v1=100
     if v2<0: v2=100
     if abs(v1-v2) > 0.3*min(v1,v2):
        return False
   
   if gr1.mDist==0 and gr1.Vls==0:
     return False
   if gr2.mDist==0 and gr2.Vls==0:
     return False   
   
   
   if projected==1:
     d = 0.5*(d1+d2)
     d1 = d
     d2 = d
   
   cl1 = cos(l1*pi/180.)
   sl1 = sin(l1*pi/180.)
   cb1 = cos(b1*pi/180.)
   sb1 = sin(b1*pi/180.)   
   
   x1 = d1 * cl1 * cb1
   y1 = d1 * sl1 * cb1
   z1 = d1 * sb1 
   
   
   cl2 = cos(l2*pi/180.)
   sl2 = sin(l2*pi/180.)
   cb2 = cos(b2*pi/180.)
   sb2 = sin(b2*pi/180.)
   
   x2 = d2 * cl2 * cb2
   y2 = d2 * sl2 * cb2
   z2 = d2 * sb2    
   
   R12 = sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
   
   if R12 < coeff*(r1t_1 + r1t_2):
     return True
 #############
   if error=='Yes':
     
     d21 = 0
     d22 = 0
     if gr2.mDistErr>0 and gr2.mDist>1:
       d21 = gr2.mDist-gr2.mDistErr*gr2.mDist
       d22 = gr2.mDist+gr2.mDistErr*gr2.mDist
     
     if d21>0:
      x2 = d21 * cl2 * cb2
      y2 = d21 * sl2 * cb2
      z2 = d21 * sb2    
      
      R12 = sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
      
      if R12 < coeff*(r1t_1 + r1t_2):
	return True       
      
     if d22>0:
      x2 = d22 * cl2 * cb2
      y2 = d22 * sl2 * cb2
      z2 = d22 * sb2    
      
      R12 = sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
      
      if R12 < coeff*(r1t_1 + r1t_2):
	return True        
 #############
   
   
   return False
   
# **************************************
   
def distance(gr1, gr2, coeff=1.0):
   
   l1  = gr1.ra
   b1  = gr1.dec
   d1    = gr1.dist
   r1t_1 = gr1.r1t
   
   l2  = gr2.ra
   b2  = gr2.dec
   d2    = gr2.dist
   r1t_2 = gr2.r1t
   
   if angle(l1, b1, l2, b2)*180./pi > 30: return None
   
   #if gr1.mDist == 0 or gr2.mDist == 0:
     #d1 = gr1.Vls / H0
     #if d1<1. : d1=1
   
     #d2 = gr2.Vls / H0
     #if d2<1. : d2=1     
     
   if gr1.mDist <0.001 or gr1.mDist <0.001:
     v1 = gr1.Vls
     v2 = gr2.Vls
     if v1<0: v1=100
     if v2<0: v2=100
     if abs(v1-v2) > 0.3*min(v1,v2):
        return None

   
   cl1 = cos(l1*pi/180.)
   sl1 = sin(l1*pi/180.)
   cb1 = cos(b1*pi/180.)
   sb1 = sin(b1*pi/180.)   
   
   x1 = d1 * cl1 * cb1
   y1 = d1 * sl1 * cb1
   z1 = d1 * sb1 
   
   
   cl2 = cos(l2*pi/180.)
   sl2 = sin(l2*pi/180.)
   cb2 = cos(b2*pi/180.)
   sb2 = sin(b2*pi/180.)
   
   x2 = d2 * cl2 * cb2
   y2 = d2 * sl2 * cb2
   z2 = d2 * sb2    
   
   R12 = sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
   
   if R12 < coeff*(r1t_1 + r1t_2):
     return R12
   else: 
     return None  

# **************************************
def readgrouplist(filename, RA=[0,360], DEC=[-90,90]):
  
   
   try:
      mytable = np.genfromtxt(filename , delimiter='|', filling_values="-100000", names=True, dtype=None )
   
   except:
      print "\n[Error] The catalog \""+filee+"\" is not available, or has a wrong format ..."
      print >> sys.stderr, "You can use \"python "+sys.argv[0]+" -h\" for help ... (or 'skyplot -h') \n"
      exit(1)
  
   id          = mytable['ID']
   flag        = mytable['flag']
   ra          = mytable['ra']
   dec         = mytable['dec']
   Ks          = mytable['Ks']
   Vls         = mytable['Vls']
   nest        = mytable['nest']
   dist        = mytable['dist']
   mDist       = mytable['mDist']
   No_Galaxies = mytable['No_Galaxies']
   sigmaP_dyn  = mytable['sigmaP_dyn']

   N = len(id)
   NoGroups = len(id[np.where(flag==2)])
   NoGals = len(id[np.where(flag!=2)])
   NoFreeGals = len(id[np.where(flag==0)])
   #print "Number of groups: ", NoGroups
   #print "Number of galaxies: ", NoGals
   #print "Number of free gals: ", NoFreeGals
  
   if NoGroups == 0:
     print "[Warning] No group found in data-base ..." 
     print "Check the input catalog and choose the right option ...\n" 
     return

   i = 0 
   if NoGroups!=0:
     while flag[i] != 2:
        i+=1

   GroupList = []

   
   while i<N:
     
     id_grp       = [id[i]]
     flag_grp     = [flag[i]]
     ra_grp       = [ra[i]]
     dec_grp      = [dec[i]]
     Ks_grp       = [Ks[i]]
     Vls_grp      = [Vls[i]]
     nest_grp     = [nest[i]]
     dist_grp     = [dist[i]]
     mDist_grp    = [mDist[i]]
     No_Galaxies_grp = [No_Galaxies[i]]
     sigmaP_dyn_grp  = [sigmaP_dyn[i]]
     
     i += 1
     while i<N and flag[i]==1: 
         
         id_grp.append(id[i])
         flag_grp.append(flag[i])
         ra_grp.append(ra[i])
         dec_grp.append(dec[i])
         Ks_grp.append(Ks[i])
         Vls_grp.append(Vls[i])
         nest_grp.append(nest[i])
         dist_grp.append(dist[i])
         mDist_grp.append(mDist[i])
         No_Galaxies_grp.append(No_Galaxies[i])
         sigmaP_dyn_grp.append(sigmaP_dyn[i])
         
         i+=1
     
     mygroup = {'pgc':id_grp, 'flag':flag_grp, 'ra':ra_grp, 'dec':dec_grp, 'Ks':Ks_grp, 'Vls':Vls_grp,  \
               'nest':nest_grp, 'dist':dist_grp, 'mDist':mDist_grp, 'No_Galaxies':No_Galaxies_grp, 'sigmaP_dyn':sigmaP_dyn_grp} 
     
     
     NewGroup = GroupNode(data=mygroup)
     
     if NewGroup.ra>=RA[0] and NewGroup.ra<RA[1] and NewGroup.dec>=DEC[0] and NewGroup.dec<DEC[1]: 
         GroupList.append(NewGroup)


   print "##########################"
   print "No. of Groups: ", len(GroupList)
   print "RA   = ", RA
   print "DEC  = ", DEC
   print "Data file loaded .... "
   print "##########################"
   print

   return GroupList
# **************************************
class GroupNode:
   
   def __init__(self, data=None):
     
     
     
      self.dist = 0
      ########self.mDist = 0
      

      self.ra         = 0
      self.dec        = 0
      self.Ks         = 0
      self.logK       = 0
      
      self.M_v2       = 0
      self.R_2t2      = 0
      self.r2t        = 0 
      self.r1t        = 0
      
      self.flag   = 0   #     
      
      if data == None:
	self.id        = 0
	self.Vls       = 0
	self.sigma     = 0 
	self.nest      = 0
	self.mDist     = 0
	self.Vls_list  = None
	self.ngal      = 0
      
        return
      
      #self.data = data   ### ????
      
      id              = data['pgc']
      flag            = data['flag']
      ra              = data['ra']
      dec             = data['dec']
      Ks              = data['Ks']
      nest            = data['nest']
      dist            = data['dist']
      mDist           = data['mDist']
      self.Vls_list   = data['Vls']
      ngal            = data['No_Galaxies']
      sigmaP_dyn      = data['sigmaP_dyn']

      
      #### ??? Radius = r1t, Mass=M_v2, R_2t2
      self.id         = id[0]
      self.ngal       = ngal[0]
      self.nest       = nest[0]
      self.flag       = flag[0]
      
      if self.ngal>1: 
        #self.Vls, self.sigma      = self.v_ave(self.Vls_list[1:])
        #self.mDist, self.mDistErr = self.dist_av(dcf2, ed)
        self.mDist = mDist[0]
        self.Vls = self.Vls_list[0]
        self.sigma = sigmaP_dyn[0]
        

        
        self.flag                 = 2
      else: 
	self.Vls      = self.Vls_list[0]
	self.sigma    = 0.
	self.mDist    = dist[0]
	self.flag     = 0
	
      self.initialize(ra, dec, Ks)
   

   # ************
   def initialize(self, ra, dec, Ks):
     
     Ltot = 0

     Mra=0; Mdec=0
     
     
     Dist_v = self.mDist

     if Dist_v==0: # or self.mDistErr > 0.10: 
       Dist_v =  self.Vls / H0
       if Dist_v<1. : Dist_v=1
     
     if self.ngal>1:
	for i in range(1, len(Ks)):
	  
	  L = 10**m_logK(Ks[i], self.Vls, distance=Dist_v)

	  Mra += L*ra[i]
	  Mdec += L*dec[i]       
	  Ltot += L

	  self.ra   = Mra/Ltot
	  self.dec  = Mdec/Ltot
	  
	  self.logK = log10(Ltot)	  
     else:
          L = 10**m_logK(Ks[0], self.Vls, distance=self.mDist)
          self.logK = log10(L)
 	  self.ra   = ra[0]
	  self.dec  = dec[0]
     


     Mk_sun = 3.28   # K-band
     M = Mk_sun - (self.logK / 0.4)

     
     self.dist =   Dist_v
     self.Ks = M + 5*log10(Dist_v) + 30 - 5
     
     self.M_v2 = Mass(10**self.logK)
     self.R_2t2 = 0.215*((self.M_v2/1.E12)**(1./3))  # Mpc
     
     self.r2t = sqrt(1.5)*self.R_2t2
     self.r1t = 3.5*self.r2t
   
   # ************
   def dist_av(self, dist):
     
     sum_dist = 0
     n = 0 
     mD=0

     for i in range(1, len(dist)):
       
       if dist[i] != 0:
	 sum_dist += dist[i]
	 n += 1 
    
     if n!=0:
       mD = sum_dist/n

     return mD
     
   
   
   
   def v_ave(self, v_list):
     
     m = 0 
     sum_v = 0
     sum_v2 = 0 
     
     for v in v_list:
       if v != 0:
	 sum_v += v
	 sum_v2 += v**2
	 m +=1
     
     if m != 0 : 
        mean_v = sum_v/m
        sigma_v = sqrt((sum_v2/m)-mean_v**2)
     else:
        mean_v = 0
        sigma_v = 0
     
     return mean_v, sigma_v
     
#################################################################
def m_logK(m, Vls, distance=0):
    
    Mk_sun = 3.28   # K-band
    
    if distance == 0: 
       distance = Vls / H0  # Mpc    
       if distance < 1:
          distance = 1 # Mpc
    
    M = m - 5*log10(distance) - 30 + 5
    logK = -0.4 * (M - Mk_sun)
    
    if Vls==0 and distance==0:
      logK = 0
    
    
    return logK
# **************************************
# If old_SuperList=None then only use GroupList
# If old_SuperList is specified, then attach GroupList to it

def SuperGroup(GroupList, old_SuperList=None):
  
  
  if old_SuperList != None:
    
    inGame = []
    Lum_mass = []
    for SuperGroup in old_SuperList:
      Lum_mass.append(SuperGroup[0].M_v2)
    indices = np.argsort(Lum_mass)
    for i in indices[::-1]: 
    #for i in indices:
      for group in old_SuperList[i][1]:
          inGame.append(group)
          
    for group in GroupList:
      if group.flag == 0 or group.flag == 2:
	inGame.append(group)
	
    #for group in GroupList:
      #group.flag = 0

  else:
  # -------------------------

    inGame = GroupList[:]
    
    Lum_mass = []
    for i in range(0, len(inGame)): Lum_mass.append(inGame[i].M_v2)
    Lum_mass = np.asarray(Lum_mass)
    indices = np.argsort(Lum_mass)
    
    tmp = []
    for i in indices[::-1]: tmp.append(inGame[i]) 
    #for i in indices: tmp.append(inGame[i]) 
    inGame = tmp
    
  # -------------------------
  
  
  SuperList = []
  N = len(inGame)
  #print N
  
  p = 0 
  while p < N-1:
    
    Players = []
    Players.append(inGame[p])
    
    q = (p+1)%N
    pivot = p
    head = Players[0]
    Bol = False
    while q!=pivot:  
      
      if touch(head, inGame[q], coeff=1.0, error='No'):
	Bol = True
	Players.append(inGame[q])
	head, Players =  Gjoint(Players)
	
	inGame.pop(q)
	N-=1
	if p>=q:  p-=1
  
	if len(Players) == 2:
	   inGame.pop(p)
	   N-=1	
	   if q>p: q = (q-1)%N    
	
	q = (q-1)%N
	pivot = q  
  
      q = (q+1)%N
    if Bol: 
      SuperList.append([head, Players])
    p+=1 
    
  return SuperList
#################################################################


def isin(id, list):
  
  for obj in list:
    
    if obj.id == id: return True
  
  return False



#################################################################
def forceMerge(SuperList, gID1, gID2):
  
  if gID1==None or gID2==None: return
  SGroup1 = None
  SGroup2 = None
  
  i=0
  while i < len(SuperList):
     objlist = SuperList[i][1]
     for objects in objlist:
       if extractPGC(objects.id, grp=True) == gID1:
          SGroup1 = objlist  
          i1=i
          i = 100000000 # break while
          break # break for
     i+=1
  
  if  SGroup1==None: 
    print 'Could not find a group for:', gID1
    return 
  
  i=0
  while i < len(SuperList):
     objlist = SuperList[i][1]
     for objects in objlist:
       if extractPGC(objects.id, grp=True) == gID2:
          SGroup2 = objlist  
          i2=i
          i = 100000000 # break while
          break # break for
     i+=1

  if  SGroup2==None: 
    print 'Could not find a group for:', gID2
    return 
  
  if i1==i2:
    print 'Both are in the same group. Nothing done ...'
    
  SuperList.pop(i1)
  if i1<i2:
    SuperList.pop(i2-1)
  else:
    SuperList.pop(i2)
  
  Players = []
  
  for group in SGroup1:
    Players.append(group)
  
  for group in SGroup2:
    Players.append(group)
    
  head, Players =  Gjoint(Players)	  
 
  SuperList.append([head, Players])
  return


#################################################################
def mergSGroup(SuperList, projected=0):
  
    NoGroups = len(SuperList)
 
    i = 0
    while i < NoGroups:
      j = i+1
      while j < NoGroups-1:
	
	if touch(SuperList[i][0], SuperList[j][0], coeff=1.0, projected=projected):
	  
	    
		  Players = []
		  
		  for group in SuperList[i][1]:
		    Players.append(group)
		  
		  for group in SuperList[j][1]:
		    Players.append(group)
		    
		  head, Players =  Gjoint(Players)	  
		  
		  SuperList[i] = [head, Players]
		  SuperList.pop(j)
		  NoGroups-=1
		  i=0
		  j=0
        j+=1
      i+=1
    
    return SuperList  
  
#################################################################

def Gjoint(Players):
   
   new_Players = []
   if len(Players)==0: return None

   if len(Players)==1: return Players[0]

   Lum_mass = []
   for i in range(0, len(Players)):
        Lum_mass.append(Players[i].M_v2)
   Lum_mass = np.asarray(Lum_mass)
   indices = np.argsort(Lum_mass)
   
   id = Players[indices[(len(indices)-1)]].id
   
   head = None
   for i in indices[::-1]:
       head = Joint(head, Players[i], ID = id+10000000000)
       new_Players.append(Players[i])
       
   return head, new_Players
 
#################################################################


def Joint(gr1, gr2, ID = 0):
  
   if gr1 == None and gr2 != None:
     return gr2
   if gr1 != None and gr2 == None:
     return gr1
   if gr1 == None and gr2 == None:
     return None  
   
   gr = GroupNode()
   gr.id = ID
   
   L1 = 10**gr1.logK
   L2 = 10**gr2.logK
   
   if L1 == 1: L1 = 0
   if L2 == 1: L2 = 0
   
   L_tot = L1 + L2
   
   gr.logK = log10(L_tot)
   
   
   d1 = gr1.dist
   d2 = gr2.dist
   
   ra1 = gr1.ra
   ra2 = gr2.ra
   dec1 = gr1.dec
   dec2 = gr2.dec

   
   ra1 = gr1.ra
   ra2 = gr2.ra
   dec1 = gr1.dec
   dec2 = gr2.dec   

   v1 = gr1.Vls
   v2 = gr2.Vls   
   gr.Vls  = (L1*v1 + L2*v2)/L_tot
   
   gr.ra, gr.dec, d     = barycenter(L1,  ra1,  dec1, d1, L2,  ra2,  dec2, d2)

     
   gr.dist = d
   Mk_sun = 3.28   # K-band
   M = Mk_sun - (gr.logK / 0.4)
   gr.Ks = M + 5*log10(gr.dist) + 30 - 5
   
   #gr.M_v2 = Mass(10**gr.logK)
   gr.M_v2 = gr1.M_v2 + gr2.M_v2   # we simply add the mas of halos at this stage 
   
   
   gr.R_2t2 = 0.215*((gr.M_v2/1.E12)**(1./3))  # Mpc
   gr.r2t = sqrt(1.5)*gr.R_2t2
   gr.r1t = 3.5*gr.r2t
   
   gr.mDist =  d
   

   gr.flag = 5
   
   
   if  gr1.ngal>1:
     gr1.flag = 4
   else: 
     gr1.flag = 3
     
   if  gr2.ngal>1:
     gr2.flag = 4
   else: 
     gr2.flag = 3     
   
   gr.ngal = gr1.ngal + gr2.ngal
   gr.Vls_list = np.asarray([gr.Vls])
   gr.Vls_list = np.concatenate((gr.Vls_list, gr1.Vls_list[1:], gr2.Vls_list[1:]))
   
   
   
   tmp, gr.sigma = gr.v_ave(gr.Vls_list[1:])
   
   
   return gr
  
#################################################################
#################################################################   
def SGroupwrite(outfile, SuperList, GroupList):
  
   NoSGroups = len(SuperList)
   
   myTable = Table()
   
   empty = []
   myTable.add_column(Column(data=empty,name='ID', dtype=np.dtype(int)))
   myTable.add_column(Column(data=empty,name='flag', dtype=np.dtype(int)))
   myTable.add_column(Column(data=empty,name='ra', format='%0.4f'))
   myTable.add_column(Column(data=empty,name='dec', format='%0.4f', length=10))


   myTable.add_column(Column(data=empty,name='Ks', format='%0.2f'))
   myTable.add_column(Column(data=empty,name='logK', format='%0.4f'))
   myTable.add_column(Column(data=empty,name='Vls', format='%0.0f'))
   myTable.add_column(Column(data=empty,name='dist', format='%0.2f'))
   myTable.add_column(Column(data=empty,name='mDist', format='%0.2f'))

   myTable.add_column(Column(data=empty,name='sigmaP_dyn', format='%0.1f'))
   myTable.add_column(Column(data=empty,name='sigmaP_lum', format='%0.1f'))
  
   myTable.add_column(Column(data=empty,name='Mv_lum', format='%1.2e'))
   myTable.add_column(Column(data=empty,name='R2t_lum', format='%0.3f'))
   myTable.add_column(Column(data=empty,name='r1t_lum', format='%0.3f'))  
   myTable.add_column(Column(data=empty,name='tX_lum', format='%1.2e'))  
   
   myTable.add_column(Column(data=empty,name='No_Galaxies',dtype=np.dtype(int)))
   myTable.add_column(Column(data=empty,name='nest', dtype=np.dtype(int)))
   
   for i in range(0, NoSGroups):  # for all groups
        meanDist = 0.
        sumDist = 0
        n = 0
        for object in SuperList[i][1]:
	  if object.mDist!=0:
	    sumDist += object.ngal*object.mDist
	    n += object.ngal
	  
	if sumDist !=0 and n != 0:
	  meanDist = sumDist/n
	    
	SuperList[i][0].mDist = meanDist
	SuperList[i][0].nest = SuperList[i][1][0].nest
     
        table_row(myTable, SuperList[i][0])
        for Group in SuperList[i][1]:
	  table_row(myTable, Group)
   for Group in GroupList:
     if Group.flag<=2: 
       table_row(myTable, Group)
   
   
   pgc = 999999999; 
   ra = 999.9999; dec=-99.99;
   gl = ra; gb = dec
   ra = ra; dec=dec
   Ty = -100000.00; B_mag=Ty
   Ks= 99.99
   logK = 99.9999
   Vls = 9999
   dcf2 = 99.99
   ed = 9.99
   Mv_dyn = 9.99E99; Mv_lum = Mv_dyn
   tX_dyn = Mv_lum; tX_lum=Mv_lum
   nest = 9999999   
   
   flag = 0
   mDist = 0
   dist = 0
   sigmaP_dyn = 0
   sigmaP_lum = 0 
   R2t_lum = 0
   r1t_lum = 0
   subGalaxies = 0
   
   myTable.add_row([pgc,flag,ra,dec,Ks,logK,Vls, dist, \
	       mDist, sigmaP_dyn, sigmaP_lum, \
	          Mv_lum, R2t_lum, r1t_lum, tX_lum, subGalaxies, nest])
   
   
   myTable.write(outfile, format='ascii.fixed_width',delimiter='|', bookend=False)
   
   ### removing the last line, (it sits o adjust the column wodths)
   command =  ["csh", "remove_lastline.csh", outfile]
   subprocess.call(command)  
   #print myTable
#################################################################
# adding a row to a table based on the input Group parameters
def table_row(myTable, Group):
  
        
        pgc    = Group.id  
        flag   = Group.flag
        ra = Group.ra
        dec = Group.dec

        
        Ks = Group.Ks
        logK = Group.logK  
        
        Vls = Group.Vls  
        dist = Group.dist
        
        mDist = Group.mDist
        
        subGalaxies = Group.ngal  
        sigmaP_dyn = Group.sigma  
        nest = Group.nest  
        
        Mv_lum = Group.M_v2
        R2t_lum = Group.R_2t2
        r1t_lum = Group.r1t
        
        sigmaP_lum = (Mv_lum / (2.0E6))**(1./3)
        tX_lum = Mpc_km*R2t_lum*1.05*sqrt(1.5)/sigmaP_lum/sqrt(2.5)

        myTable.add_row([pgc,flag,ra,dec,Ks,logK,Vls, dist, \
	       mDist, sigmaP_dyn, sigmaP_lum, \
	          Mv_lum, R2t_lum, r1t_lum, tX_lum, subGalaxies, nest])

#################################################################
def reset(GroupList):
  
  for grp in GroupList:
    if grp.flag == 4: grp.flag = 2
    if grp.flag == 3: grp.flag = 0
  
  return
#################################################################

def head_list(SuperList):
  
  list = []
  for SuperGroup in SuperList:
    list.append([SuperGroup[0],[]])
  
  
  return list
#################################################################
def force_add(SuperList, GroupList, SGID=None , obj_id=None):

  if SGID==None or obj_id==None: return
 
  i = 0
  sgid = None
  if SGID != None:
    while i < len(SuperList):
     objlist = SuperList[i][1]
     for objects in objlist:
       if extractPGC(objects.id, grp=True) == SGID:
          sgid = SuperList[i][0].id
          i = 1000000000 # break while
          #print 'found it'
          break # break for
     i+=1   
  
  if  sgid == None: 
    #print 'Error ...'
    return 
  
  obj = None
  for grp in GroupList:
    if extractPGC(grp.id, grp=True) == obj_id:
      obj = grp

  
  if  obj == None: 
    #print 'Error ...'
    return 
  
  #print 'flag: ', obj.flag
  remove_obj(SuperList, GroupList, obj_id, grp=True)
  
  i = 0
  bol = True
  while i < len(SuperList) and bol:
      if SuperList[i][0].id == sgid:
	#print sgid
	obj_list = SuperList[i][1]
	obj_list.append(obj)
	SuperList.pop(i)
	header, obj_list = Gjoint(obj_list)
	SuperList.append([header, obj_list])
	bol = False
	#print 'done'
      i+=1



#################################################################
def remove_obj(SuperList, GroupList, id, grp=False, gal=True):  # removing obj from SG-list

  i = 0
  bol = True
  while i < len(SuperList) and bol:
    j = 0 
    while j < len(SuperList[i][1]) and bol:
      
      if not grp and gal:  # it's a galaxy ID
	myID = SuperList[i][1][j].id
      if grp:  # it is a group ID
	myID = extractPGC(SuperList[i][1][j].id, grp=True)
      
      if myID == id:
	#print 'removed ....'
	if SuperList[i][1][j].flag == 3: SuperList[i][1][j].flag = 0
	if SuperList[i][1][j].flag == 4: SuperList[i][1][j].flag = 2
	SuperList[i][1].pop(j)
	obj_list = SuperList[i][1]
	SuperList.pop(i)
	if len(obj_list) > 1:
	  header, obj_list = Gjoint(obj_list)
	  SuperList.append([header, obj_list])
	else:
	  if obj_list[0].flag == 3: obj_list[0].flag = 0
	  if obj_list[0].flag == 4: obj_list[0].flag = 2
	  
	bol = False
      j+=1
    i+=1
  
  
#################################################################




#################################################################
def populate(heads, GroupList, except_id=None):
  
  #print except_id
  SuperList = []
  i = 0
  if except_id != None:
    while i < len(heads):
     objlist = heads[i][1]
     for objects in objlist:
       if extractPGC(objects.id, grp=True) == except_id:
          #print "TEST !!!: ", except_id, heads[i][0].id
          i = 1000000000 # break while
          break # break for
     i+=1
    
  
  for grp in GroupList:
    
    if grp.flag == 2 or grp.flag == 0:
    
	indx = []
	dist = []
	for i in range(len(heads)):
	  
	  d = distance(heads[i][0], grp, coeff=1.0)
	  if d != None and (except_id==None or extractPGC(heads[i][0].id, supergrp=True) != except_id):
	    indx.append(i)
	    dist.append(d)
	
	indx = np.asarray(indx)
	dist = np.asarray(dist)

	if len(dist) > 0:
	    indices = np.argsort(dist)
	    #if len(indices) > 1:
	    p = indx[indices][0]
	      #p = p[0]
	    #else:
	      #p = indx[indices][0]
	    
	    heads[p][1].append(grp)


  for i in range(len(heads)):
    
    if len(heads[i][1]) > 1:
      
      header, grouplist = Gjoint(heads[i][1])
      SuperList.append([header, grouplist])
  
  return  SuperList

#################################################################
# removing objects from a supergroup
def trimGroup_vel(SuperList, id, Vmin=-100000, Vmax=100000):

  SGroup = None
  i=0
  while i < len(SuperList):
     objlist = SuperList[i][1]
     for objects in objlist:
       if extractPGC(objects.id, grp=True) == id:
          SGroup = objlist  
          SuperList.pop(i)
          i = 100000000 # break while

          break # break for
     i+=1
  
  if  SGroup==None: return 
  else: 
    
    newGroup = []
    for p in range(0, len(SGroup)):
        obj = SGroup[p]

        if obj.Vls < Vmax and obj.Vls > Vmin:
	  newGroup.append(obj)
	else:

	   if obj.flag == 4: obj.flag = 2
	   if obj.flag == 3: obj.flag = 0
    
    if len(newGroup) > 0:
      
      header, newGroup = Gjoint(newGroup)
      SuperList.append([header, newGroup])
      
#################################################################
def destryoySgroup(SuperList, SGID=2557):
  
  if SGID==None: return

  i = 0
  sgid = None
  if SGID != None:
    while i < len(SuperList):
     objlist = SuperList[i][1]
     for objects in objlist:
       if extractPGC(objects.id, grp=True) == SGID:
          sgid = SuperList[i][0].id
          obj_list = SuperList[i][1]
          SuperList.pop(i)
          i = 1000000000 # break while
          #print 'found it'
          break # break for
     i+=1   

  if  sgid == None: 
    #print 'Error ...'
    return 
  
  if len(obj_list) > 1:
    for obj in obj_list:
      print "tst: ", obj.id
      if obj.flag == 3: obj.flag = 0
      if obj.flag == 4: obj.flag = 2
  
  return
#################################################################
def formSgroup(SuperList, GroupList, idlist=None):
  
  if idlist==None:
    print "Warning: No ID list has been provided :( "
    return
  
  for id in idlist:
    destryoySgroup(SuperList, SGID=id)
  
  obj_list = []
  for id in idlist:
    obj = None
    for grp in GroupList:
      if grp.flag==2:
	myID = extractPGC(grp.id, grp=True)
      else:
	myID = grp.id
      
      if myID==id:
	remove_obj(SuperList, GroupList, id, grp=True)  
	obj = grp
	break
     
    if obj!=None:
      obj_list.append(obj)
      
  
  header, obj_list = Gjoint(obj_list)
  SuperList.append([header, obj_list])

#################################################################
def Make_SuperGroup(GroupList):


   print 
   print "[NEW] working on SuperGrouping ... "+str(len(GroupList))+' Groups'
    
   SuperList = SuperGroup(GroupList)
   SuperList = SuperGroup(GroupList, old_SuperList=SuperList)
   
   
   print "Populating ..."
   for qp in range(3):
    reset(GroupList)
    heads = head_list(SuperList)
    SuperList = populate(heads, GroupList)
   
   print "Mergeing ...."
   for qp in range(10):
      mergSGroup(SuperList)       

   ### Populate all Sgroups except 13418 (it's big, to avoid capturing everything)
   #SuperList = populate(SuperList, GroupList, except_id=13418)   
   ## except_id is the ID of one of the SGroups members, it must be the PGC id of the main 
   ## group galaxy
   
   return SuperList
#################################################################
def Super_Long_Patch(ra_del, ra_overlap, SGroups, Groups, DEC=[0,90]):
  
  print len(SGroups), len(Groups)
          
  for ra_i in np.arange(ra_del,360,ra_del):
    print '.......................................'
    print "Revisiting Patch .... RA,DEC = ", [ra_i-ra_overlap,ra_i+ra_overlap], DEC
    patch_GroupList = Remove_Patch(SGroups, Groups, RA=[ra_i-ra_overlap,ra_i+ra_overlap], DEC=DEC)
    SGroups += Make_SuperGroup(patch_GroupList)
    Groups   += patch_GroupList


  print '.......................................'
  print "Revisiting Patch .... RA,DEC = ", [360-ra_overlap,0+ra_overlap], DEC
  patch_GroupList = Remove_Patch(SGroups, Groups, RA=[360-ra_overlap,0+ra_overlap], DEC=DEC, RAA=True)
  SGroups += Make_SuperGroup(patch_GroupList)
  Groups   += patch_GroupList    

  

#################################################################
#################################################################
def Remove_Patch(SGroups, Groups, RA=[0,10], DEC=[0,10], RAA=False):
    
    patch_GroupList = []
    
    N_grps = len(Groups)
    N_Sgrps = len(SGroups)
    
    i = 0
    while i<N_grps:
        grp = Groups[i]
        
        Bol = False
        if RAA and (grp.ra>=RA[0] or grp.ra<=RA[1]) and grp.dec>=DEC[0] and grp.dec<DEC[1]:
            Bol = True
        
        if Bol or (grp.ra>=RA[0] and grp.ra<RA[1] and grp.dec>=DEC[0] and grp.dec<DEC[1]):
            if grp.flag == 2 or grp.flag == 0:
                Groups.pop(i)
                patch_GroupList.append(grp)
                i-=1
                N_grps-=1
        i+=1    
        
    remove_ids = [] 
    i = 0   
    while i<N_Sgrps:
        
        
        Sgroup = SGroups[i]
        Bol = False
        if RAA and (Sgroup[0].ra>=RA[0] or Sgroup[0].ra<=RA[1]) and Sgroup[0].dec>=DEC[0] and Sgroup[0].dec<DEC[1]:
            Bol = True


        if Bol or (Sgroup[0].ra>=RA[0] and Sgroup[0].ra<RA[1] and Sgroup[0].dec>=DEC[0] and Sgroup[0].dec<DEC[1]):
            SGroups.pop(i)
            i-=1
            N_Sgrps-=1
            
            for grp in Sgroup[1]:
                remove_ids.append(grp.id)
        i+=1
                
                
    N_grps = len(Groups)
    i = 0
    while i<N_grps:
        grp = Groups[i]
        if grp.id in remove_ids:
            patch_GroupList.append(Groups.pop(i))
            i-=1
            N_grps-=1
        i+=1             
    
    for grp in patch_GroupList:
      if grp.flag == 4: grp.flag = 2
      if grp.flag == 3: grp.flag = 0

    return patch_GroupList

#################################################################
#################################################################
def Equatorial_Patch(ra_del, dec, SGroups, Groups):
    
  for ra_i in np.arange(20,360-ra_del,ra_del):
      ra_f = ra_i+ra_del
      print '.......................................'
      print "Equatorial Ring RA,DEC= ", [ra_i, ra_f], [dec-20,dec+20]
      patch_GroupList = Remove_Patch(SGroups, Groups, RA=[ra_i, ra_f], DEC=[dec-20,dec+20])
      SGroups += Make_SuperGroup(patch_GroupList)
      Groups   += patch_GroupList
  
  print '.......................................'
  print "Equatorial Ring RA,DEC= ", [ra_f, 20], [dec-20,dec+20]
  patch_GroupList = Remove_Patch(SGroups, Groups, RA=[ra_f, 20], DEC=[dec-20,dec+20], RAA=True)
  SGroups += Make_SuperGroup(patch_GroupList)
  Groups   += patch_GroupList    
  
  
  
        
#################################################################

if __name__ == '__main__':

   
   table = np.genfromtxt('M_L_curve_v2.csv' , delimiter=',', filling_values=0, names=True, dtype=None)
   Lumin_in   = 10**table['log10_L']
   Mass_out   = 10**table['log10_M']
  
   LtoM_func = interpolate.interp1d(Lumin_in, Mass_out)

 

   
   filename = "halo.0.100.v01.group" 
   #filename = 'simul.0.100.2D.v01.group.tmp4'
   #filename = "simul.MagK16.d5.d100.2D.v01.group"
   
   print filename
   
   ra_del = 45.
   ra_overlap = 22.5
  
   Dec_lst = [-90,-45,0,45,90]
   #Dec_lst = [-90,0,90]
  
  
  #############################################################
     
   GroupLists = []
   
   for ra_i in np.arange(0,360,ra_del):
     ra_f = ra_i+ra_del
     for i in range(len(Dec_lst)-1):
        GroupLists.append(readgrouplist(filename, RA=[ra_i, ra_f], DEC=[Dec_lst[i],Dec_lst[i+1]]))


  #############################################################
  
   Groups = []
   for GroupList in GroupLists:
      Groups += GroupList 
   print  ' Total # of Groups: ', len(Groups)
   
   
   SuperG_lists = []
   N = len(GroupLists)
   for i in range(N):
      GroupList = GroupLists[i]
      print
      print "Patch # "+ str(i+1)+'/'+str(N)
      SuperG_list = Make_SuperGroup(GroupList)
      SuperG_lists.append(SuperG_list)   

########################################################################
        
   SGroups = []
   for SGlist in SuperG_lists:
      SGroups += SGlist
  
   Groups = []
   for GroupList in GroupLists:
      Groups += GroupList
  
########################################################################

   for i in range(len(Dec_lst)-1):
      Super_Long_Patch(ra_del, ra_overlap, SGroups, Groups, DEC=[Dec_lst[i],Dec_lst[i+1]])
   


   for i in range(1,len(Dec_lst)-1):
      Equatorial_Patch(ra_del, Dec_lst[i], SGroups, Groups)


   print '.......................................'
   print "North Cap "
   patch_GroupList = Remove_Patch(SGroups, Groups, RA=[0,360], DEC=[60,90])
   SGroups += Make_SuperGroup(patch_GroupList)
   Groups   += patch_GroupList


   print '.......................................'
   print "South Cap "
   patch_GroupList = Remove_Patch(SGroups, Groups, RA=[0,360], DEC=[-90,-60])
   SGroups += Make_SuperGroup(patch_GroupList)
   Groups   += patch_GroupList



#####################################################################

   print "Final Mergeing ...."
   for qp in range(10):
      mergSGroup(SGroups)  
#####################################################################

   outFile = 'halo.0.100.v01.supergroup'   # 'tmp.supergroup2'
   #outFile = 'simul.MagK16.d5.d100.2D.v01.supergroup' 
   print "\n total super groups #:", len(SGroups)
   SGroupwrite(outFile, SGroups, Groups)
   print "Created: "+outFile
   print    

