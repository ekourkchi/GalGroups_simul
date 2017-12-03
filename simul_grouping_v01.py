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
import astropy.stats.funcs as st

# **************************************
# Global Variables
# Physical Constants
# **************************************
H0 = 75.           # hubble constant
sglV = 102.8806    # M87 center - super galactic longitude
sgbV = -2.3479     # M87 center - super galactic latitude
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
  
  if L_k==0:
    return 0
  
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
   
def distance(l1,b1,d1,l2,b2,d2):
    
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
   
   return R12
  

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
  
   if d1==0 or d2==0:
       
       dd1 = 1.
       dd2 = 1.
   else: 
       dd1 = d1
       dd2 = d2
 
   
   cl1 = cos(l1*pi/180.)
   sl1 = sin(l1*pi/180.)
   cb1 = cos(b1*pi/180.)
   sb1 = sin(b1*pi/180.)
   
   x1 = dd1 * cl1 * cb1
   y1 = dd1 * sl1 * cb1
   z1 = dd1 * sb1
   
   cl2 = cos(l2*pi/180.)
   sl2 = sin(l2*pi/180.)
   cb2 = cos(b2*pi/180.)
   sb2 = sin(b2*pi/180.)
   
   x2 = dd2 * cl2 * cb2
   y2 = dd2 * sl2 * cb2
   z2 = dd2 * sb2   

   if L1==0 and L2==0:
     x = 0.5*(x1+x2)
     y = 0.5*(y1+y2)
     z = 0.5*(z1+z2)
   else: 
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

def galJoint(galNode1, galNode2, ID):
   
   if galNode1==None and galNode2!=None:
     return galNode2
   if galNode2==None and galNode1!=None:
     return galNode1
   if galNode2==None and galNode1==None:
     return None
   

   
   L1 = 10**galNode1.logK
   L2 = 10**galNode2.logK
   
   if L1 == 1: L1 = 0
   if L2 == 1: L2 = 0
   
   L_tot = L1 + L2
   
   try:
     logK_tot = log10(L_tot)
   except:
     logK_tot = 0
     
     
   
   
   d1 = galNode1.dist
   d2 = galNode2.dist
     
   ra1 = galNode1.ra
   ra2 = galNode2.ra
   dec1 = galNode1.dec
   dec2 = galNode2.dec   
      
   
   
   
   v1 = galNode1.Vls
   v2 = galNode2.Vls
   
   if d1==0 and d2==0:
     d1 = 0.
     d2 = 0
     d_final=0.
   elif d1!=0 and d2==0:
     d2 = galNode2.Vls/H0
     d_final=0
   elif d2!=0 and d1==0:
     d1 = galNode1.Vls/H0
     d_final=0
   else:
     d_final=1
     
     

   ra, dec, d = barycenter(L1, ra1, dec1, d1, L2, ra2, dec2, d2)
   
   if d_final==0:
      d = 0
      dm = 0
   else:
      dm = 5*log(d)+25


   n1 = galNode1.subGalaxies
   n2 = galNode2.subGalaxies   
   Vls = (n1*galNode1.v_av + n2*galNode2.v_av) / (n1+n2)
   if galNode1.v_av == 0 and galNode2.v_av != 0 :
     Vls = galNode2.v_av
     if L1 == 0:
       galNode1.logK = m_logK(galNode1.MagK)
       
   if galNode1.v_av != 0 and galNode2.v_av == 0 :
     Vls = galNode1.v_av   
     if L2 == 0:
       galNode2.logK = m_logK(galNode2.MagK)

   newNode = GalxyNode(ID, ra, dec, Vls, d, 0, 0, -1)
   newNode.logK = logK_tot
   newNode.ra = ra
   newNode.dec = dec

   
   newNode.subGalaxies = n1 + n2
   newNode.level = max([galNode1.level, galNode2.level]) + 1
   
   galNode1.topLevel = newNode.level
   galNode2.topLevel = newNode.level
   
   newNode.sumDist = galNode1.sumDist + galNode2.sumDist
   
   
   meanDist = newNode.sumDist/newNode.subGalaxies



    
   
   
   if galNode1.mDist != 0 and galNode2.mDist != 0:
     if L1 > L2: 
       newNode.mDist = galNode1.mDist
     else:
       newNode.mDist = galNode2.mDist

       
   elif galNode1.mDist != 0:
     newNode.mDist = galNode1.mDist    
   elif galNode2.mDist != 0:
     newNode.mDist = galNode2.mDist
   else:
     newNode.mDist = meanDist

   newNode.dist = newNode.mDist 
   
   # Brighter galaxy is the left child
   if L1 >= L2:
      newNode.left = galNode1
      newNode.right = galNode2
      newNode.nest = galNode1.nest
   else:
      newNode.left = galNode2
      newNode.right = galNode1
      newNode.nest = galNode2.nest
   

   newNode.v_av = Vls
   newNode.v2_av = (n1*galNode1.v2_av + n2*galNode2.v2_av) / (n1+n2)
   if galNode1.v2_av == 0 and galNode2.v2_av != 0 :
     newNode.v2_av = galNode2.v2_av
   if galNode1.v2_av != 0 and galNode2.v2_av == 0 :
     newNode.v2_av = galNode1.v2_av   
   
   
   
   if (newNode.v2_av - newNode.v_av**2) > 0:
      newNode.sigma =  sqrt(newNode.v2_av - newNode.v_av**2) 
   
   
   newNode.R_theta = Theta_max(newNode)

   if newNode.sigma == 0:
     sig = 1
   else:
     sig = newNode.sigma
   

   mass = Mass(L_tot)
   newNode.M_v2 = mass
   if np.isinf(newNode.M_v2): 
	    print '[ERROR] Infinity M_v2 ... Exit !'
	    sys.exit()
   newNode.R_2t2 = 0.215*((mass/1.E12)**(1./3))  # Mpc
   newNode.r2t = sqrt(1.5)*newNode.R_2t2
   newNode.r1t = 3.5*newNode.r2t

   return newNode
# **************************************
# The calls definition of a galaxy
# each node contains all esseential property of a galaxy
# when gaalxies get connected along a tree, the new nodes
# are defnided as new entities, which are new galaxies in
# the context of this code
class GalxyNode:
  
  # Class constructor
  def __init__(self, galid, ra, dec, Vls, dist, MagB, MagK, grp_id):
    
    
    self.id   = galid
    self.Vls  = Vls
    
    if dist==0:
        self.Ks   = MagK + 5*log10(0.00001) + 30 - 5
    else:
        self.Ks   = MagK + 5*log10(dist) + 30 - 5
    self.ra   = ra
    self.dec  = dec
    self.dist = dist
    self.MagB   = MagB
    self.MagK   = MagK
    self.grp_id = grp_id

   
    self.subGalaxies = 1
    self.level = 0
    self.Rgroup = 0.  
    self.R_theta = 0.
    self.nest = galid
    self.sigma = 0. 
    self.v_av = Vls
    self.v2_av = Vls**2
    
    # 0: if the galaxy is NOT in a group
    # 1: if the galaxy falls in a group
    self.inGroup = 0.
    
    self.left = None
    self.right = None
    self.flag = 0
    
    self.mDist = dist
    
    self.logK = m_logK(self.MagK)
 
    self.M_v2 = Mass(10**self.logK)
    self.R_2t2 = 0.215*((self.M_v2/1.E12)**(1./3))  # Mpc
    
    self.r2t = sqrt(1.5)*self.R_2t2
    self.r1t = 3.5*self.r2t

    if dist!=0:
        self.sumDist = 1. * dist
    else:
        self.sumDist = 0.
        
  def setMeanDist(self, Dist, GRP_vel = 0, dist=0):
    
      self.mDist = Dist
      if GRP_vel == 0 : 
	vel = self.Vls
      else: vel = GRP_vel
      self.logK = m_logK(self.MagK)
      self.M_v2 = Mass(10**self.logK)
      self.R_2t2 = 0.215*((self.M_v2/1.E12)**(1./3))  # Mpc
      self.r2t = sqrt(1.5)*self.R_2t2
      self.r1t = 3.5*self.r2t
      
    

# **************************************
# Given the apparent magnitude of a galaxy, m
# it returns the absolute luminosity 
# Vls is the radial velocity of the galaxy
# So, the distance to a galaxy is estimated 
# based on its velocity
def m_logK(MagK):
    
    Mk_sun = 3.28   # K-band

    logK = -0.4 * (MagK - Mk_sun)

    return logK
# **************************************
# It returns all leaf-nodes at the bottom of a gaalxy tree.
# root: is the root of the tree
# It works recursively, so for large trees, it might 
# exhaust the recursion stack, and the codes crashe
# so some exception handling is required in this case
def NodeLeaf(root):
  list = [root]
  NodeLeafCore(root, list)
  return list


#def NodeLeafCore(root, list):  
  #if root.left == None:
    #list.append(root)
  #else:
    #NodeLeafCore(root.left, list)
    #NodeLeafCore(root.right, list)



def NodeLeafCore(root, list): 
    
    stack = [root]
    while len(stack)>0:
        parent = stack.pop()
        if parent.left == None:
            list.append(parent)
        else:
            stack.append(parent.left)
            stack.append(parent.right)
    
################################################################
def readgalList(table, VLS=[0,4000], d=[0,100], RA=[0,360], DEC=[-90,90], MagK_limit=0):

  grp_id = table['grp_id']
  Vls    = table['v_rad']
  dist   = table['d_kpc']/1000.
  ra     =  table['ra']
  dec    = table['dec']
  MagB   = table['MagB']
  MagK   = table['MagK']

   
  N_galaxies = len(grp_id)

  print "Reading the data file .... "
  galList = []
  
  for i in range(N_galaxies):
    
    if dist[i]==0:
        ra[i]  = 0.00001
        dec[i] = 0.00001
        Vls[i] = 0.00001
               
    galid = 1000000+i
    if Vls[i]>=VLS[0] and Vls[i]<VLS[1] and dist[i]>=d[0] and dist[i]<d[1] and MagK[i]<=MagK_limit:         
       if ra[i]>=RA[0] and ra[i]<RA[1] and dec[i]>=DEC[0] and dec[i]<DEC[1]:
           
	   node = GalxyNode(galid, ra[i], dec[i], Vls[i], dist[i], MagB[i], MagK[i], grp_id[i])
	   galList.append(node)
	   if Vls[i]<0: print "Warning, Negative velocity, galid: ", galid
           
  
  print "##########################"
  print "No. of galaxies: ", len(galList)
  print "VLS  = ", VLS
  print 'Dist = ', d
  print "RA   = ", RA
  print "DEC  = ", DEC
  print "MagK_limit : ", MagK_limit
  print "Data file loaded .... "
  print "##########################"
  print
  return galList
################################################################  
# b is galactic latitude [deg]
# 2*theta + sin(2*theta) = Pi * sin(b)

def theta(b):
  
  if b == 90.: return pi/2.
  if b == -90.: return -pi/2.
  
  b = b*pi/180.
  theta0 = b
  theta1 = theta0-(2.*theta0 + sin(2.*theta0) - pi * sin(b))/(2.+2.*cos(2.*theta0))
  
  while abs(theta1-theta0) > 0.01:
    theta0 = theta1
    theta1 = theta0-(2.*theta0 + sin(2.*theta0) - pi * sin(b))/(2.+2.*cos(2.*theta0))
  
  return theta1
  
################################################################ 

  
################################################################# 
def groupWrite(outfile, G_list, galList):
  

  NoGroups = len(G_list)

  myTable = Table()
    
    
  empty = []
  myTable.add_column(Column(data=empty,name='ID', dtype=np.dtype(int)))
  myTable.add_column(Column(data=empty,name='grpID', dtype=np.dtype(int)))
  myTable.add_column(Column(data=empty,name='flag', dtype=np.dtype(int)))
  myTable.add_column(Column(data=empty,name='ra', format='%0.4f'))
  myTable.add_column(Column(data=empty,name='dec', format='%0.4f', length=10))

  myTable.add_column(Column(data=empty,name='MagK', format='%0.2f'))
  myTable.add_column(Column(data=empty,name='Ks', format='%0.2f'))
  myTable.add_column(Column(data=empty,name='logK', format='%0.4f'))
  myTable.add_column(Column(data=empty,name='Vls', format='%0.0f'))

  myTable.add_column(Column(data=empty,name='dist', format='%0.2f'))


  myTable.add_column(Column(data=empty,name='mDist', format='%0.2f'))

  myTable.add_column(Column(data=empty,name='R_theta', format='%0.5f'))
  myTable.add_column(Column(data=empty,name='sigmaP_dyn', format='%0.1f'))
  myTable.add_column(Column(data=empty,name='sigmaP_lum', format='%0.1f'))
  
  myTable.add_column(Column(data=empty,name='Mv_dyn', format='%1.2e'))
  myTable.add_column(Column(data=empty,name='Mv_lum', format='%1.2e'))
  myTable.add_column(Column(data=empty,name='Rg_angular', format='%0.3f'))
  myTable.add_column(Column(data=empty,name='Rg_dyn', format='%0.3f'))
  myTable.add_column(Column(data=empty,name='R2t_dyn', format='%0.3f'))
  myTable.add_column(Column(data=empty,name='R2t_lum', format='%0.3f'))
  myTable.add_column(Column(data=empty,name='r2t_lum', format='%0.3f'))
  myTable.add_column(Column(data=empty,name='tX_dyn', format='%1.2e'))
  myTable.add_column(Column(data=empty,name='tX_lum', format='%1.2e'))  
  myTable.add_column(Column(data=empty,name='No_Galaxies',dtype=np.dtype(int)))
  myTable.add_column(Column(data=empty,name='nest', dtype=np.dtype(int)))


  
  
  #print "# of all groups: ", NoGroups
  
  for i in range(0, NoGroups):  # for all groups
    
   
    if G_list[i][0].Vls <= 3500:
       dist = G_list[i][0].mDist
       if dist == 0 : dist = G_list[i][0].Vls / H0
       if dist<1: dist = 1
       Rg =  dist * Rg_radian(G_list[i][1:]) # old version
       
       Rg_angular = dist * G_list[i][0].R_theta


       ra_lst = []
       dec_lst = []
       v_galaxies = []
       for gal in G_list[i][1:]:
          ra_lst.append(gal.ra)
          dec_lst.append(gal.dec)
          v_galaxies.append(gal.Vls)
       
       
       mean = st.biweight_location(v_galaxies)
       Rg_bi = biweight_Rg(ra_lst,dec_lst, mean, dist)
       Rg = Rg_bi    # bi-weight version
       #print 'Whisle ...:', pgc, dist , Rg_bi, Rg
       G_list[i][0].sigma  = st.biweight_midvariance(v_galaxies)

       for j in range(0, len(G_list[i])):  # for all galaxies
        
        
        galaxy = G_list[i][j]
        flag = galaxy.flag
        ID = galaxy.id  
        grp_id = galaxy.grp_id
        ra = galaxy.ra  
        dec = galaxy.dec
 
        Vls = galaxy.Vls  
        
        
        
        logK = galaxy.logK  
        Ks = galaxy.Ks
        dist = galaxy.dist

        mDist = galaxy.mDist
        
        if flag == 2:
              distance = mDist
        else: 
              distance = dist
              
        if distance==0:
            distance = Vls / H0

        MagK = galaxy.Ks - 5*log10(distance) - 30 + 5
        
        
        subGalaxies = G_list[i][0].subGalaxies  
        R_theta = G_list[i][0].R_theta  
        sigmaP_dyn = G_list[i][0].sigma  
        nest = G_list[i][0].nest  
        Mv_lum = galaxy.M_v2
        R2t_lum = galaxy.R_2t2
        r2t_lum = galaxy.r2t
        
        
        sigmaP_lum = (Mv_lum / (2.0E6))**(1./3)
        tX_lum = Mpc_km*R2t_lum*1.05*sqrt(1.5)/sigmaP_lum/sqrt(2.5)
        if j != 0: 
	  Rg = 0  # dynamical virial radius
	  R_theta = 0
	  Rg_angular = 0
        Mv_dyn = (1.E9 * (2.5*pi/G/2.) * (Mpc_km*Rg) * sigmaP_dyn**2)/M_sun  # solar mass
        
        if sigmaP_dyn == 0: 
	  tX_dyn = 0 
	else:
          tX_dyn = Mpc_km*0.5*pi*Rg/sigmaP_dyn/sqrt(2.5)
        R2t_dyn = (Rg*pi*0.5)/1.05/sqrt(1.5)
       
        myTable.add_row([ID,grp_id,flag,ra,dec,MagK,Ks,logK,Vls, dist, \
	       mDist, R_theta, sigmaP_dyn, sigmaP_lum, \
	       Mv_dyn, Mv_lum, Rg_angular, Rg, R2t_dyn, R2t_lum, r2t_lum, tX_dyn, tX_lum, subGalaxies, nest])
	

  
  # writing individual galaxies
  for galaxy in galList:
    


      if galaxy.inGroup <= 0. and galaxy.Vls<=3500:
	flag = galaxy.inGroup
	galaxy.flag = flag
        ID = galaxy.id  
        grp_id = galaxy.grp_id
        ra = galaxy.ra  
        dec = galaxy.dec
        Vls = galaxy.Vls  
        logK = galaxy.logK  
        Ks = galaxy.Ks
        dist = galaxy.dist
        
        
        if dist==0:
            distance = 0.00001
        else:
            distance = dist
        MagK = galaxy.Ks - 5*log10(distance) - 30 + 5
        
        mDist = galaxy.dist

        
        subGalaxies = galaxy.subGalaxies  
        R_theta = galaxy.R_theta  
        sigmaP_dyn = galaxy.sigma  
        nest = galaxy.id  
        Mv_lum = galaxy.M_v2
        R2t_lum = galaxy.R_2t2
        
        
        sigmaP_lum = (Mv_lum / (2.0E6))**(1./3)
        tX_lum = Mpc_km*R2t_lum*1.05*sqrt(1.5)/sigmaP_lum/sqrt(2.5)
        Rg = 0  # dynamical virial radius
        Mv_dyn = 0  # solar mass
        tX_dyn = 0
        R2t_dyn = 0
        Rg_angular = 0
        R_theta = 0 
       
        myTable.add_row([ID,grp_id,flag,ra,dec,MagK,Ks,logK,Vls, dist, \
	       mDist, R_theta, sigmaP_dyn, sigmaP_lum, \
	       Mv_dyn, Mv_lum, Rg_angular, Rg, R2t_dyn, R2t_lum, r2t_lum, tX_dyn, tX_lum, subGalaxies, nest])

	
  
  ID = 999999999; 
  ra = 999.9999; dec=-99.99;
  Ks= 99.99
  MagK = -99.99
  logK = 99.9999
  Vls = 9999
  dist = 99.99
  mDist = 99.99
  Mv_dyn = 9.99E99; Mv_lum = Mv_dyn
  tX_dyn = Mv_lum; tX_lum=Mv_lum
  nest = 9999999
  subGalaxies = 9999


  
  myTable.add_row([ID,grp_id,flag,ra,dec,MagK,Ks,logK,Vls, dist, \
	       mDist, R_theta, sigmaP_dyn, sigmaP_lum, \
	       Mv_dyn, Mv_lum, Rg_angular, Rg, R2t_dyn, R2t_lum, r2t_lum, tX_dyn, tX_lum, subGalaxies, nest])
  
  
  myTable.write(outfile, format='ascii.fixed_width',delimiter='|', bookend=False)
  
  # removing the last line, (it sits o adjust the column wodths)
  command =  ["csh", "remove_lastline.csh", outfile]
  subprocess.call(command)
################################################################# 
################################################################# 

def Theta_max(root):
  
  if root.left == None: return 0
  
  galList = NodeLeaf(root)
  
  N = len(galList)
  theta = np.zeros(N-1) 
  
  for i in range(1,N):
      
      theta[i-1] = angle(root.ra, root.dec, galList[i].ra, galList[i].dec)
  
  return np.max(theta)


################################################################# 
# gets a list of galaxies and returns Rg in terms of radian
def Rg_radian(galList):
  
  N = len(galList)
  sum = 0.
  BOL = True
  Rg = 0 
  n = 0 
  for i in range(0,N-1):
    for j in range(i+1,N):
      
      distij = angle(galList[i].ra, galList[i].dec, galList[j].ra, galList[j].dec)
      if distij !=0 : 
         sum += 1. / distij
         n+=1
      #else: BOL = False    # If one of the didtnaces is zero, therre is change to have duplicates 
      
  if sum != 0:
    Rg = n / sum
  
  return Rg


################################################################# 
 
def biweight_Rg(ra, dec, v, d):
   
   N = len(ra)
   dist_inverse = [] 
   for i in range(0,N-1):
    for j in range(i+1,N):
      distij = angle(ra[i], dec[i], ra[j], dec[j])
      if distij !=0 : 
         dist_inverse.append(  distij)
   
   dist_inverse = np.asarray(dist_inverse)
   n = len(dist_inverse)
   if n!=0:
     Med=np.median(dist_inverse)
     Rh_1 = Med
     #for p in range(0,10):
        #Rh_1 = st.biweight_location(dist_inverse, M=Med)
        #Med = Rh_1
     
     #Rg_radian =  N*N/(n*st.biweight_location(dist_inverse))
     #Rg_radian =  N*N/dist_inverse.sum()
     Rg_radian = Rh_1#2*N*Rh_1/(N-1)
     #print dist_inverse
     #print N, '  ' ,st.biweight_location(dist_inverse), '  ', dist_inverse.sum()/n
   else: 
     Rg_radian = float('nan')
   
   if d == 0 : d = v / H0
   return  d * Rg_radian
################################################################# 

        
#################################################################
def mergGroup(G_list, DDD=False):
    

    
    NoGroups = len(G_list)
 
    i = 0
    while i < NoGroups:
      j = i+1
      while j < NoGroups-1:
	
        n1 = G_list[i][0].subGalaxies
        n2 = G_list[j][0].subGalaxies

	d1 = G_list[i][0].dist # Vls/H0
	d2 = G_list[j][0].dist # Vls/H0
	#if d1 < 1 : d1 = 1 
	#if d2 < 1 : d2 = 1
	r1 = (180.*atan(G_list[i][0].R_2t2/d1)/pi)
	r2 = (180.*atan(G_list[j][0].R_2t2/d2)/pi)


	ang12 = (180./pi)*angle(G_list[i][0].ra, G_list[i][0].dec, G_list[j][0].ra, G_list[j][0].dec)
	
	
	d1 = G_list[i][0].mDist
	e1 = 0.05*d1

	d2 = G_list[j][0].mDist
	e2 = 0.05*d2

	delt = abs(d1-d2)
	
	v1 = G_list[i][0].Vls
	v2 = G_list[j][0].Vls

	
	
	# using dynamical velocity dispersions calculated based on the luminosities
        sig1 = (G_list[i][0].M_v2 / (2.0E6))**(1./3)
	sig2 = (G_list[j][0].M_v2 / (2.0E6))**(1./3)
	

	n1 = G_list[i][0].subGalaxies
	n2 = G_list[j][0].subGalaxies
	Bol = False
	

#######################################################################
        Sigquad = sqrt(sig1**2+sig2**2)
#######################################################################

        

	      
#######################################################################
      
	#if ang12 <= 1.1*(r1+r2) and max(r1,r2)<6:
	  #if abs(v1-v2) <= max(sig1,sig2) and (min(r1,r2))**3 < 0.2*(max(r1,r2))**3:
	       #Bol = True    
	       
	if ang12 <= 1.0*(r1+r2) and max(r1,r2)<6 and min(n1,n2)<5:
	  if abs(v1-v2) <= 2.0*max(sig1,sig2):
	      Bol = True
	      
	if ang12 <= 0.6*(r1+r2) and max(r1,r2)<6:
	  if abs(v1-v2) <= 2*Sigquad:
	      Bol = True

	      
	if ang12 <= 1.0*max(r1,r2) and max(r1,r2)<6:
	  if abs(v1-v2) <= max(2*Sigquad,3.0*max(sig1,sig2)):
	      Bol = True
	
	      
	if ang12 <= 1.0*max(r1,r2) and max(r1,r2)<6:
	  if d1!=0 and d2!=0 and delt < min(e1, e2):
	      Bol = True
	      
	
	  
	# one group completely projected on another one (e.g. smal in big)
	if ang12 <= 1.0*(max(r1,r2)-min(r1,r2)) and max(r1,r2)<6:
	  if abs(v1-v2) <= max(2*Sigquad,3.0*max(sig1,sig2)):
	      Bol = True	      
#######################################################################
        grp1 = G_list[i][0]
        grp2 = G_list[j][0]
        if DDD==True:
            Bol = False
            separation = distance(grp1.ra,grp1.dec,grp1.dist,grp2.ra,grp2.dec,grp2.dist)
            
        if DDD==True and separation<1.00*(grp1.r2t+grp2.r2t):
            Bol = True

#######################################################################
	  
        if Bol:  # merging
	  
	  newGroup = []
	  
	  for p in range(1, len(G_list[i])):
	    newGroup.append(G_list[i][p])
	  
	  for q in range(1, len(G_list[j])):
	    newGroup.append(G_list[j][q])
	  
 
	  root = LoosgalListJoint(newGroup, grID = 200000000)
	  G_list[i] = NodeLeaf(root)

	  
	
	  G_list.pop(j) 
	  NoGroups-=1
	  i=0
	  j=0
        j+=1
      i+=1
    
    return G_list

 #################################################################
 # Function name: addGalGroup
 #
 # This function tries to find those galaxies which do not fall 
 # into any Group. If the galaxy is close enough to the center of any group,
 # according to its radial velocity (i.e. within 2*sigma of the group radial velocity) 
 # it might get absorbed by that group...
 #################################################################

def addGalGroup(G_list, galList, DDD=False):
  

    singles = []
    for galaxy in galList:
       if galaxy.inGroup == 0:
	  singles.append([galaxy, -1, -10000])  # [galaxy, new.group.id, angular.separation]
    
    N = len(G_list)

    
    for entity in singles:
            ra = entity[0].ra
            dec = entity[0].dec
    
	    for i in range(len(G_list)): # group_indices[0]:
    
		    ang12 = (180./pi)*angle(G_list[i][0].ra, G_list[i][0].dec, ra, dec)
		    d = G_list[i][0].dist # Vls / H0
		    #if d < 1: d = 1 
		    r = (180.*atan(G_list[i][0].R_2t2/(d))/pi)	      
		    

		    sig_p = (G_list[i][0].M_v2 / (2.0E6))**(1./3)
		    

		    if ang12 <= 2*r or DDD==True:

                      join = False
                      
		      if DDD==True or (ang12 <= 1.01*(180./pi)*atan(1.0*G_list[i][0].R_2t2/d) and abs(entity[0].Vls - G_list[i][0].Vls) < 2.0*sig_p):
			  d1 = G_list[i][0].mDist
			  e1 = d1 * 0.05
			  d2 = entity[0].dist
			  e2 = d2 * 0.05
			  delt = abs(d1-d2)
			  
  
			  if ang12 > 1.01*(180./pi)*atan(1.0*G_list[i][0].R_2t2/d) and d1!=0 and d2!=0 and delt<max(e1,e2):
			    join = True
			  if d1==0 or d2==0:
			    join = True
			  if ang12 <= 1.01*(180./pi)*atan(1.0*G_list[i][0].R_2t2/d):
			    join = True
			  
			  if DDD==True: join = False
			  if DDD==True and distance(G_list[i][0].ra,G_list[i][0].dec,G_list[i][0].dist,entity[0].ra,entity[0].dec,entity[0].dist)<1.00*G_list[i][0].r2t:
                              join=True
                              
                              
                              
			  if join==True and  sig_p > entity[2]:
			    entity[1] = i
			    entity[2] = sig_p
		  
            if entity[1] > -1: # if the galaxy is absorbed (check the id of the target group)
              entity[0].inGroup = 1
              G_list[entity[1]].append(entity[0])
    for i in range(N):
       G_list[i].pop(0)
       root  = LoosgalListJoint(G_list[i], grID = 300000000)
       G_list[i] = NodeLeaf(root) 
    return G_list


#################################################################
def Remove_Patch(G_list, galList, RA=[0,10], DEC=[0,10], RAA=False):
    
    patch_galList = []
    
    N_gals = len(galList)
    N_grps = len(G_list)
    
    i = 0
    while i<N_gals:
        galaxy = galList[i]
        
        Bol = False
        if RAA and (galaxy.ra>=RA[0] or galaxy.ra<=RA[1]) and galaxy.dec>=DEC[0] and galaxy.dec<DEC[1]:
            Bol = True
        
        if Bol or (galaxy.ra>=RA[0] and galaxy.ra<RA[1] and galaxy.dec>=DEC[0] and galaxy.dec<DEC[1]):
            if galaxy.inGroup == 0:
                galList.pop(i)
                patch_galList.append(galaxy)
                i-=1
                N_gals-=1
        i+=1    
        
    remove_ids = [] 
    i = 0   
    while i<N_grps:
        
        
        group = G_list[i]
        Bol = False
        if RAA and (group[0].ra>=RA[0] or group[0].ra<=RA[1]) and group[0].dec>=DEC[0] and group[0].dec<DEC[1]:
            Bol = True


        if Bol or (group[0].ra>=RA[0] and group[0].ra<RA[1] and group[0].dec>=DEC[0] and group[0].dec<DEC[1]):
            G_list.pop(i)
            i-=1
            N_grps-=1
            
            for galaxy in group[1:]:
                remove_ids.append(galaxy.id)
        i+=1
                
                
    N_gals = len(galList)
    i = 0
    while i<N_gals:
        galaxy = galList[i]
        if galaxy.id in remove_ids:
            patch_galList.append(galList.pop(i))
            i-=1
            N_gals-=1
        i+=1             
    
    for galaxy in patch_galList:
        galaxy.inGroup = 0
        galaxy.nest = galaxy.id
        galaxy.mDist = galaxy.dist
            
    return patch_galList

#################################################################
# Sometimes when groups grow the order of galaxy absorbation must be different:
def IndividualGal_modifier(G_list, galList, DDD=False):

  tmp_galList = []
  N = len(galList)
  for galaxy in galList:
    if galaxy.inGroup > 0:
      galaxy.inGroup = 0
      tmp_galList.append(galaxy)
    

  new_G_list=[]
  for Group in G_list:
    Group[1].inGroup = 1  # first galaxy of group, dominant one
    new_G_list.append([Group[0], Group[1]])
  

  Lum_mass = []
  for i in range(0, len(new_G_list)): Lum_mass.append(new_G_list[i][0].M_v2) # look at the mass of each group
  Lum_mass = np.asarray(Lum_mass)
  indices = np.argsort(Lum_mass)
  
  NEWgpLIST = []
  for i in indices[::-1]:
            
            
            grp = new_G_list[i][0]
            ra = grp.ra
            dec = grp.dec
            r2t = grp.r2t
            dist = grp.dist 
    

	    j = 0 
	    N = len(tmp_galList)
	    while j < N: 
                     
                  gal = tmp_galList[j]

	          ang12 = (180./pi)*angle(ra, dec, gal.ra, gal.dec)
                  d = grp.dist # Vls / H0
                  #if d < 1: d = 1 
                  r = (180.*atan(0.95*grp.R_2t2/(d))/pi)   
                  sig_p = (grp.M_v2 / (2.0E6))**(1./3)
                  
                  join = False
                  
                  if DDD==False and ang12 <= r and gal.inGroup == 0 and abs(grp.Vls-gal.Vls) <= 2.0*sig_p:
                      join = True
                  
                  if DDD==True and gal.inGroup == 0 and distance(ra, dec, dist, gal.ra, gal.dec, gal.dist)<1.05*r2t:
                      join = True
                      
                  
                  
                  if join:
			new_G_list[i].append(gal)
			gal.inGroup = 1
		   	tmp_galList.pop(j)
			j -= 1
			N -= 1
	          j += 1

    
	    new_G_list[i].pop(0)
	    if len(new_G_list[i]) > 1: 
	        grp =  LoosgalListJoint(new_G_list[i], grID = 400000000)
	        NEWgpLIST.append(NodeLeaf(grp))
	    else:
	        new_G_list[i][0].inGroup = 0
       
  return NEWgpLIST


#################################################################
## Trying to find linkages when both galaxies that have distances

def find_pair_dist(galList_org, DDD=False):
  
  galList = []
  for galaxy in galList_org:
    if galaxy.inGroup == 0 and galaxy.Vls<=4000: 
	galList.append(galaxy)  
  
  individuals = []
  Lum_mass = []
  for i in range(0, len(galList)): Lum_mass.append(galList[i].M_v2)
  Lum_mass = np.asarray(Lum_mass)
  indices = np.argsort(Lum_mass)
  for i in indices[::-1]: individuals.append(galList[i])  
  
  NEWgpLIST =[]
  p = 0 
  N = len(individuals)
  while p < N-1:
    #print p, N
    galaxies = []
    galaxies.append(individuals[p])

    q = (p+1)%N
    pivot = p
    grp = galaxies[0]
    Bol2 = False
    
    if grp.dist==0:
      p+=1
      continue
    
    while q!=pivot:
      
      if individuals[q].dist==0:
	q = (q+1)%N
	continue
      
      ang12 = (180./pi)*angle(grp.ra, grp.dec, individuals[q].ra, individuals[q].dec)
      
      d = grp.dist
      if d == 0:
	q = (q+1)%N
	continue
      
      
      mass = Mass(10**m_logK(grp.MagK))
      R_2t2 = 0.215*((mass/1.E12)**(1./3))

      thet = (180.*atan(1.0*R_2t2/d)/pi)

      coeff = 1.3
      
      Bol = False

      if DDD==False and abs(individuals[q].dist - d)/d < 0.1 and ang12 < coeff*thet :
          Bol =True
      
      if DDD==True and distance(grp.ra,grp.dec,grp.dist,individuals[q].ra,individuals[q].dec,individuals[q].dist) < 1.05*grp.r2t: 
          Bol = True
      
      
      if Bol:
            
            mass = Mass(10**m_logK(individuals[q].MagK))
            R_2t2_2 = 0.215*((mass/1.E12)**(1./3))
            
            Bol2 = True

	    galaxies.append(individuals[q])
	    grp =  LoosgalListJoint(galaxies, grID = 400000000)
	    individuals.pop(q)

	    N-=1
	    
	    if p>=q:  p-=1
	    
	    
	    if len(galaxies) == 2:
	      individuals.pop(p)
	      N-=1	
	      if q>p: q = (q-1)%N
	      break
            
            
	    q = (q-1)%N
	    pivot = q 
	    
      q = (q+1)%N
    if Bol2: 
      NEWgpLIST.append(NodeLeaf(grp))
    p+=1
    
    for i in range(0, len(NEWgpLIST)):
      for j in range(1, len(NEWgpLIST[i])):
         NEWgpLIST[i][j].inGroup = 1
  
  
  return NEWgpLIST  


#################################################################

def IndividualGal(galList_org, pairs=False, DDD=False):
  
  
  galList = []
  for galaxy in galList_org:
    if galaxy.inGroup == 0 and galaxy.Vls<=4000: 
	galList.append(galaxy)   
  

  individuals = []
  Lum_mass = []
  for i in range(0, len(galList)): Lum_mass.append(galList[i].M_v2)
  Lum_mass = np.asarray(Lum_mass)
  indices = np.argsort(Lum_mass)
  for i in indices[::-1]: individuals.append(galList[i])  
  
  NEWgpLIST =[]
  p = 0 
  N = len(individuals)
  while p < N-1:
    #print p, N
    galaxies = []
    galaxies.append(individuals[p])
    #individuals.pop(p)
    #N-=1
    q = (p+1)%N
    pivot = p
    grp = galaxies[0]
    Bol2 = False
    while q!=pivot:
      
      ang12 = (180./pi)*angle(grp.ra, grp.dec, individuals[q].ra, individuals[q].dec)
      
      d = grp.mDist # Vls/H0
      #if d <=1: d =1

      thet = (180.*atan(1.0*grp.R_2t2/d)/pi)

      
      coeff = 1.0 # 1.05
      delta_v = abs(grp.Vls-individuals[q].Vls)
      sig_p = (grp.M_v2 / (2.0E6))**(1./3)
      
      test = False
      if pairs:
	sig_sig = 100000
	if len(galaxies) == 1:
	  d = individuals[q].dist # Vls/H0
	  #if d <=1: d =1
	  thet +=  (180.*atan(individuals[q].R_2t2/d)/pi)
	  L1 = 10**grp.logK
	  L2 = 10**individuals[q].logK
	  L_tot = L1 + L2
	  mass = Mass(L_tot)
	  sig_sig = (mass / (2.0E6))**(1./3)
	  if delta_v < 2.0*sig_sig: 
	    test = True
	
      
      
      if True:
	
        Bol = False
	if DDD==False and ang12 <= coeff*thet and (test==True or delta_v <= 2.0*sig_p):
            Bol = True
        
        if DDD==True and distance(grp.ra,grp.dec,grp.dist,individuals[q].ra,individuals[q].dec,individuals[q].dist) < grp.r2t: 
            Bol = True
        
            
        
        
        
        if Bol:

	    Bol2 = True
	    galaxies.append(individuals[q])
	    grp =  LoosgalListJoint(galaxies, grID = 400000000)
	    individuals.pop(q)

	    N-=1
	    
	    if p>=q:  p-=1
	    
	    
	    if len(galaxies) == 2:
	      individuals.pop(p)
	      N-=1	
	      if q>p: q = (q-1)%N
            
            
	    q = (q-1)%N
	    pivot = q 
	    
      q = (q+1)%N
    if Bol2: 
      NEWgpLIST.append(NodeLeaf(grp))
    p+=1
    
    for i in range(0, len(NEWgpLIST)):
      for j in range(1, len(NEWgpLIST[i])):
         NEWgpLIST[i][j].inGroup = 1
   
  return NEWgpLIST

#################################################################

def LoosgalListJoint(galList, grID = 500000000):
   
   if len(galList)==0: return None
   
   Lum_mass = []
   for i in range(len(galList)):
       Lum_mass.append(galList[i].logK)
   Lum_mass = np.asarray(Lum_mass)
   indices = np.argsort(Lum_mass, kind='mergesort')
   indices = indices[::-1]
       
   new_galList = []
   for i in indices:
      new_galList.append(galList[i])
   
   galList = new_galList 
   
   root = None
   for i in indices:
       root = galJoint(root, galList[i], grID + galList[i].id)
   
   sumDist = 0
   V_mean = 0 

   for galaxy in galList:
    sumDist += galaxy.sumDist
    V_mean += galaxy.Vls
   
   V_mean /= len(galList)
   mDist = sumDist / len(galList)

   for gal in galList:
      gal.logK = m_logK(gal.MagK)
      gal.mDist = mDist
      
   root = None
   for i in indices:
      root = galJoint(root, galList[i], grID + galList[i].id)   
   
   
   
   return root

#################################################################

def Theta_gr(head, galList):
  
  
  
  N = len(galList)
  theta = np.zeros(N-1) 
  
  for i in range(1,N):
      
      theta[i-1] = angle(head.l, head.b, galList[i].l, galList[i].b)
  
  return np.max(theta)




#################################################################
def extractPGC(id, grp=False, supergrp=False):
  
  if not grp and not supergrp:
    return id
  
  
  if grp:
    pgc = int(id)%100000000
  
  
  if supergrp:
    grp = int(id)%10000000000
    pgc = int(grp)%100000000
  
  return pgc
#################################################################

def group_moderator(G_list):
  
  new_G_list = []
  for group in G_list:
      
       Lum_mass = []
       for i in range(len(group)):
           Lum_mass.append(group[i].logK)
       Lum_mass = np.asarray(Lum_mass)
       indices = np.argsort(Lum_mass, kind='mergesort')
       indices = indices[::-1]
       
       new_group = []
       for i in indices:
           new_group.append(group[i])

       new_G_list.append(new_group)
  
  G_list = new_G_list
  
  for group in new_G_list:
	  groupHeader = group[0]
	  groupHeader.flag = 2
	  
	  
	  groupHeader.id = group[1].id + 900000000
	  groupHeader.nest = group[1].id

          sum_v = 0.
          sum_v2 = 0.
       
          n_v = 0
          for galaxy in group[1:]:
              
            galaxy.nest = group[1].id

	    if galaxy.Vls !=0:
	      sum_v += galaxy.Vls
	      sum_v2 += galaxy.Vls**2
	      n_v += 1
	  
	  if n_v != 0:
	    v_av = sum_v / n_v
	    
	    groupHeader.Vls = v_av
	    groupHeader.sigma = sqrt(sum_v2/n_v - v_av**2) 
	    

	  else:
	    groupHeader.Vls = 0
	    groupHeader.sigma = 0  # in Vls
	  
	  
	  meanDist = 0.
	  sumDist = 0
	  N= 0 
	  for galaxy in group[1:]:
	    if galaxy.dist != 0:
	            sumDist += galaxy.dist
		    N += 1
	  
  
          meanDist = sumDist/N

          groupHeader.mDist = meanDist
          groupHeader.dist  = meanDist

	  L_tot = 0

	  for galaxy in group[1:]:
	    galaxy.flag = 1
	    galaxy.setMeanDist(meanDist, GRP_vel = groupHeader.Vls, dist=galaxy.dist)
	    L_tot += 10**galaxy.logK

          groupHeader.logK = log10(L_tot)
          Dist_v =  groupHeader.dist # Vls / H0
          #if Dist_v<1. : Dist_v=1
          Mk_sun = 3.28   # K-band
          M = Mk_sun - (groupHeader.logK / 0.4)
          groupHeader.Ks = M + 5*log10(Dist_v) + 30 - 5
          
          groupHeader.M_v2 = Mass(10**groupHeader.logK)
          if np.isinf(groupHeader.M_v2): 
	    print '[ERROR] Modification ...'
	    print Mass(10**10.6453)
	    print 10**groupHeader.logK
	    sys.exit()
          groupHeader.R_2t2 = 0.215*((groupHeader.M_v2/1.E12)**(1./3))  # Mpc
          groupHeader.r2t = sqrt(1.5)*groupHeader.R_2t2
          groupHeader.r1t = 3.5*groupHeader.r2t
  
  return new_G_list

#################################################################
def Make_Group(galList, DDD=False):
    
        print 
        print "[NEW] working on grouping ... "+str(len(galList))+' glaxies'
        
        if len(galList)==0:
            print "[Warning] No Group ... "+str(len(galList))+' glaxies'
            return []
        
	##  Finding groups based on the criteria
        G_list = []
        
        for qp in range(3):
	   G_list += IndividualGal(galList, DDD=DDD)
	   G_list += find_pair_dist(galList, DDD=DDD)
	
	
	print "working on addGalGroup"
	#### fixed
        G_list = addGalGroup(G_list, galList, DDD=DDD)
        
        G_list = group_moderator(G_list)

        
        
        if True:
          print "working on IndividualGal_modifier"
          for pq in range(3): 
	    G_list = IndividualGal_modifier(G_list, galList, DDD=DDD)
	    G_list += find_pair_dist(galList)
	  
	  
	  print "working on addGalGroup2"  
	  for pq in range(3): 
	    G_list = addGalGroup(G_list, galList, DDD=DDD)
	  
	  print "working on merge"  
          for qp in range(10):
            mergGroup(G_list, DDD=DDD)
          
          G_list = group_moderator(G_list)
          mergGroup(G_list, DDD=DDD)
          
          print "working on addGalGroup3" 
	  for pq in range(3): 
	    G_list = addGalGroup(G_list, galList, DDD=DDD)  
	  
	  print "working on merge2"
          for qp in range(3):
            mergGroup(G_list, DDD=DDD)
                       
          G_list = group_moderator(G_list)
        
        print 'No. of groups: ', len(G_list)
        
        return G_list
    
#################################################################
def Equatorial_Patch(ra_del, dec, Groups, Gals):
    
  for ra_i in np.arange(20,360-ra_del,ra_del):
      ra_f = ra_i+ra_del
      print '.......................................'
      print "Equatorial Ring RA,DEC= ", [ra_i, ra_f], [dec-10,dec+10]
      patch_galList = Remove_Patch(Groups, Gals, RA=[ra_i, ra_f], DEC=[dec-10,dec+10])
      Groups += Make_Group(patch_galList)
      Gals   += patch_galList
  
  print '.......................................'
  print "Equatorial Ring RA,DEC= ", [ra_f, 20], [dec-10,dec+10]
  patch_galList = Remove_Patch(Groups, Gals, RA=[ra_f, 20], DEC=[dec-10,dec+10], RAA=True)
  Groups += Make_Group(patch_galList)
  Gals   += patch_galList          
#################################################################
def Long_Patch(ra_del, ra_overlap, Groups, Gals, DEC=[0,90]):
              
  for ra_i in np.arange(ra_del,360,ra_del):
    print '.......................................'
    print "Revisiting Patch .... RA,DEC = ", [ra_i-ra_overlap,ra_i+ra_overlap], DEC
    patch_galList = Remove_Patch(Groups, Gals, RA=[ra_i-ra_overlap,ra_i+ra_overlap], DEC=DEC)
    Groups += Make_Group(patch_galList)
    Gals   += patch_galList


  print '.......................................'
  print "Revisiting Patch .... RA,DEC = ", [360-ra_overlap,0+ra_overlap], DEC
  patch_galList = Remove_Patch(Groups, Gals, RA=[360-ra_overlap,0+ra_overlap], DEC=DEC, RAA=True)
  Groups += Make_Group(patch_galList)
  Gals   += patch_galList    

  print len(Groups), len(Gals)
#################################################################

if __name__ == '__main__':
  
  if len(sys.argv) != 2:
      print 'Give an argument (2D or 3D): '
      sys.exit()      
  if sys.argv[1] == '2D':
      _3D = False
  elif sys.argv[1] == '3D':
      _3D = True
  else:
      print 'Give an argument (2D or 3D): '
      sys.exit()
        
  # Originally made by M_L_ratio_curve_v2.py
  table = np.genfromtxt('M_L_curve_v2.csv' , delimiter=',', filling_values=0, names=True, dtype=None)
  Lumin_in   = 10**table['log10_L']
  Mass_out   = 10**table['log10_M']
  
  LtoM_func = interpolate.interp1d(Lumin_in, Mass_out)
  
  inFile = 'simple_mock-14.0.csv'
  #inFile = 'simple_mock.csv'
  
  table = np.genfromtxt(inFile , delimiter=',', filling_values=0, names=True, dtype=None)
  
  
  d_min = 5.
  d_max = 100.

  ra_del = 30.
  ra_overlap = 10
  
  Dec_lst = [-90,-60,-25,0,25,60,90]
  #Dec_lst = [-90,0,90]
  #############################################################
  galLists = []
  
  for ra_i in np.arange(0,360,ra_del):
    ra_f = ra_i+ra_del
    for i in range(len(Dec_lst)-1):
       galLists.append(readgalList(table, VLS=[-4000,4000], d=[d_min, d_max], MagK_limit=-16, RA=[ra_i, ra_f], DEC=[Dec_lst[i],Dec_lst[i+1]]))
  
  #galLists.append(readgalList(table, VLS=[-4000,4000], d=[d_min, d_max], MagK_limit=-16, RA=[100,102], DEC=[-47,-45]))
  #############################################################
  
  Gals = []
  for gals in galLists:
      Gals += gals 
  print  ' Total # of galaxies: ', len(Gals)
  
  
  G_lists = []
  N = len(galLists)
  for i in range(N):
      galList = galLists[i]
      print
      print "Patch # "+ str(i+1)+'/'+str(N)
      G_list = Make_Group(galList)
      G_lists.append(G_list)

  
########################################################################
        
  Groups = []
  for Glist in G_lists:
      Groups += Glist
  
  Gals = []
  for gals in galLists:
      Gals += gals
  
########################################################################

  Long_Patch(ra_del, ra_overlap, Groups, Gals, DEC=[0,90])
  Long_Patch(ra_del, ra_overlap, Groups, Gals, DEC=[-90,0])
  
  for i in range(1,len(Dec_lst)-1):
     Equatorial_Patch(ra_del, Dec_lst[i], Groups, Gals)


          

  print '.......................................'
  print "North Cap "
  patch_galList = Remove_Patch(Groups, Gals, RA=[0,360], DEC=[60,90])
  Groups += Make_Group(patch_galList)
  Gals   += patch_galList


  print '.......................................'
  print "South Cap "
  patch_galList = Remove_Patch(Groups, Gals, RA=[0,360], DEC=[-90,-60])
  Groups += Make_Group(patch_galList)
  Gals   += patch_galList

          
#####################################################################          
          
  print "# of Groups, Gals", len(Groups), len(Gals)   
### MErging Sub Group-lists  #####################################################################


  print "working on FINAL merge ... "
  for qp in range(10):
      mergGroup(Groups, DDD=_3D)
  print "working on addGalGroup3" 
  for pq in range(3): 
      Groups = addGalGroup(Groups, Gals, DDD=_3D)    
  Groups = group_moderator(Groups)
  
  print "# of Groups, Gals", len(Groups), len(Gals) 
#####################################################################


  print "working on writing ..."
  outFile = 'simul.MagK16.d5.d100.'+sys.argv[1]+'.v01.group'
  #outFile = 'test.group'
  groupWrite(outFile, Groups, Gals)
  print "Created: "+outFile
  print 

        

  
