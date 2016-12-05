#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 13 09:45:04 2016

This code allows me to calculate all combinations
of variable leaf parameters. I will input these
combinations of parameters into my photosynthesis
model in order to perform a sensitivity analysis
of all variables. The combinations of parameters
allows me to test covariance (dependencies) among
parameters. 

Note: Using all combinations of parameters
generates WAY too many plots. I want to do
two combinations at a time (with all other
parameters constant). In other words,
I want to try (A) 2 increasing parameters,
(B) 2 decreasing parameters, and (C) 1 parameter
that increases and another that decreases and 
vice versa. 

@author: Katherine
"""

##CHECK WUE FOR VARIABLE IJ--WHY DOES WUE DECREASE???--might want to mess around with range values (0.1-0.5 for example)##

import numpy as np
import itertools as it

#---------------Functions---------------#   

def merge_dicts(*dict_args):
    '''
    Given any number of dicts, shallow copy and merge into a new dict,
    precedence goes to key value pairs in latter dicts.
    '''
    result = {}
    for dict_int in dict_args:
        for dictionary in dict_int:
            result.update(dictionary)
    return result

    
    
#---------------Create Dictionary For All Variable Parameter Inputs---------------#  

##mean trait assemblage parameters for dry meadow, moist meadow, wet meadow, snowbed 

#fellfield
s_f=0.015 #specific leaf area (m2/g)
nm_f=0.0253  #leaf N (gN/gC)
flnr_f=0.55 #fraction of leaf nitrogen in rubisco (g leaf N/g rubisco)
chl_f=200 #chlorophyll content (umol chl/m2)
tl_f=36.0 #temp of leaf (K)-->f(leaf area,angle,height)
vwc_f=0.10 #soil volumetric water content (cm3/cm3)
ij_f=0.92 #index of light penetration (low leaf area and high leaf angle) (unitless)-->f(leaf area & angle)

#dry meadow
s_d=0.017 #specific leaf area (m2/g)
nm_d=0.0313 #leaf N (gN/gC)
flnr_d=0.6 #fraction of leaf nitrogen in rubisco (g leaf N/g rubisco)
chl_d=300 #chlorophyll content (umol chl/m2)
tl_d=32.0 #temp of leaf (K)-->f(leaf area,angle,height)
vwc_d=0.16  #soil volumetric water content (cm3/cm3)
ij_d=0.94 #index of light penetration (low medium leaf area and high medium leaf angle) (unitless)-->f(leaf area & angle)

#moist meadow
s_m=0.019 #specific leaf area (m2/g)
nm_m=0.0372 #leaf N (gN/gC)
flnr_m=0.65 #fraction of leaf nitrogen in rubisco (g leaf N/g rubisco)
chl_m=400 #chlorophyll content (umol chl/m2)
tl_m=28.0 #temp of leaf (K)-->f(leaf area,angle,height)
vwc_m=0.22 #soil volumetric water content (cm3/cm3)
ij_m=0.96 #index of light penetration (medium leaf area and medium leaf angle) (unitless)-->f(leaf area & angle)

#wet meadow
s_w=0.022 #specific leaf area (m2/g)
nm_w=0.0435#leaf N (gN/gC)
flnr_w=0.7 #fraction of leaf nitrogen in rubisco (g leaf N/g rubisco)
chl_w=500 #chlorophyll content (umol chl/m2)
tl_w=24.0 #temp of leaf (K)-->f(leaf area,angle,height)
vwc_w=0.28 #soil volumetric water content (cm3/cm3)
ij_w=0.98 #index of light penetration (high medium leaf area and low medium leaf angle) (unitless)-->f(leaf area & angle)

#snowbed
s_s=0.025 #specific leaf area (m2/g)
nm_s=0.0495 #leaf N (gN/gC)
flnr_s=0.75 #fraction of leaf nitrogen in rubisco (g leaf N/g rubisco)
chl_s=600 #chlorophyll content (umol chl/m2)
tl_s=20.0 #temp of leaf (K)-->f(leaf area,angle,height)
vwc_s=0.34 #soil volumetric water content (cm3/cm3)
ij_s=1.0 #index of light penetration (high leaf area and low leaf angle) (unitless)-->f(leaf area & angle)

##lambda (wue)
#lambda_h=10000.0*(1.0/10.0**6) #high (umol CO2/umol H20)
#lambda_hm=7600.0*(1.0/10.0**6) #high medium (umol CO2/umol H20)
#lambda_m=5200.0*(1.0/10.0**6) #medium (umol CO2/umol H20)
#lambda_lm=2800.0*(1.0/10.0**6) #low medium (umol CO2/umol H20)
#lambda_l=400.0*(1.0/10.0**6) #low (umol CO2/umol H20)

#ball berry parameter (sensitivity of stomatal conductance)
m_f=24.0 #(unitless)
m_d=20.0 #(unitless)
m_m=16.0 #(unitless)
m_w=12.0 #(unitless)
m_s=8.0 #(unitless)

std_s=np.std([s_f,s_d,s_m,s_w,s_s])*0.1
std_nm=np.std([nm_f,nm_d,nm_m,nm_w,nm_s])*0.2
std_chl=np.std([chl_f,chl_d,chl_m,chl_w,chl_s])*0.1
std_tl=np.std([tl_f,tl_d,tl_m,tl_w,tl_s])*0.1
std_m=np.std([m_f,m_d,m_m,m_w,m_s])*0.1
std_ij=np.std([ij_f,ij_d,ij_m,ij_w,ij_s])*0.1

mean_s=np.mean([s_f,s_d,s_m,s_w,s_s])
mean_nm=np.mean([nm_f,nm_d,nm_m,nm_w,nm_s])
mean_chl=np.mean([chl_f,chl_d,chl_m,chl_w,chl_s])
mean_tl=np.mean([tl_f,tl_d,tl_m,tl_w,tl_s])
mean_m=np.mean([m_f,m_d,m_m,m_w,m_s])
mean_ij=np.mean([ij_f,ij_d,ij_m,ij_w,ij_s])

s_range=s_m-min(np.random.normal(mean_s,std_s,100000))
nm_range=nm_m-min(np.random.normal(mean_nm,std_nm,100000))
chl_range=chl_m-min(np.random.normal(mean_chl,std_chl,100000))
tl_range=tl_m-min(np.random.normal(mean_tl,std_tl,100000))
m_range=m_m-min(np.random.normal(mean_m,std_m,100000))
ij_range=ij_m-min(np.random.normal(mean_ij,std_ij,100000))
vwc_range=0.0

    

###plasticity of traits
#s_range=np.abs(s_s-s_w)*0.1
#nm_range=np.abs(nm_s-nm_w)*0.1
#chl_range=np.abs(chl_s-chl_w)*0.1
#tl_range=np.abs(tl_s-tl_w)*0.1
#vwc_range=np.abs(vwc_s-vwc_w)*0.0 #assume no range in this value (not a trait)
#m_range=np.abs(m_s-m_w)*0.1
#ij_range=np.abs(ij_s-ij_w)*0.1

##---In the Below Dictionary, sla, leaf n, and fraction of leaf n in rub are the only proportionally dependent variables---##

#variable s
#dict=[
#      [{'s_f':np.linspace(s_f-s_range,s_f+s_range,2)},
#       {'nm_f':np.linspace(nm_m-nm_range,nm_m+nm_range,2)},
#       {'m_f1':np.zeros(shape=2)+m_m-m_range},
#       {'m_f2':np.zeros(shape=2)+m_m+m_range},
#       {'chl_f1':np.zeros(shape=2)+chl_m-chl_range},
#       {'chl_f2':np.zeros(shape=2)+chl_m+chl_range},
#       {'tl_f1':np.zeros(shape=2)+tl_m-tl_range},
#       {'tl_f2':np.zeros(shape=2)+tl_m+tl_range},
#       {'vwc_f1':np.zeros(shape=2)+vwc_m-vwc_range},
#       {'vwc_f2':np.zeros(shape=2)+vwc_m+vwc_range},
#       {'ij_f1':np.zeros(shape=2)+ij_m-ij_range},
#       {'ij_f2':np.zeros(shape=2)+ij_m+ij_range},
#        ],
#      
#      [{'s_d':np.linspace(s_d-s_range,s_d+s_range,2)},
#       {'nm_d':np.linspace(nm_m-nm_range,nm_m+nm_range,2)},
#       {'m_d1':np.zeros(shape=2)+m_m-m_range},
#       {'m_d2':np.zeros(shape=2)+m_m+m_range},
#       {'chl_d1':np.zeros(shape=2)+chl_m-chl_range},
#       {'chl_d2':np.zeros(shape=2)+chl_m+chl_range},
#       {'tl_d1':np.zeros(shape=2)+tl_m-tl_range},
#       {'tl_d2':np.zeros(shape=2)+tl_m+tl_range},
#       {'vwc_d1':np.zeros(shape=2)+vwc_m-vwc_range},
#       {'vwc_d2':np.zeros(shape=2)+vwc_m+vwc_range},
#       {'ij_d1':np.zeros(shape=2)+ij_m-ij_range},
#       {'ij_d2':np.zeros(shape=2)+ij_m+ij_range},
#        ],
#        
#      [{'s_m':np.linspace(s_m-s_range,s_m+s_range,2)},
#       {'nm_m':np.linspace(nm_m-nm_range,nm_m+nm_range,2)},
#       {'m_m1':np.zeros(shape=2)+m_m-m_range},
#       {'m_m2':np.zeros(shape=2)+m_m+m_range},
#       {'chl_m1':np.zeros(shape=2)+chl_m-chl_range},
#       {'chl_m2':np.zeros(shape=2)+chl_m+chl_range},
#       {'tl_m1':np.zeros(shape=2)+tl_m-tl_range},
#       {'tl_m2':np.zeros(shape=2)+tl_m+tl_range},
#       {'vwc_m1':np.zeros(shape=2)+vwc_m-vwc_range},
#       {'vwc_m2':np.zeros(shape=2)+vwc_m+vwc_range},
#       {'ij_m1':np.zeros(shape=2)+ij_m-ij_range},
#       {'ij_m2':np.zeros(shape=2)+ij_m+ij_range},
#        ],
#       
#      [{'s_w':np.linspace(s_w-s_range,s_w+s_range,2)},
#       {'nm_w':np.linspace(nm_m-nm_range,nm_m+nm_range,2)},
#       {'m_w1':np.zeros(shape=2)+m_m-m_range},
#       {'m_w2':np.zeros(shape=2)+m_m+m_range},
#       {'chl_w1':np.zeros(shape=2)+chl_m-chl_range},
#       {'chl_w2':np.zeros(shape=2)+chl_m+chl_range},
#       {'tl_w1':np.zeros(shape=2)+tl_m-tl_range},
#       {'tl_w2':np.zeros(shape=2)+tl_m+tl_range},
#       {'vwc_w1':np.zeros(shape=2)+vwc_m-vwc_range},
#       {'vwc_w2':np.zeros(shape=2)+vwc_m+vwc_range},
#       {'ij_w1':np.zeros(shape=2)+ij_m-ij_range},
#       {'ij_w2':np.zeros(shape=2)+ij_m+ij_range},
#        ], 
#      
#      [{'s_s':np.linspace(s_s-s_range,s_s+s_range,2)},
#       {'nm_s':np.linspace(nm_m-nm_range,nm_m+nm_range,2)},
#       {'m_s1':np.zeros(shape=2)+m_m-m_range},
#       {'m_s2':np.zeros(shape=2)+m_m+m_range},
#       {'chl_s1':np.zeros(shape=2)+chl_m-chl_range},
#       {'chl_s2':np.zeros(shape=2)+chl_m+chl_range},
#       {'tl_s1':np.zeros(shape=2)+tl_m-tl_range},
#       {'tl_s2':np.zeros(shape=2)+tl_m+tl_range},
#       {'vwc_s1':np.zeros(shape=2)+vwc_m-vwc_range},
#       {'vwc_s2':np.zeros(shape=2)+vwc_m+vwc_range},
#       {'ij_s1':np.zeros(shape=2)+ij_m-ij_range},
#       {'ij_s2':np.zeros(shape=2)+ij_m+ij_range},
#        ]]


#variable nm
#dict=[
#      [{'s_f':np.linspace(s_m-s_range,s_m+s_range,2)},
#       {'nm_f':np.linspace(nm_f-nm_range,nm_f+nm_range,2)},
#       {'m_f1':np.zeros(shape=2)+m_m-m_range},
#       {'m_f2':np.zeros(shape=2)+m_m+m_range},
#       {'chl_f1':np.zeros(shape=2)+chl_m-chl_range},
#       {'chl_f2':np.zeros(shape=2)+chl_m+chl_range},
#       {'tl_f1':np.zeros(shape=2)+tl_m-tl_range},
#       {'tl_f2':np.zeros(shape=2)+tl_m+tl_range},
#       {'vwc_f1':np.zeros(shape=2)+vwc_m-vwc_range},
#       {'vwc_f2':np.zeros(shape=2)+vwc_m+vwc_range},
#       {'ij_f1':np.zeros(shape=2)+ij_m-ij_range},
#       {'ij_f2':np.zeros(shape=2)+ij_m+ij_range},
#        ],
#      
#      [{'s_d':np.linspace(s_m-s_range,s_m+s_range,2)},
#       {'nm_d':np.linspace(nm_d-nm_range,nm_d+nm_range,2)},
#       {'m_d1':np.zeros(shape=2)+m_m-m_range},
#       {'m_d2':np.zeros(shape=2)+m_m+m_range},
#       {'chl_d1':np.zeros(shape=2)+chl_m-chl_range},
#       {'chl_d2':np.zeros(shape=2)+chl_m+chl_range},
#       {'tl_d1':np.zeros(shape=2)+tl_m-tl_range},
#       {'tl_d2':np.zeros(shape=2)+tl_m+tl_range},
#       {'vwc_d1':np.zeros(shape=2)+vwc_m-vwc_range},
#       {'vwc_d2':np.zeros(shape=2)+vwc_m+vwc_range},
#       {'ij_d1':np.zeros(shape=2)+ij_m-ij_range},
#       {'ij_d2':np.zeros(shape=2)+ij_m+ij_range},
#        ],
#        
#      [{'s_m':np.linspace(s_m-s_range,s_m+s_range,2)},
#       {'nm_m':np.linspace(nm_m-nm_range,nm_m+nm_range,2)},
#       {'m_m1':np.zeros(shape=2)+m_m-m_range},
#       {'m_m2':np.zeros(shape=2)+m_m+m_range},
#       {'chl_m1':np.zeros(shape=2)+chl_m-chl_range},
#       {'chl_m2':np.zeros(shape=2)+chl_m+chl_range},
#       {'tl_m1':np.zeros(shape=2)+tl_m-tl_range},
#       {'tl_m2':np.zeros(shape=2)+tl_m+tl_range},
#       {'vwc_m1':np.zeros(shape=2)+vwc_m-vwc_range},
#       {'vwc_m2':np.zeros(shape=2)+vwc_m+vwc_range},
#       {'ij_m1':np.zeros(shape=2)+ij_m-ij_range},
#       {'ij_m2':np.zeros(shape=2)+ij_m+ij_range},
#        ],
#       
#      [{'s_w':np.linspace(s_m-s_range,s_m+s_range,2)},
#       {'nm_w':np.linspace(nm_w-nm_range,nm_w+nm_range,2)},
#       {'m_w1':np.zeros(shape=2)+m_m-m_range},
#       {'m_w2':np.zeros(shape=2)+m_m+m_range},
#       {'chl_w1':np.zeros(shape=2)+chl_m-chl_range},
#       {'chl_w2':np.zeros(shape=2)+chl_m+chl_range},
#       {'tl_w1':np.zeros(shape=2)+tl_m-tl_range},
#       {'tl_w2':np.zeros(shape=2)+tl_m+tl_range},
#       {'vwc_w1':np.zeros(shape=2)+vwc_m-vwc_range},
#       {'vwc_w2':np.zeros(shape=2)+vwc_m+vwc_range},
#       {'ij_w1':np.zeros(shape=2)+ij_m-ij_range},
#       {'ij_w2':np.zeros(shape=2)+ij_m+ij_range},
#        ], 
#      
#      [{'s_s':np.linspace(s_m-s_range,s_m+s_range,2)},
#       {'nm_s':np.linspace(nm_s-nm_range,nm_s+nm_range,2)},
#       {'m_s1':np.zeros(shape=2)+m_m-m_range},
#       {'m_s2':np.zeros(shape=2)+m_m+m_range},
#       {'chl_s1':np.zeros(shape=2)+chl_m-chl_range},
#       {'chl_s2':np.zeros(shape=2)+chl_m+chl_range},
#       {'tl_s1':np.zeros(shape=2)+tl_m-tl_range},
#       {'tl_s2':np.zeros(shape=2)+tl_m+tl_range},
#       {'vwc_s1':np.zeros(shape=2)+vwc_m-vwc_range},
#       {'vwc_s2':np.zeros(shape=2)+vwc_m+vwc_range},
#       {'ij_s1':np.zeros(shape=2)+ij_m-ij_range},
#       {'ij_s2':np.zeros(shape=2)+ij_m+ij_range},
#        ]]

#variable nm/s
#dict=[
#      [{'s_f':np.linspace(s_f-s_range,s_f+s_range,2)},
#       {'nm_f':np.linspace(nm_f-nm_range,nm_f+nm_range,2)},
#       {'m_f1':np.zeros(shape=2)+m_m-m_range},
#       {'m_f2':np.zeros(shape=2)+m_m+m_range},
#       {'chl_f1':np.zeros(shape=2)+chl_m-chl_range},
#       {'chl_f2':np.zeros(shape=2)+chl_m+chl_range},
#       {'tl_f1':np.zeros(shape=2)+tl_m-tl_range},
#       {'tl_f2':np.zeros(shape=2)+tl_m+tl_range},
#       {'vwc_f1':np.zeros(shape=2)+vwc_m-vwc_range},
#       {'vwc_f2':np.zeros(shape=2)+vwc_m+vwc_range},
#       {'ij_f1':np.zeros(shape=2)+ij_m-ij_range},
#       {'ij_f2':np.zeros(shape=2)+ij_m+ij_range},
#        ],
#      
#      [{'s_d':np.linspace(s_d-s_range,s_d+s_range,2)},
#       {'nm_d':np.linspace(nm_d-nm_range,nm_d+nm_range,2)},
#       {'m_d1':np.zeros(shape=2)+m_m-m_range},
#       {'m_d2':np.zeros(shape=2)+m_m+m_range},
#       {'chl_d1':np.zeros(shape=2)+chl_m-chl_range},
#       {'chl_d2':np.zeros(shape=2)+chl_m+chl_range},
#       {'tl_d1':np.zeros(shape=2)+tl_m-tl_range},
#       {'tl_d2':np.zeros(shape=2)+tl_m+tl_range},
#       {'vwc_d1':np.zeros(shape=2)+vwc_m-vwc_range},
#       {'vwc_d2':np.zeros(shape=2)+vwc_m+vwc_range},
#       {'ij_d1':np.zeros(shape=2)+ij_m-ij_range},
#       {'ij_d2':np.zeros(shape=2)+ij_m+ij_range},
#        ],
#        
#      [{'s_m':np.linspace(s_m-s_range,s_m+s_range,2)},
#       {'nm_m':np.linspace(nm_m-nm_range,nm_m+nm_range,2)},
#       {'m_m1':np.zeros(shape=2)+m_m-m_range},
#       {'m_m2':np.zeros(shape=2)+m_m+m_range},
#       {'chl_m1':np.zeros(shape=2)+chl_m-chl_range},
#       {'chl_m2':np.zeros(shape=2)+chl_m+chl_range},
#       {'tl_m1':np.zeros(shape=2)+tl_m-tl_range},
#       {'tl_m2':np.zeros(shape=2)+tl_m+tl_range},
#       {'vwc_m1':np.zeros(shape=2)+vwc_m-vwc_range},
#       {'vwc_m2':np.zeros(shape=2)+vwc_m+vwc_range},
#       {'ij_m1':np.zeros(shape=2)+ij_m-ij_range},
#       {'ij_m2':np.zeros(shape=2)+ij_m+ij_range},
#        ],
#       
#      [{'s_w':np.linspace(s_w-s_range,s_w+s_range,2)},
#       {'nm_w':np.linspace(nm_w-nm_range,nm_w+nm_range,2)},
#       {'m_w1':np.zeros(shape=2)+m_m-m_range},
#       {'m_w2':np.zeros(shape=2)+m_m+m_range},
#       {'chl_w1':np.zeros(shape=2)+chl_m-chl_range},
#       {'chl_w2':np.zeros(shape=2)+chl_m+chl_range},
#       {'tl_w1':np.zeros(shape=2)+tl_m-tl_range},
#       {'tl_w2':np.zeros(shape=2)+tl_m+tl_range},
#       {'vwc_w1':np.zeros(shape=2)+vwc_m-vwc_range},
#       {'vwc_w2':np.zeros(shape=2)+vwc_m+vwc_range},
#       {'ij_w1':np.zeros(shape=2)+ij_m-ij_range},
#       {'ij_w2':np.zeros(shape=2)+ij_m+ij_range},
#        ], 
#      
#      [{'s_s':np.linspace(s_s-s_range,s_s+s_range,2)},
#       {'nm_s':np.linspace(nm_s-nm_range,nm_s+nm_range,2)},
#       {'m_s1':np.zeros(shape=2)+m_m-m_range},
#       {'m_s2':np.zeros(shape=2)+m_m+m_range},
#       {'chl_s1':np.zeros(shape=2)+chl_m-chl_range},
#       {'chl_s2':np.zeros(shape=2)+chl_m+chl_range},
#       {'tl_s1':np.zeros(shape=2)+tl_m-tl_range},
#       {'tl_s2':np.zeros(shape=2)+tl_m+tl_range},
#       {'vwc_s1':np.zeros(shape=2)+vwc_m-vwc_range},
#       {'vwc_s2':np.zeros(shape=2)+vwc_m+vwc_range},
#       {'ij_s1':np.zeros(shape=2)+ij_m-ij_range},
#       {'ij_s2':np.zeros(shape=2)+ij_m+ij_range},
#        ]]

##variable m  
#dict=[
#      [{'s_f':np.linspace(s_m-s_range,s_m+s_range,2)},
#       {'nm_f':np.linspace(nm_m-nm_range,nm_m+nm_range,2)},
#       {'m_f1':np.zeros(shape=2)+m_f-m_range},
#       {'m_f2':np.zeros(shape=2)+m_f+m_range},
#       {'chl_f1':np.zeros(shape=2)+chl_m-chl_range},
#       {'chl_f2':np.zeros(shape=2)+chl_m+chl_range},
#       {'tl_f1':np.zeros(shape=2)+tl_m-tl_range},
#       {'tl_f2':np.zeros(shape=2)+tl_m+tl_range},
#       {'vwc_f1':np.zeros(shape=2)+vwc_m-vwc_range},
#       {'vwc_f2':np.zeros(shape=2)+vwc_m+vwc_range},
#       {'ij_f1':np.zeros(shape=2)+ij_m-ij_range},
#       {'ij_f2':np.zeros(shape=2)+ij_m+ij_range},
#        ],
#      
#      [{'s_d':np.linspace(s_m-s_range,s_m+s_range,2)},
#       {'nm_d':np.linspace(nm_m-nm_range,nm_m+nm_range,2)},
#       {'m_d1':np.zeros(shape=2)+m_d-m_range},
#       {'m_d2':np.zeros(shape=2)+m_d+m_range},
#       {'chl_d1':np.zeros(shape=2)+chl_m-chl_range},
#       {'chl_d2':np.zeros(shape=2)+chl_m+chl_range},
#       {'tl_d1':np.zeros(shape=2)+tl_m-tl_range},
#       {'tl_d2':np.zeros(shape=2)+tl_m+tl_range},
#       {'vwc_d1':np.zeros(shape=2)+vwc_m-vwc_range},
#       {'vwc_d2':np.zeros(shape=2)+vwc_m+vwc_range},
#       {'ij_d1':np.zeros(shape=2)+ij_m-ij_range},
#       {'ij_d2':np.zeros(shape=2)+ij_m+ij_range},
#        ],
#        
#      [{'s_m':np.linspace(s_m-s_range,s_m+s_range,2)},
#       {'nm_m':np.linspace(nm_m-nm_range,nm_m+nm_range,2)},
#       {'m_m1':np.zeros(shape=2)+m_m-m_range},
#       {'m_m2':np.zeros(shape=2)+m_m+m_range},
#       {'chl_m1':np.zeros(shape=2)+chl_m-chl_range},
#       {'chl_m2':np.zeros(shape=2)+chl_m+chl_range},
#       {'tl_m1':np.zeros(shape=2)+tl_m-tl_range},
#       {'tl_m2':np.zeros(shape=2)+tl_m+tl_range},
#       {'vwc_m1':np.zeros(shape=2)+vwc_m-vwc_range},
#       {'vwc_m2':np.zeros(shape=2)+vwc_m+vwc_range},
#       {'ij_m1':np.zeros(shape=2)+ij_m-ij_range},
#       {'ij_m2':np.zeros(shape=2)+ij_m+ij_range},
#        ],
#       
#      [{'s_w':np.linspace(s_m-s_range,s_m+s_range,2)},
#       {'nm_w':np.linspace(nm_m-nm_range,nm_m+nm_range,2)},
#       {'m_hm1':np.zeros(shape=2)+m_w-m_range},
#       {'m_hm2':np.zeros(shape=2)+m_w+m_range},
#       {'chl_w1':np.zeros(shape=2)+chl_m-chl_range},
#       {'chl_w2':np.zeros(shape=2)+chl_m+chl_range},
#       {'tl_w1':np.zeros(shape=2)+tl_m-tl_range},
#       {'tl_w2':np.zeros(shape=2)+tl_m+tl_range},
#       {'vwc_w1':np.zeros(shape=2)+vwc_m-vwc_range},
#       {'vwc_w2':np.zeros(shape=2)+vwc_m+vwc_range},
#       {'ij_w1':np.zeros(shape=2)+ij_m-ij_range},
#       {'ij_w2':np.zeros(shape=2)+ij_m+ij_range},
#        ], 
#      
#      [{'s_s':np.linspace(s_m-s_range,s_m+s_range,2)},
#       {'nm_s':np.linspace(nm_m-nm_range,nm_m+nm_range,2)},
#       {'m_s1':np.zeros(shape=2)+m_s-m_range},
#       {'m_s2':np.zeros(shape=2)+m_s+m_range},
#       {'chl_s1':np.zeros(shape=2)+chl_m-chl_range},
#       {'chl_s2':np.zeros(shape=2)+chl_m+chl_range},
#       {'tl_s1':np.zeros(shape=2)+tl_m-tl_range},
#       {'tl_s2':np.zeros(shape=2)+tl_m+tl_range},
#       {'vwc_s1':np.zeros(shape=2)+vwc_m-vwc_range},
#       {'vwc_s2':np.zeros(shape=2)+vwc_m+vwc_range},
#       {'ij_s1':np.zeros(shape=2)+ij_m-ij_range},
#       {'ij_s2':np.zeros(shape=2)+ij_m+ij_range},
#        ]]
 
#variable chl        
#dict=[
#      [{'s_f':np.linspace(s_m-s_range,s_m+s_range,2)},
#       {'nm_f':np.linspace(nm_m-nm_range,nm_m+nm_range,2)},
#       {'m_f1':np.zeros(shape=2)+m_m-m_range},
#       {'m_f2':np.zeros(shape=2)+m_m+m_range},
#       {'chl_f1':np.zeros(shape=2)+chl_f-chl_range},
#       {'chl_f2':np.zeros(shape=2)+chl_f+chl_range},
#       {'tl_f1':np.zeros(shape=2)+tl_m-tl_range},
#       {'tl_f2':np.zeros(shape=2)+tl_m+tl_range},
#       {'vwc_f1':np.zeros(shape=2)+vwc_m-vwc_range},
#       {'vwc_f2':np.zeros(shape=2)+vwc_m+vwc_range},
#       {'ij_f1':np.zeros(shape=2)+ij_m-ij_range},
#       {'ij_f2':np.zeros(shape=2)+ij_m+ij_range},
#        ],
#      
#      [{'s_d':np.linspace(s_m-s_range,s_m+s_range,2)},
#       {'nm_d':np.linspace(nm_m-nm_range,nm_m+nm_range,2)},
#       {'m_d1':np.zeros(shape=2)+m_m-m_range},
#       {'m_d2':np.zeros(shape=2)+m_m+m_range},
#       {'chl_d1':np.zeros(shape=2)+chl_d-chl_range},
#       {'chl_d2':np.zeros(shape=2)+chl_d+chl_range},
#       {'tl_d1':np.zeros(shape=2)+tl_m-tl_range},
#       {'tl_d2':np.zeros(shape=2)+tl_m+tl_range},
#       {'vwc_d1':np.zeros(shape=2)+vwc_m-vwc_range},
#       {'vwc_d2':np.zeros(shape=2)+vwc_m+vwc_range},
#       {'ij_d1':np.zeros(shape=2)+ij_m-ij_range},
#       {'ij_d2':np.zeros(shape=2)+ij_m+ij_range},
#        ],
#        
#      [{'s_m':np.linspace(s_m-s_range,s_m+s_range,2)},
#       {'nm_m':np.linspace(nm_m-nm_range,nm_m+nm_range,2)},
#       {'m_m1':np.zeros(shape=2)+m_m-m_range},
#       {'m_m2':np.zeros(shape=2)+m_m+m_range},
#       {'chl_m1':np.zeros(shape=2)+chl_m-chl_range},
#       {'chl_m2':np.zeros(shape=2)+chl_m+chl_range},
#       {'tl_m1':np.zeros(shape=2)+tl_m-tl_range},
#       {'tl_m2':np.zeros(shape=2)+tl_m+tl_range},
#       {'vwc_m1':np.zeros(shape=2)+vwc_m-vwc_range},
#       {'vwc_m2':np.zeros(shape=2)+vwc_m+vwc_range},
#       {'ij_m1':np.zeros(shape=2)+ij_m-ij_range},
#       {'ij_m2':np.zeros(shape=2)+ij_m+ij_range},
#        ],
#       
#      [{'s_w':np.linspace(s_m-s_range,s_m+s_range,2)},
#       {'nm_w':np.linspace(nm_m-nm_range,nm_m+nm_range,2)},
#       {'m_w1':np.zeros(shape=2)+m_m-m_range},
#       {'m_w2':np.zeros(shape=2)+m_m+m_range},
#       {'chl_w1':np.zeros(shape=2)+chl_w-chl_range},
#       {'chl_w2':np.zeros(shape=2)+chl_w+chl_range},
#       {'tl_w1':np.zeros(shape=2)+tl_m-tl_range},
#       {'tl_w2':np.zeros(shape=2)+tl_m+tl_range},
#       {'vwc_w1':np.zeros(shape=2)+vwc_m-vwc_range},
#       {'vwc_w2':np.zeros(shape=2)+vwc_m+vwc_range},
#       {'ij_w1':np.zeros(shape=2)+ij_m-ij_range},
#       {'ij_w2':np.zeros(shape=2)+ij_m+ij_range},
#        ], 
#      
#      [{'s_s':np.linspace(s_m-s_range,s_m+s_range,2)},
#       {'nm_s':np.linspace(nm_m-nm_range,nm_m+nm_range,2)},
#       {'m_s1':np.zeros(shape=2)+m_m-m_range},
#       {'m_s2':np.zeros(shape=2)+m_m+m_range},
#       {'chl_s1':np.zeros(shape=2)+chl_s-chl_range},
#       {'chl_s2':np.zeros(shape=2)+chl_s+chl_range},
#       {'tl_s1':np.zeros(shape=2)+tl_m-tl_range},
#       {'tl_s2':np.zeros(shape=2)+tl_m+tl_range},
#       {'vwc_s1':np.zeros(shape=2)+vwc_m-vwc_range},
#       {'vwc_s2':np.zeros(shape=2)+vwc_m+vwc_range},
#       {'ij_s1':np.zeros(shape=2)+ij_m-ij_range},
#       {'ij_s2':np.zeros(shape=2)+ij_m+ij_range},
#        ]]

#variable tl
#dict=[
#      [{'s_f':np.linspace(s_m-s_range,s_m+s_range,2)},
#       {'nm_f':np.linspace(nm_m-nm_range,nm_m+nm_range,2)},
#       {'m_f1':np.zeros(shape=2)+m_m-m_range},
#       {'m_f2':np.zeros(shape=2)+m_m+m_range},
#       {'chl_f1':np.zeros(shape=2)+chl_m-chl_range},
#       {'chl_f2':np.zeros(shape=2)+chl_m+chl_range},
#       {'tl_f1':np.zeros(shape=2)+tl_f-tl_range},
#       {'tl_f2':np.zeros(shape=2)+tl_f+tl_range},
#       {'vwc_f1':np.zeros(shape=2)+vwc_m-vwc_range},
#       {'vwc_f2':np.zeros(shape=2)+vwc_m+vwc_range},
#       {'ij_f1':np.zeros(shape=2)+ij_m-ij_range},
#       {'ij_f2':np.zeros(shape=2)+ij_m+ij_range},
#        ],
#      
#      [{'s_d':np.linspace(s_m-s_range,s_m+s_range,2)},
#       {'nm_d':np.linspace(nm_m-nm_range,nm_m+nm_range,2)},
#       {'m_d1':np.zeros(shape=2)+m_m-m_range},
#       {'m_d2':np.zeros(shape=2)+m_m+m_range},
#       {'chl_d1':np.zeros(shape=2)+chl_m-chl_range},
#       {'chl_d2':np.zeros(shape=2)+chl_m+chl_range},
#       {'tl_d1':np.zeros(shape=2)+tl_d-tl_range},
#       {'tl_d2':np.zeros(shape=2)+tl_d+tl_range},
#       {'vwc_d1':np.zeros(shape=2)+vwc_m-vwc_range},
#       {'vwc_d2':np.zeros(shape=2)+vwc_m+vwc_range},
#       {'ij_d1':np.zeros(shape=2)+ij_m-ij_range},
#       {'ij_d2':np.zeros(shape=2)+ij_m+ij_range},
#        ],
#        
#      [{'s_m':np.linspace(s_m-s_range,s_m+s_range,2)},
#       {'nm_m':np.linspace(nm_m-nm_range,nm_m+nm_range,2)},
#       {'m_m1':np.zeros(shape=2)+m_m-m_range},
#       {'m_m2':np.zeros(shape=2)+m_m+m_range},
#       {'chl_m1':np.zeros(shape=2)+chl_m-chl_range},
#       {'chl_m2':np.zeros(shape=2)+chl_m+chl_range},
#       {'tl_m1':np.zeros(shape=2)+tl_m-tl_range},
#       {'tl_m2':np.zeros(shape=2)+tl_m+tl_range},
#       {'vwc_m1':np.zeros(shape=2)+vwc_m-vwc_range},
#       {'vwc_m2':np.zeros(shape=2)+vwc_m+vwc_range},
#       {'ij_m1':np.zeros(shape=2)+ij_m-ij_range},
#       {'ij_m2':np.zeros(shape=2)+ij_m+ij_range},
#        ],
#       
#      [{'s_w':np.linspace(s_m-s_range,s_m+s_range,2)},
#       {'nm_w':np.linspace(nm_m-nm_range,nm_m+nm_range,2)},
#       {'m_w1':np.zeros(shape=2)+m_m-m_range},
#       {'m_w2':np.zeros(shape=2)+m_m+m_range},
#       {'chl_w1':np.zeros(shape=2)+chl_m-chl_range},
#       {'chl_w2':np.zeros(shape=2)+chl_m+chl_range},
#       {'tl_w1':np.zeros(shape=2)+tl_w-tl_range},
#       {'tl_w2':np.zeros(shape=2)+tl_w+tl_range},
#       {'vwc_w1':np.zeros(shape=2)+vwc_m-vwc_range},
#       {'vwc_w2':np.zeros(shape=2)+vwc_m+vwc_range},
#       {'ij_w1':np.zeros(shape=2)+ij_m-ij_range},
#       {'ij_w2':np.zeros(shape=2)+ij_m+ij_range},
#        ], 
#      
#      [{'s_s':np.linspace(s_m-s_range,s_m+s_range,2)},
#       {'nm_s':np.linspace(nm_m-nm_range,nm_m+nm_range,2)},
#       {'m_s1':np.zeros(shape=2)+m_m-m_range},
#       {'m_s2':np.zeros(shape=2)+m_m+m_range},
#       {'chl_s1':np.zeros(shape=2)+chl_m-chl_range},
#       {'chl_s2':np.zeros(shape=2)+chl_m+chl_range},
#       {'tl_s1':np.zeros(shape=2)+tl_s-tl_range},
#       {'tl_s2':np.zeros(shape=2)+tl_s+tl_range},
#       {'vwc_s1':np.zeros(shape=2)+vwc_m-vwc_range},
#       {'vwc_s2':np.zeros(shape=2)+vwc_m+vwc_range},
#       {'ij_s1':np.zeros(shape=2)+ij_m-ij_range},
#       {'ij_s2':np.zeros(shape=2)+ij_m+ij_range},
#        ]]

        
# 
##variable vwc
#dict=[
#      [{'s_f':np.linspace(s_m-s_range,s_m+s_range,2)},
#       {'nm_f':np.linspace(nm_m-nm_range,nm_m+nm_range,2)},
#       {'m_f1':np.zeros(shape=2)+m_m-m_range},
#       {'m_f2':np.zeros(shape=2)+m_m+m_range},
#       {'chl_f1':np.zeros(shape=2)+chl_m-chl_range},
#       {'chl_f2':np.zeros(shape=2)+chl_m+chl_range},
#       {'tl_f1':np.zeros(shape=2)+tl_m-tl_range},
#       {'tl_f2':np.zeros(shape=2)+tl_m+tl_range},
#       {'vwc_f1':np.zeros(shape=2)+vwc_f-vwc_range},
#       {'vwc_f2':np.zeros(shape=2)+vwc_f+vwc_range},
#       {'ij_f1':np.zeros(shape=2)+ij_m-ij_range},
#       {'ij_f2':np.zeros(shape=2)+ij_m+ij_range},
#        ],
#      
#      [{'s_d':np.linspace(s_m-s_range,s_m+s_range,2)},
#       {'nm_d':np.linspace(nm_m-nm_range,nm_m+nm_range,2)},
#       {'m_d1':np.zeros(shape=2)+m_m-m_range},
#       {'m_d2':np.zeros(shape=2)+m_m+m_range},
#       {'chl_d1':np.zeros(shape=2)+chl_m-chl_range},
#       {'chl_d2':np.zeros(shape=2)+chl_m+chl_range},
#       {'tl_d1':np.zeros(shape=2)+tl_m-tl_range},
#       {'tl_d2':np.zeros(shape=2)+tl_m+tl_range},
#       {'vwc_d1':np.zeros(shape=2)+vwc_d-vwc_range},
#       {'vwc_d2':np.zeros(shape=2)+vwc_d+vwc_range},
#       {'ij_d1':np.zeros(shape=2)+ij_m-ij_range},
#       {'ij_d2':np.zeros(shape=2)+ij_m+ij_range},
#        ],
#        
#      [{'s_m':np.linspace(s_m-s_range,s_m+s_range,2)},
#       {'nm_m':np.linspace(nm_m-nm_range,nm_m+nm_range,2)},
#       {'m_m1':np.zeros(shape=2)+m_m-m_range},
#       {'m_m2':np.zeros(shape=2)+m_m+m_range},
#       {'chl_m1':np.zeros(shape=2)+chl_m-chl_range},
#       {'chl_m2':np.zeros(shape=2)+chl_m+chl_range},
#       {'tl_m1':np.zeros(shape=2)+tl_m-tl_range},
#       {'tl_m2':np.zeros(shape=2)+tl_m+tl_range},
#       {'vwc_m1':np.zeros(shape=2)+vwc_m-vwc_range},
#       {'vwc_m2':np.zeros(shape=2)+vwc_m+vwc_range},
#       {'ij_m1':np.zeros(shape=2)+ij_m-ij_range},
#       {'ij_m2':np.zeros(shape=2)+ij_m+ij_range},
#        ],
#       
#      [{'s_w':np.linspace(s_m-s_range,s_m+s_range,2)},
#       {'nm_w':np.linspace(nm_m-nm_range,nm_m+nm_range,2)},
#       {'m_w1':np.zeros(shape=2)+m_m-m_range},
#       {'m_w2':np.zeros(shape=2)+m_m+m_range},
#       {'chl_w1':np.zeros(shape=2)+chl_m-chl_range},
#       {'chl_w2':np.zeros(shape=2)+chl_m+chl_range},
#       {'tl_w1':np.zeros(shape=2)+tl_m-tl_range},
#       {'tl_w2':np.zeros(shape=2)+tl_m+tl_range},
#       {'vwc_w1':np.zeros(shape=2)+vwc_w-vwc_range},
#       {'vwc_w2':np.zeros(shape=2)+vwc_w+vwc_range},
#       {'ij_w1':np.zeros(shape=2)+ij_m-ij_range},
#       {'ij_w2':np.zeros(shape=2)+ij_m+ij_range},
#        ], 
#      
#      [{'s_s':np.linspace(s_m-s_range,s_m+s_range,2)},
#       {'nm_s':np.linspace(nm_m-nm_range,nm_m+nm_range,2)},
#       {'m_s1':np.zeros(shape=2)+m_m-m_range},
#       {'m_s2':np.zeros(shape=2)+m_m+m_range},
#       {'chl_s1':np.zeros(shape=2)+chl_m-chl_range},
#       {'chl_s2':np.zeros(shape=2)+chl_m+chl_range},
#       {'tl_s1':np.zeros(shape=2)+tl_m-tl_range},
#       {'tl_s2':np.zeros(shape=2)+tl_m+tl_range},
#       {'vwc_s1':np.zeros(shape=2)+vwc_s-vwc_range},
#       {'vwc_s2':np.zeros(shape=2)+vwc_s+vwc_range},
#       {'ij_s1':np.zeros(shape=2)+ij_m-ij_range},
#       {'ij_s2':np.zeros(shape=2)+ij_m+ij_range},
#        ]]       

 
#variable ij
#dict=[
#      [{'s_f':np.linspace(s_m-0.1*s_m,s_m+0.1*s_m,2)},  ##STOPPED HERE
#       {'nm_f':np.linspace(nm_m-nm_range,nm_m+nm_range,2)},
#       {'m_f1':np.zeros(shape=2)+m_m-m_range},
#       {'m_f2':np.zeros(shape=2)+m_m+m_range},
#       {'chl_f1':np.zeros(shape=2)+chl_m-chl_range},
#       {'chl_f2':np.zeros(shape=2)+chl_m+chl_range},
#       {'tl_f1':np.zeros(shape=2)+tl_m-tl_range},
#       {'tl_f2':np.zeros(shape=2)+tl_m+tl_range},
#       {'vwc_f1':np.zeros(shape=2)+vwc_m-vwc_range},
#       {'vwc_f2':np.zeros(shape=2)+vwc_m+vwc_range},
#       {'ij_f1':np.zeros(shape=2)+ij_f-ij_range},
#       {'ij_f2':np.zeros(shape=2)+ij_f+ij_range},
#        ],
#      
#      [{'s_d':np.linspace(s_m-s_range,s_m+s_range,2)},
#       {'nm_d':np.linspace(nm_m-nm_range,nm_m+nm_range,2)},
#       {'m_d1':np.zeros(shape=2)+m_m-m_range},
#       {'m_d2':np.zeros(shape=2)+m_m+m_range},
#       {'chl_d1':np.zeros(shape=2)+chl_m-chl_range},
#       {'chl_d2':np.zeros(shape=2)+chl_m+chl_range},
#       {'tl_d1':np.zeros(shape=2)+tl_m-tl_range},
#       {'tl_d2':np.zeros(shape=2)+tl_m+tl_range},
#       {'vwc_d1':np.zeros(shape=2)+vwc_m-vwc_range},
#       {'vwc_d2':np.zeros(shape=2)+vwc_m+vwc_range},
#       {'ij_d1':np.zeros(shape=2)+ij_d-ij_range},
#       {'ij_d2':np.zeros(shape=2)+ij_d+ij_range},
#        ],
#        
#      [{'s_m':np.linspace(s_m-s_range,s_m+s_range,2)},
#       {'nm_m':np.linspace(nm_m-nm_range,nm_m+nm_range,2)},
#       {'m_m1':np.zeros(shape=2)+m_m-m_range},
#       {'m_m2':np.zeros(shape=2)+m_m+m_range},
#       {'chl_m1':np.zeros(shape=2)+chl_m-chl_range},
#       {'chl_m2':np.zeros(shape=2)+chl_m+chl_range},
#       {'tl_m1':np.zeros(shape=2)+tl_m-tl_range},
#       {'tl_m2':np.zeros(shape=2)+tl_m+tl_range},
#       {'vwc_m1':np.zeros(shape=2)+vwc_m-vwc_range},
#       {'vwc_m2':np.zeros(shape=2)+vwc_m+vwc_range},
#       {'ij_m1':np.zeros(shape=2)+ij_m-ij_range},
#       {'ij_m2':np.zeros(shape=2)+ij_m+ij_range},
#        ],
#       
#      [{'s_w':np.linspace(s_m-s_range,s_m+s_range,2)},
#       {'nm_w':np.linspace(nm_m-nm_range,nm_m+nm_range,2)},
#       {'m_w1':np.zeros(shape=2)+m_m-m_range},
#       {'m_w2':np.zeros(shape=2)+m_m+m_range},
#       {'chl_w1':np.zeros(shape=2)+chl_m-chl_range},
#       {'chl_w2':np.zeros(shape=2)+chl_m+chl_range},
#       {'tl_w1':np.zeros(shape=2)+tl_m-tl_range},
#       {'tl_w2':np.zeros(shape=2)+tl_m+tl_range},
#       {'vwc_w1':np.zeros(shape=2)+vwc_m-vwc_range},
#       {'vwc_w2':np.zeros(shape=2)+vwc_m+vwc_range},
#       {'ij_w1':np.zeros(shape=2)+ij_w-ij_range},
#       {'ij_w2':np.zeros(shape=2)+ij_w+ij_range},
#        ], 
#      
#      [{'s_s':np.linspace(s_m-s_range,s_m+s_range,2)},
#       {'nm_s':np.linspace(nm_m-nm_range,nm_m+nm_range,2)},
#       {'m_s1':np.zeros(shape=2)+m_m-m_range},
#       {'m_s2':np.zeros(shape=2)+m_m+m_range},
#       {'chl_s1':np.zeros(shape=2)+chl_m-chl_range},
#       {'chl_s2':np.zeros(shape=2)+chl_m+chl_range},
#       {'tl_s1':np.zeros(shape=2)+tl_m-tl_range},
#       {'tl_s2':np.zeros(shape=2)+tl_m+tl_range},
#       {'vwc_s1':np.zeros(shape=2)+vwc_m-vwc_range},
#       {'vwc_s2':np.zeros(shape=2)+vwc_m+vwc_range},
#       {'ij_s1':np.zeros(shape=2)+ij_s-ij_range},
#       {'ij_s2':np.zeros(shape=2)+ij_s+ij_range},
#        ]]       


#traits vary for each community

dict=[
      [{'s_f':np.linspace(s_f-s_range,s_f+s_range,2)},
       {'nm_f':np.linspace(nm_f-nm_range,nm_f+nm_range,2)},
       {'m_f1':np.zeros(shape=2)+m_f-m_range},
       {'m_f2':np.zeros(shape=2)+m_f+m_range},
       {'chl_f1':np.zeros(shape=2)+chl_f-chl_range},
       {'chl_f2':np.zeros(shape=2)+chl_f+chl_range},
       {'tl_f1':np.zeros(shape=2)+tl_f-tl_range},
       {'tl_f2':np.zeros(shape=2)+tl_f+tl_range},
       {'vwc_f1':np.zeros(shape=2)+vwc_m-vwc_range},
       {'vwc_f2':np.zeros(shape=2)+vwc_m+vwc_range},
       {'ij_f1':np.zeros(shape=2)+ij_f-ij_range},
       {'ij_f2':np.zeros(shape=2)+ij_f+ij_range},
        ],
      
      [{'s_d':np.linspace(s_d-s_range,s_d+s_range,2)},
       {'nm_d':np.linspace(nm_d-nm_range,nm_d+nm_range,2)},
       {'m_d1':np.zeros(shape=2)+m_d-m_range},
       {'m_d2':np.zeros(shape=2)+m_d+m_range},
       {'chl_d1':np.zeros(shape=2)+chl_d-chl_range},
       {'chl_d2':np.zeros(shape=2)+chl_d+chl_range},
       {'tl_d1':np.zeros(shape=2)+tl_d-tl_range},
       {'tl_d2':np.zeros(shape=2)+tl_d+tl_range},
       {'vwc_d1':np.zeros(shape=2)+vwc_m-vwc_range},
       {'vwc_d2':np.zeros(shape=2)+vwc_m+vwc_range},
        {'ij_d1':np.zeros(shape=2)+ij_d-ij_range},
        {'ij_d2':np.zeros(shape=2)+ij_d+ij_range},
        ],
        
      [{'s_m':np.linspace(s_m-s_range,s_m+s_range,2)},
       {'nm_m':np.linspace(nm_m-nm_range,nm_m+nm_range,2)},
       {'m_m1':np.zeros(shape=2)+m_m-m_range},
       {'m_m2':np.zeros(shape=2)+m_m+m_range},
       {'chl_m1':np.zeros(shape=2)+chl_m-chl_range},
       {'chl_m2':np.zeros(shape=2)+chl_m+chl_range},
       {'tl_m1':np.zeros(shape=2)+tl_m-tl_range},
       {'tl_m2':np.zeros(shape=2)+tl_m+tl_range},
       {'vwc_m1':np.zeros(shape=2)+vwc_m-vwc_range},
       {'vwc_m2':np.zeros(shape=2)+vwc_m+vwc_range},
       {'ij_m1':np.zeros(shape=2)+ij_m-ij_range},
       {'ij_m2':np.zeros(shape=2)+ij_m+ij_range},
        ],
       
      [{'s_w':np.linspace(s_w-s_range,s_w+s_range,2)},
       {'nm_w':np.linspace(nm_w-nm_range,nm_w+nm_range,2)},
       {'m_w1':np.zeros(shape=2)+m_w-m_range},
       {'m_w2':np.zeros(shape=2)+m_w+m_range},
       {'chl_w1':np.zeros(shape=2)+chl_w-chl_range},
       {'chl_w2':np.zeros(shape=2)+chl_w+chl_range},
       {'tl_w1':np.zeros(shape=2)+tl_w-tl_range},
       {'tl_w2':np.zeros(shape=2)+tl_w+tl_range},
       {'vwc_w1':np.zeros(shape=2)+vwc_m-vwc_range},
       {'vwc_w2':np.zeros(shape=2)+vwc_m+vwc_range},
       {'ij_w1':np.zeros(shape=2)+ij_w-ij_range},
       {'ij_w2':np.zeros(shape=2)+ij_w+ij_range},
        ], 
      
      [{'s_s':np.linspace(s_s-s_range,s_s+s_range,2)},
       {'nm_s':np.linspace(nm_s-nm_range,nm_s+nm_range,2)},
       {'m_s1':np.zeros(shape=2)+m_s-m_range},
       {'m_s2':np.zeros(shape=2)+m_s+m_range},
       {'chl_s1':np.zeros(shape=2)+chl_s-chl_range},
       {'chl_s2':np.zeros(shape=2)+chl_s+chl_range},
       {'tl_s1':np.zeros(shape=2)+tl_s-tl_range},
       {'tl_s2':np.zeros(shape=2)+tl_s+tl_range},
       {'vwc_s1':np.zeros(shape=2)+vwc_m-vwc_range},
       {'vwc_s2':np.zeros(shape=2)+vwc_m+vwc_range},
       {'ij_s1':np.zeros(shape=2)+ij_s-ij_range},
       {'ij_s2':np.zeros(shape=2)+ij_s+ij_range},
        ]] 
#---------------Define Number of Variable Parameters---------------#  

num_params=7

leaf_params=[]

#---------------Run For Loop of Different Dictionary Values---------------# 

for i in range(len(dict)):

#---------------Make Dictionaries that Contain All Combinations of Parameter Inputs---------------#  

    dict_combs_raw=[]
    for x in it.combinations(dict[i],num_params):
        dict_combs_raw+=[x]

    #merge dictionaries 
    dict_combs=[]
    for x in range(len(dict_combs_raw)):
        dict_combs+=[merge_dicts(dict_combs_raw[x])]
 

#---------------Get Rid of Repeated Parameters (i.e. s_in & s) in Parameter Combination Dictionaries---------------#  
                 
    dict_combs_new=[]

    for i in range(len(dict_combs)):
        keys=dict_combs[i].keys()
        key_param=[]
        for ii in range(len(keys)):
            key_param+=[keys[ii][0]]
        if len(np.unique(key_param))<num_params:
            continue
        else:
            dict_combs_new+=[dict_combs[i]]

#---------------Rename All Parameters in Parameter Combination Dictionaries to Model Parameter Names---------------#  

    params_int=[]


    for i in range(len(dict_combs_new)):
        params=dict_combs_new[i].copy()
        
        #rename sla
        if 's_f' in params.keys():
            params['s']=params.pop('s_f')        
        if 's_d' in params.keys():
            params['s']=params.pop('s_d')
        if 's_m' in params.keys():
            params['s']=params.pop('s_m') 
        if 's_w' in params.keys():
            params['s']=params.pop('s_w') 
        if 's_s' in params.keys():
            params['s']=params.pop('s_s')               
        
        #rename leaf n
        if 'nm_f' in params.keys():
            params['nm']=params.pop('nm_f')        
        if 'nm_d' in params.keys():
            params['nm']=params.pop('nm_d')
        if 'nm_m' in params.keys():
            params['nm']=params.pop('nm_m')
        if 'nm_w' in params.keys():
            params['nm']=params.pop('nm_w')
        if 'nm_s' in params.keys():
            params['nm']=params.pop('nm_s')            
        

            
        #rename ball berry slope parameter
        if 'm_f1' in params.keys():
            params['m']=params.pop('m_f1')
        if 'm_d1' in params.keys():
            params['m']=params.pop('m_d1')
        if 'm_m1' in params.keys():
            params['m']=params.pop('m_m1')
        if 'm_w1' in params.keys():
            params['m']=params.pop('m_w1')            
        if 'm_s1' in params.keys():
            params['m']=params.pop('m_s1')  

        if 'm_f2' in params.keys():
            params['m']=params.pop('m_f2')
        if 'm_d2' in params.keys():
            params['m']=params.pop('m_d2')
        if 'm_m2' in params.keys():
            params['m']=params.pop('m_m2')
        if 'm_w2' in params.keys():
            params['m']=params.pop('m_w2')            
        if 'm_s2' in params.keys():
            params['m']=params.pop('m_s2')              

        
        #rename chlorophyll content
        if 'chl_f1' in params.keys():
            params['chl']=params.pop('chl_f1')        
        if 'chl_d1' in params.keys():
            params['chl']=params.pop('chl_d1')
        if 'chl_m1' in params.keys():
            params['chl']=params.pop('chl_m1')
        if 'chl_w1' in params.keys():
            params['chl']=params.pop('chl_w1')
        if 'chl_s1' in params.keys():
            params['chl']=params.pop('chl_s1')            

        if 'chl_f2' in params.keys():
            params['chl']=params.pop('chl_f2')            
        if 'chl_d2' in params.keys():
            params['chl']=params.pop('chl_d2')
        if 'chl_m2' in params.keys():
            params['chl']=params.pop('chl_m2')
        if 'chl_w2' in params.keys():
            params['chl']=params.pop('chl_w2')
        if 'chl_s2' in params.keys():
            params['chl']=params.pop('chl_s2')              
 
        
        #rename temp of leaf
        if 'tl_f1' in params.keys():
            params['tl']=params.pop('tl_f1')        
        if 'tl_d1' in params.keys():
            params['tl']=params.pop('tl_d1')
        if 'tl_m1' in params.keys():
            params['tl']=params.pop('tl_m1')
        if 'tl_w1' in params.keys():
            params['tl']=params.pop('tl_w1')
        if 'tl_s1' in params.keys():
            params['tl']=params.pop('tl_s1')            

        if 'tl_f2' in params.keys():
            params['tl']=params.pop('tl_f2')            
        if 'tl_d2' in params.keys():
            params['tl']=params.pop('tl_d2')
        if 'tl_m2' in params.keys():
            params['tl']=params.pop('tl_m2')
        if 'tl_w2' in params.keys():
            params['tl']=params.pop('tl_w2')
        if 'tl_s2' in params.keys():
            params['tl']=params.pop('tl_s2')              

          
        #rename volumetric water content 
        if 'vwc_f1' in params.keys():
            params['vwc']=params.pop('vwc_f1')        
        if 'vwc_d1' in params.keys():
            params['vwc']=params.pop('vwc_d1')
        if 'vwc_m1' in params.keys():
            params['vwc']=params.pop('vwc_m1')
        if 'vwc_w1' in params.keys():
            params['vwc']=params.pop('vwc_w1')
        if 'vwc_s1' in params.keys():
            params['vwc']=params.pop('vwc_s1')            

        if 'vwc_f2' in params.keys():
            params['vwc']=params.pop('vwc_f2')            
        if 'vwc_d2' in params.keys():
            params['vwc']=params.pop('vwc_d2')
        if 'vwc_m2' in params.keys():
            params['vwc']=params.pop('vwc_m2')
        if 'vwc_w2' in params.keys():
            params['vwc']=params.pop('vwc_w2')
        if 'vwc_s2' in params.keys():
            params['vwc']=params.pop('vwc_s2')              

                
        
        #rename leaf area/angle index
        if 'ij_f1' in params.keys():
            params['ij']=params.pop('ij_f1')            
        if 'ij_d1' in params.keys():
            params['ij']=params.pop('ij_d1')
        if 'ij_m1' in params.keys():
            params['ij']=params.pop('ij_m1')
        if 'ij_w1' in params.keys():
            params['ij']=params.pop('ij_w1')
        if 'ij_s1' in params.keys():
            params['ij']=params.pop('ij_s1')    
                
        if 'ij_f2' in params.keys():
            params['ij']=params.pop('ij_f2')            
        if 'ij_d2' in params.keys():
            params['ij']=params.pop('ij_d2')
        if 'ij_m2' in params.keys():
            params['ij']=params.pop('ij_m2')
        if 'ij_w2' in params.keys():
            params['ij']=params.pop('ij_w2')
        if 'ij_s2' in params.keys():
            params['ij']=params.pop('ij_s2') 
                  
        
        params_int+=[params]       
    
    leaf_params+=params_int