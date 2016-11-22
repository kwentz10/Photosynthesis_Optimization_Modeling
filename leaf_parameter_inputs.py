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
nm_f=((s_f*(100.0**2))*0.077+20.25)/1000.0 #leaf N (gN/gC)
flnr_f=0.55 #fraction of leaf nitrogen in rubisco (g leaf N/g rubisco)
chl_f=200 #chlorophyll content (umol chl/m2)
tl_f=35.0+273.15 #temp of leaf (K)-->f(leaf area,angle,height)
vwc_f=0.10 #soil volumetric water content (cm3/cm3)
ij_f=0.92 #index of light penetration (low leaf area and high leaf angle) (unitless)-->f(leaf area & angle)

#dry meadow
s_d=0.017 #specific leaf area (m2/g)
nm_d=((s_d*(100.0**2))*0.077+20.25)/1000.0 #leaf N (gN/gC)
flnr_d=0.6 #fraction of leaf nitrogen in rubisco (g leaf N/g rubisco)
chl_d=300 #chlorophyll content (umol chl/m2)
tl_d=33.0+273.15 #temp of leaf (K)-->f(leaf area,angle,height)
vwc_d=0.12  #soil volumetric water content (cm3/cm3)
ij_d=0.94 #index of light penetration (low medium leaf area and high medium leaf angle) (unitless)-->f(leaf area & angle)

#moist meadow
s_m=0.019 #specific leaf area (m2/g)
nm_m=((s_m*(100.0**2))*0.077+20.25)/1000.0 #leaf N (gN/gC)
flnr_m=0.65 #fraction of leaf nitrogen in rubisco (g leaf N/g rubisco)
chl_m=400 #chlorophyll content (umol chl/m2)
tl_m=31.0+273.15 #temp of leaf (K)-->f(leaf area,angle,height)
vwc_m=0.15 #soil volumetric water content (cm3/cm3)
ij_m=0.96 #index of light penetration (medium leaf area and medium leaf angle) (unitless)-->f(leaf area & angle)

#wet meadow
s_w=0.022 #specific leaf area (m2/g)
nm_w=((s_w*(100.0**2))*0.077+20.25)/1000.0 #leaf N (gN/gC)
flnr_w=0.7 #fraction of leaf nitrogen in rubisco (g leaf N/g rubisco)
chl_w=500 #chlorophyll content (umol chl/m2)
tl_w=29.0+273.15 #temp of leaf (K)-->f(leaf area,angle,height)
vwc_w=0.27 #soil volumetric water content (cm3/cm3)
ij_w=0.98 #index of light penetration (high medium leaf area and low medium leaf angle) (unitless)-->f(leaf area & angle)

#snowbed
s_s=0.025 #specific leaf area (m2/g)
nm_s=((s_s*(100.0**2))*0.077+20.25)/1000.0 #leaf N (gN/gC)
flnr_s=0.75 #fraction of leaf nitrogen in rubisco (g leaf N/g rubisco)
chl_s=600 #chlorophyll content (umol chl/m2)
tl_s=27.0+273.15 #temp of leaf (K)-->f(leaf area,angle,height)
vwc_s=0.4 #soil volumetric water content (cm3/cm3)
ij_s=1.0 #index of light penetration (high leaf area and low leaf angle) (unitless)-->f(leaf area & angle)

#lambda (wue)
lambda_h=10000.0*(1.0/10.0**6) #high (umol CO2/umol H20)
lambda_hm=7600.0*(1.0/10.0**6) #high medium (umol CO2/umol H20)
lambda_m=5200.0*(1.0/10.0**6) #medium (umol CO2/umol H20)
lambda_lm=2800.0*(1.0/10.0**6) #low medium (umol CO2/umol H20)
lambda_l=400.0*(1.0/10.0**6) #low (umol CO2/umol H20)

#ball berry parameter (sensitivity of stomatal conductance)
m_h=25.0 #(unitless)
m_hm=20.0 #(unitless)
m_m=15.0 #(unitless)
m_lm=10.0 #(unitless)
m_l=5.0 #(unitless)

##plasticity of traits
s_range=0.0005
nm_range=0.0005
flnr_range=0.01
chl_range=75
tl_range=0.01
lambda_range=0.00001
vwc_range=0.001
m_range=0.5

##---In the Below Dictionary, sla, leaf n, and fraction of leaf n in rub are the only proportionally dependent variables---##
dict=[
      [{'s_f':np.linspace(s_f-s_range,s_f+s_range,3)},
       {'nm_f':np.linspace(nm_f-nm_range,nm_f+nm_range,3)},
#       {'flnr_f':np.linspace(flnr_f-flnr_range,flnr_f+flnr_range,3)},
#       {'lambda_l1':np.zeros(shape=3)+lambda_m-lambda_range},
#       {'lambda_l3':np.zeros(shape=3)+lambda_m+lambda_range},
       {'m_l1':np.zeros(shape=3)+m_h-m_range},
       {'m_l3':np.zeros(shape=3)+m_h+m_range},
       {'chl_f1':np.zeros(shape=3)+chl_f-chl_range},
       {'chl_f3':np.zeros(shape=3)+chl_f+chl_range},
       {'tl_f1':np.zeros(shape=3)+tl_f-tl_range},
       {'tl_f3':np.zeros(shape=3)+tl_f+tl_range},
       {'vwc_f1':np.zeros(shape=3)+vwc_f-vwc_range},
       {'vwc_f3':np.zeros(shape=3)+vwc_f+vwc_range},
        {'ij_f':np.zeros(shape=3)+ij_f},
        ],
      
      [{'s_d':np.linspace(s_d-s_range,s_d+s_range,3)},
       {'nm_d':np.linspace(nm_d-nm_range,nm_d+nm_range,3)},
#       {'flnr_d':np.linspace(flnr_d-flnr_range,flnr_d+flnr_range,3)},
#       {'lambda_lm1':np.zeros(shape=3)+lambda_m-lambda_range},
#       {'lambda_lm3':np.zeros(shape=3)+lambda_m+lambda_range},
       {'m_lm1':np.zeros(shape=3)+m_hm-m_range},
       {'m_lm3':np.zeros(shape=3)+m_hm+m_range},
       {'chl_d1':np.zeros(shape=3)+chl_d-chl_range},
       {'chl_d3':np.zeros(shape=3)+chl_d+chl_range},
       {'tl_d1':np.zeros(shape=3)+tl_d-tl_range},
       {'tl_d3':np.zeros(shape=3)+tl_d+tl_range},
       {'vwc_d1':np.zeros(shape=3)+vwc_d-vwc_range},
       {'vwc_d3':np.zeros(shape=3)+vwc_d+vwc_range},
        {'ij_d':np.zeros(shape=3)+ij_d},
        ],
        
      [{'s_m':np.linspace(s_m-s_range,s_m+s_range,3)},
       {'nm_m':np.linspace(nm_m-nm_range,nm_m+nm_range,3)},
#       {'flnr_m':np.linspace(flnr_m-flnr_range,flnr_m+flnr_range,3)},
#       {'lambda_m1':np.zeros(shape=3)+lambda_m-lambda_range},
#       {'lambda_m3':np.zeros(shape=3)+lambda_m+lambda_range},
       {'m_m1':np.zeros(shape=3)+m_m-m_range},
       {'m_m3':np.zeros(shape=3)+m_m+m_range},
       {'chl_m1':np.zeros(shape=3)+chl_m-chl_range},
       {'chl_m3':np.zeros(shape=3)+chl_m+chl_range},
       {'tl_m1':np.zeros(shape=3)+tl_m-tl_range},
       {'tl_m3':np.zeros(shape=3)+tl_m+tl_range},
       {'vwc_m1':np.zeros(shape=3)+vwc_m-vwc_range},
       {'vwc_m3':np.zeros(shape=3)+vwc_m+vwc_range},
        {'ij_m':np.zeros(shape=3)+ij_m},
        ],
       
      [{'s_w':np.linspace(s_w-s_range,s_w+s_range,3)},
       {'nm_w':np.linspace(nm_w-nm_range,nm_w+nm_range,3)},
#       {'flnr_w':np.linspace(flnr_w-flnr_range,flnr_w+flnr_range,3)},
#       {'lambda_hm1':np.zeros(shape=3)+lambda_m-lambda_range},
#       {'lambda_hm3':np.zeros(shape=3)+lambda_m+lambda_range},
       {'m_hm1':np.zeros(shape=3)+m_lm-m_range},
       {'m_hm3':np.zeros(shape=3)+m_lm+m_range},
       {'chl_w1':np.zeros(shape=3)+chl_w-chl_range},
       {'chl_w3':np.zeros(shape=3)+chl_w+chl_range},
       {'tl_w1':np.zeros(shape=3)+tl_w-tl_range},
       {'tl_w3':np.zeros(shape=3)+tl_w+tl_range},
       {'vwc_w1':np.zeros(shape=3)+vwc_w-vwc_range},
       {'vwc_w3':np.zeros(shape=3)+vwc_w+vwc_range},
        {'ij_w':np.zeros(shape=3)+ij_w},
        ], 
      
      [{'s_s':np.linspace(s_s-s_range,s_s+s_range,3)},
       {'nm_s':np.linspace(nm_s-nm_range,nm_s+nm_range,3)},
#       {'flnr_s':np.linspace(flnr_s-flnr_range,flnr_s+flnr_range,3)},
#       {'lambda_h1':np.zeros(shape=3)+lambda_m-lambda_range},
#       {'lambda_h3':np.zeros(shape=3)+lambda_m+lambda_range},
       {'m_h1':np.zeros(shape=3)+m_l-m_range},
       {'m_h3':np.zeros(shape=3)+m_l+m_range},
       {'chl_s1':np.zeros(shape=3)+chl_s-chl_range},
       {'chl_s3':np.zeros(shape=3)+chl_s+chl_range},
       {'tl_s1':np.zeros(shape=3)+tl_s-tl_range},
       {'tl_s3':np.zeros(shape=3)+tl_s+tl_range},
       {'vwc_s1':np.zeros(shape=3)+vwc_s-vwc_range},
       {'vwc_s3':np.zeros(shape=3)+vwc_s+vwc_range},
        {'ij_s':np.zeros(shape=3)+ij_s},
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
         
        #rename fraction of n in rubisco
        if 'flnr_f' in params.keys():
            params['flnr']=params.pop('flnr_f')        
        if 'flnr_d' in params.keys():
            params['flnr']=params.pop('flnr_d')
        if 'flnr_m' in params.keys():
            params['flnr']=params.pop('flnr_m')
        if 'flnr_w' in params.keys():
            params['flnr']=params.pop('flnr_w')
        if 'flnr_s' in params.keys():
            params['flnr']=params.pop('flnr_s')
        
        #rename wue  
        if 'lambda_l1' in params.keys():
            params['lamb']=params.pop('lambda_l1')
        if 'lambda_lm1' in params.keys():
            params['lamb']=params.pop('lambda_lm1')
        if 'lambda_m1' in params.keys():
            params['lamb']=params.pop('lambda_m1')
        if 'lambda_hm1' in params.keys():
            params['lamb']=params.pop('lambda_hm1')            
        if 'lambda_h1' in params.keys():
            params['lamb']=params.pop('lambda_h1') 

        if 'lambda_l2' in params.keys():
            params['lamb']=params.pop('lambda_l2')
        if 'lambda_lm2' in params.keys():
            params['lamb']=params.pop('lambda_lm2')
        if 'lambda_m2' in params.keys():
            params['lamb']=params.pop('lambda_m2')
        if 'lambda_hm2' in params.keys():
            params['lamb']=params.pop('lambda_hm2')            
        if 'lambda_h2' in params.keys():
            params['lamb']=params.pop('lambda_h2')              
     
        if 'lambda_l3' in params.keys():
            params['lamb']=params.pop('lambda_l3')
        if 'lambda_lm3' in params.keys():
            params['lamb']=params.pop('lambda_lm3')
        if 'lambda_m3' in params.keys():
            params['lamb']=params.pop('lambda_m3')
        if 'lambda_hm3' in params.keys():
            params['lamb']=params.pop('lambda_hm3')            
        if 'lambda_h3' in params.keys():
            params['lamb']=params.pop('lambda_h3')  
            
        #rename ball berry slope parameter
        if 'm_l1' in params.keys():
            params['m']=params.pop('m_l1')
        if 'm_lm1' in params.keys():
            params['m']=params.pop('m_lm1')
        if 'm_m1' in params.keys():
            params['m']=params.pop('m_m1')
        if 'm_hm1' in params.keys():
            params['m']=params.pop('m_hm1')            
        if 'm_h1' in params.keys():
            params['m']=params.pop('m_h1')  

        if 'm_l2' in params.keys():
            params['m']=params.pop('m_l2')
        if 'm_lm2' in params.keys():
            params['m']=params.pop('m_lm2')
        if 'm_m2' in params.keys():
            params['m']=params.pop('m_m2')
        if 'm_hm2' in params.keys():
            params['m']=params.pop('m_hm2')            
        if 'm_h2' in params.keys():
            params['m']=params.pop('m_h2')              
     
        if 'm_l3' in params.keys():
            params['m']=params.pop('m_l3')
        if 'm_lm3' in params.keys():
            params['m']=params.pop('m_lm3')
        if 'm_m3' in params.keys():
            params['m']=params.pop('m_m3')
        if 'm_hm3' in params.keys():
            params['m']=params.pop('m_hm3')            
        if 'm_h3' in params.keys():
            params['m']=params.pop('m_h3')  
        
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
        
        if 'chl_f3' in params.keys():
            params['chl']=params.pop('chl_f3')
        if 'chl_d3' in params.keys():
            params['chl']=params.pop('chl_d3')
        if 'chl_m3' in params.keys():
            params['chl']=params.pop('chl_m3')
        if 'chl_w3' in params.keys():
            params['chl']=params.pop('chl_w3')
        if 'chl_s3' in params.keys():
            params['chl']=params.pop('chl_s3')   
        
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

        if 'tl_f3' in params.keys():
            params['tl']=params.pop('tl_f3')            
        if 'tl_d3' in params.keys():
            params['tl']=params.pop('tl_d3')
        if 'tl_m3' in params.keys():
            params['tl']=params.pop('tl_m3')
        if 'tl_w3' in params.keys():
            params['tl']=params.pop('tl_w3')
        if 'tl_s3' in params.keys():
            params['tl']=params.pop('tl_s3')   
          
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

        if 'vwc_f3' in params.keys():
            params['vwc']=params.pop('vwc_f3')            
        if 'vwc_d3' in params.keys():
            params['vwc']=params.pop('vwc_d3')
        if 'vwc_m3' in params.keys():
            params['vwc']=params.pop('vwc_m3')
        if 'vwc_w3' in params.keys():
            params['vwc']=params.pop('vwc_w3')
        if 'vwc_s3' in params.keys():
            params['vwc']=params.pop('vwc_s3')            
        
        #rename leaf area/angle index
        if 'ij_f' in params.keys():
            params['ij']=params.pop('ij_f')            
        if 'ij_d' in params.keys():
            params['ij']=params.pop('ij_d')
        if 'ij_m' in params.keys():
            params['ij']=params.pop('ij_m')
        if 'ij_w' in params.keys():
            params['ij']=params.pop('ij_w')
        if 'ij_s' in params.keys():
            params['ij']=params.pop('ij_s')    
        
        params_int+=[params]       
    
    leaf_params+=params_int