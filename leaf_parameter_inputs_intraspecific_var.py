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
nm_f=0.0253 #leaf N (gN/gC)
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
nm_w=0.0435 #leaf N (gN/gC)
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

##plasticity of traits

var_plas=[0.6,1.0,1.4] #variation in leaf trait plasticity

leaf_params_plas=[] #empty list of leaf_params for each type of leaf trait plasiticty

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

for j in var_plas:
    
    s_range=np.random.normal(mean_s,std_s,100000)
    nm_range=np.random.normal(mean_nm,std_nm,100000)
    chl_range=np.random.normal(mean_chl,std_chl,100000)
    tl_range=np.random.normal(mean_tl,std_tl,100000)
    m_range=np.random.normal(mean_m,std_m,100000)
    ij_range=np.random.normal(mean_ij,std_ij*j,100000)
    
    s_f_range=np.array([s_f-(s_m-min(s_range)),s_f+(s_m-min(s_range))])
    nm_f_range=np.array([nm_f-(nm_m-min(nm_range)),nm_f+(nm_m-min(nm_range))])
    chl_f_range=np.array([chl_f-(chl_m-min(chl_range)),chl_f+(chl_m-min(chl_range))])
    tl_f_range=np.array([tl_f-(tl_m-min(tl_range)),tl_f+(tl_m-min(tl_range))])
    m_f_range=np.array([m_f-(m_m-min(m_range)),m_f+(m_m-min(m_range))])
    ij_f_range=np.array([ij_f-(ij_m-min(ij_range)),ij_f+(ij_m-min(ij_range))])

    s_d_range=np.array([s_d-(s_m-min(s_range)),s_d+(s_m-min(s_range))])
    nm_d_range=np.array([nm_d-(nm_m-min(nm_range)),nm_d+(nm_m-min(nm_range))])
    chl_d_range=np.array([chl_d-(chl_m-min(chl_range)),chl_d+(chl_m-min(chl_range))])
    tl_d_range=np.array([tl_d-(tl_m-min(tl_range)),tl_d+(tl_m-min(tl_range))])
    m_d_range=np.array([m_d-(m_m-min(m_range)),m_d+(m_m-min(m_range))])
    ij_d_range=np.array([ij_d-(ij_m-min(ij_range)),ij_d+(ij_m-min(ij_range))])

    s_m_range=np.array([s_m-(s_m-min(s_range)),s_m+(s_m-min(s_range))])
    nm_m_range=np.array([nm_m-(nm_m-min(nm_range)),nm_m+(nm_m-min(nm_range))])
    chl_m_range=np.array([chl_m-(chl_m-min(chl_range)),chl_m+(chl_m-min(chl_range))])
    tl_m_range=np.array([tl_m-(tl_m-min(tl_range)),tl_m+(tl_m-min(tl_range))])
    m_m_range=np.array([m_m-(m_m-min(m_range)),m_m+(m_m-min(m_range))])
    ij_m_range=np.array([ij_m-(ij_m-min(ij_range)),ij_m+(ij_m-min(ij_range))])

    s_w_range=np.array([s_w-(s_m-min(s_range)),s_w+(s_m-min(s_range))])
    nm_w_range=np.array([nm_w-(nm_m-min(nm_range)),nm_w+(nm_m-min(nm_range))])
    chl_w_range=np.array([chl_w-(chl_m-min(chl_range)),chl_w+(chl_m-min(chl_range))])
    tl_w_range=np.array([tl_w-(tl_m-min(tl_range)),tl_w+(tl_m-min(tl_range))])
    m_w_range=np.array([m_w-(m_m-min(m_range)),m_w+(m_m-min(m_range))])
    ij_w_range=np.array([ij_w-(ij_m-min(ij_range)),ij_w+(ij_m-min(ij_range))])
    
    s_s_range=np.array([s_s-(s_m-min(s_range)),s_s+(s_m-min(s_range))])
    nm_s_range=np.array([nm_s-(nm_m-min(nm_range)),nm_s+(nm_m-min(nm_range))])
    chl_s_range=np.array([chl_s-(chl_m-min(chl_range)),chl_s+(chl_m-min(chl_range))])
    tl_s_range=np.array([tl_s-(tl_m-min(tl_range)),tl_s+(tl_m-min(tl_range))])
    m_s_range=np.array([m_s-(m_m-min(m_range)),m_s+(m_m-min(m_range))])
    ij_s_range=np.array([ij_s-(ij_m-min(ij_range)),ij_s+(ij_m-min(ij_range))])

##---In the Below Dictionary, sla, leaf n, and fraction of leaf n in rub are the only proportionally dependent variables---##  


#traits vary for each community

    dict=[
      [{'s_f':np.linspace(min(s_f_range),max(s_f_range),2)},
       {'nm_f':np.linspace(min(nm_f_range),max(nm_f_range),2)},
       {'m_f1':np.zeros(shape=2)+min(m_f_range)},
       {'m_f2':np.zeros(shape=2)+max(m_f_range)},
       {'chl_f1':np.zeros(shape=2)+min(chl_f_range)},
       {'chl_f2':np.zeros(shape=2)+max(chl_f_range)},
       {'tl_f1':np.zeros(shape=2)+min(tl_f_range)},
       {'tl_f2':np.zeros(shape=2)+max(tl_f_range)},
       {'vwc_f1':np.zeros(shape=2)+vwc_m},
       {'vwc_f2':np.zeros(shape=2)+vwc_m},
       {'ij_f1':np.zeros(shape=2)+min(ij_f_range)},
       {'ij_f2':np.zeros(shape=2)+max(ij_f_range)},
        ],
      
      [{'s_d':np.linspace(min(s_d_range),max(s_d_range),2)},
       {'nm_d':np.linspace(min(nm_d_range),max(nm_d_range),2)},
       {'m_d1':np.zeros(shape=2)+min(m_d_range)},
       {'m_d2':np.zeros(shape=2)+max(m_d_range)},
       {'chl_d1':np.zeros(shape=2)+min(chl_d_range)},
       {'chl_d2':np.zeros(shape=2)+max(chl_d_range)},
       {'tl_d1':np.zeros(shape=2)+min(tl_d_range)},
       {'tl_d2':np.zeros(shape=2)+max(tl_d_range)},
       {'vwc_d1':np.zeros(shape=2)+vwc_m},
       {'vwc_d2':np.zeros(shape=2)+vwc_m},
        {'ij_d1':np.zeros(shape=2)+min(ij_d_range)},
        {'ij_d2':np.zeros(shape=2)+max(ij_d_range)},
        ],
        
      [{'s_m':np.linspace(min(s_m_range),max(s_m_range),2)},
       {'nm_m':np.linspace(min(nm_m_range),max(nm_m_range),2)},
       {'m_m1':np.zeros(shape=2)+min(m_m_range)},
       {'m_m2':np.zeros(shape=2)+max(m_m_range)},
       {'chl_m1':np.zeros(shape=2)+min(chl_m_range)},
       {'chl_m2':np.zeros(shape=2)+max(chl_m_range)},
       {'tl_m1':np.zeros(shape=2)+min(tl_m_range)},
       {'tl_m2':np.zeros(shape=2)+max(tl_m_range)},
       {'vwc_m1':np.zeros(shape=2)+vwc_m},
       {'vwc_m2':np.zeros(shape=2)+vwc_m},
       {'ij_m1':np.zeros(shape=2)+min(ij_m_range)},
       {'ij_m2':np.zeros(shape=2)+max(ij_m_range)},
        ],
       
      [{'s_w':np.linspace(min(s_w_range),max(s_w_range),2)},
       {'nm_w':np.linspace(min(nm_w_range),max(nm_w_range),2)},
       {'m_w1':np.zeros(shape=2)+min(m_w_range)},
       {'m_w2':np.zeros(shape=2)+max(m_w_range)},
       {'chl_w1':np.zeros(shape=2)+min(chl_w_range)},
       {'chl_w2':np.zeros(shape=2)+max(chl_w_range)},
       {'tl_w1':np.zeros(shape=2)+min(tl_w_range)},
       {'tl_w2':np.zeros(shape=2)+max(tl_w_range)},
       {'vwc_w1':np.zeros(shape=2)+vwc_m},
       {'vwc_w2':np.zeros(shape=2)+vwc_m},
       {'ij_w1':np.zeros(shape=2)+min(ij_w_range)},
       {'ij_w2':np.zeros(shape=2)+max(ij_w_range)},
        ], 
      
      [{'s_s':np.linspace(min(s_s_range),max(s_s_range),2)},
       {'nm_s':np.linspace(min(nm_s_range),max(nm_s_range),2)},
       {'m_s1':np.zeros(shape=2)+min(m_s_range)},
       {'m_s2':np.zeros(shape=2)+max(m_s_range)},
       {'chl_s1':np.zeros(shape=2)+min(chl_s_range)},
       {'chl_s2':np.zeros(shape=2)+max(chl_s_range)},
       {'tl_s1':np.zeros(shape=2)+min(tl_s_range)},
       {'tl_s2':np.zeros(shape=2)+max(tl_s_range)},
       {'vwc_s1':np.zeros(shape=2)+vwc_m},
       {'vwc_s2':np.zeros(shape=2)+vwc_m},
       {'ij_s1':np.zeros(shape=2)+min(ij_s_range)},
       {'ij_s2':np.zeros(shape=2)+max(ij_s_range)},
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
    
    leaf_params_plas+=[leaf_params]