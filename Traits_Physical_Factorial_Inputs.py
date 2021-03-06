#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May 24 20:31:57 2017

@author: Katherine
"""

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



#note to self: alternative relationship between sla and N: 0.61*log(sla)+0.55=log(N)

#log(5)=1.6

#exp(1.6)=5

#np.exp(0.61*log(sla)+0.55)=N


    
#---------------Create Dictionary For All Variable Parameter Inputs---------------#  

ekc_lw=80500.0*0.8 #Activation energy for K of CO2 (J mol-1) VARIABLE
ekc_md=80500.0 #Activation energy for K of CO2 (J mol-1) VARIABLE
ekc_hg=80500.0*1.2 #Activation energy for K of CO2 (J mol-1) VARIABLE

eko_lw=14500.0*0.8 #Activation energy for K of O2 (J mol-1) VARIABLE
eko_md=14500.0 #Activation energy for K of O2 (J mol-1) VARIABLE
eko_hg=14500.0*1.2 #Activation energy for K of O2 (J mol-1) VARIABLE


etau_lw=-29000.0*0.8  #Activation energy for tau (???) (J mol-1) VARIABLE
etau_md=-29000.0  #Activation energy for tau (???) (J mol-1) VARIABLE
etau_hg=-29000.0*1.2  #Activation energy for tau (???) (J mol-1) VARIABLE


ev_lw=55000.0*0.8 #Activation energy for carboxylation (J mol-1) VARIABLE
ev_md=55000.0 #Activation energy for carboxylation (J mol-1) VARIABLE
ev_hg=55000.0*1.2 #Activation energy for carboxylation (J mol-1) VARIABLE

ej_lw=55000.0*0.8 #Activation energy for electron transport (J mol-1) VARIABLE
ej_md=55000.0 #Activation energy for electron transport (J mol-1) VARIABLE
ej_hg=55000.0*1.2 #Activation energy for electron transport (J mol-1) VARIABLE

ra_lw=20.7*0.8 #specific rubisco activity (umol CO2/g Rub s) VARIABLE
ra_md=20.7 #specific rubisco activity (umol CO2/g Rub s) VARIABLE
ra_hg=20.7*1.2 #specific rubisco activity (umol CO2/g Rub s) VARIABLE

flnr_lw=0.1*0.8 #fraction of leaf nitrogen in rubisco (g N Rub/g N leaf) VARIABLE
flnr_md=0.1 #fraction of leaf nitrogen in rubisco (g N Rub/g N leaf) VARIABLE
flnr_hg=0.1*1.2 #fraction of leaf nitrogen in rubisco (g N Rub/g N leaf) VARIABLE

rh_lw=0.5*0.8 #relative humidity (kPa/kPa) VARIABLE
rh_md=0.5 #relative humidity (kPa/kPa) VARIABLE
rh_hg=0.5*1.2 #relative humidity (kPa/kPa) VARIABLE


ko25_lw=30000*0.8 #Michaelis-Menten kinetic coefficient for oxygen at 25 C(Pa) VARIABLE
ko25_md=30000 #Michaelis-Menten kinetic coefficient for oxygen at 25 C(Pa) VARIABLE
ko25_hg=30000*1.2 #Michaelis-Menten kinetic coefficient for oxygen at 25 C(Pa) VARIABLE


kc25_lw=30*0.8 #Michaelis-Menten kinetic coefficient for carbon dioxide at 25 C (Pa) VARIABLE
kc25_md=30 #Michaelis-Menten kinetic coefficient for carbon dioxide at 25 C (Pa) VARIABLE
kc25_hg=30*1.2 #Michaelis-Menten kinetic coefficient for carbon dioxide at 25 C (Pa) VARIABLE


g0_lw=0.002*0.8 #Ball-Berry stomatal conductance intercept parameter (mol H2O/m2s) VARIABLE
g0_md=0.002 #Ball-Berry stomatal conductance intercept parameter (mol H2O/m2s) VARIABLE
g0_hg=0.002*1.2 #Ball-Berry stomatal conductance intercept parameter (mol H2O/m2s) VARIABLE


m_lw=9.0*0.8 #ball-berry parameter (unitless) VARIABLE
m_md=9.0 #ball-berry parameter (unitless) VARIABLE
m_hg=9.0*1.2 #ball-berry parameter (unitless) VARIABLE


u_lw=5.0*0.8 #windspeed (m/s) VARIABLE
u_md=5.0 #windspeed (m/s) VARIABLE
u_hg=5.0*1.2 #windspeed (m/s) VARIABLE

qeff_lw=0.32*0.8 #leaf quantum yield, electrons VARIABLE
qeff_md=0.32 #leaf quantum yield, electrons VARIABLE
qeff_hg=0.32*1.2 #leaf quantum yield, electrons VARIABLE

PAR_lw=2000*0.8 #photosynthetic active radiation (umol/m2s) VARIABLE
PAR_md=2000 #photosynthetic active radiation (umol/m2s) VARIABLE
PAR_hg=2000*1.2 #photosynthetic active radiation (umol/m2s) VARIABLE


jm_lw=2.68*0.8 #slope coefficient  VARIABLE
jm_md=2.68 #slope coefficient  VARIABLE
jm_hg=2.68*1.2 #slope coefficient  VARIABLE


vwc_min_lw=0.08*0.8 #minimum soil water content for photosynthesis to occur (permanent wilting point) (cm3/cm3)  VARIABLE
vwc_min_md=0.08 #minimum soil water content for photosynthesis to occur (permanent wilting point) (cm3/cm3)  VARIABLE
vwc_min_hg=0.08*1.2 #minimum soil water content for photosynthesis to occur (permanent wilting point) (cm3/cm3)  VARIABLE

vwc_max_lw=0.68*0.8 #maximum soil water content where increases in soil water do not affect photosynthesis (field capacity?) (cm3/cm3) VARIABLE
vwc_max_md=0.68 #maximum soil water content where increases in soil water do not affect photosynthesis (field capacity?) (cm3/cm3) VARIABLE
vwc_max_hg=0.68*1.2 #maximum soil water content where increases in soil water do not affect photosynthesis (field capacity?) (cm3/cm3) VARIABLE


q_lw=0.2*0.8 #parameter for soil water affect on photosynthesis (unitless) VARIABLE
q_md=0.2 #parameter for soil water affect on photosynthesis (unitless) VARIABLE
q_hg=0.2*1.2 #parameter for soil water affect on photosynthesis (unitless) VARIABLE



##mean trait assemblage parameters for dry meadow, moist meadow, wet meadow, snowbed 

#fellfield

#chl_f1=275.4669 #chlorophyll content (umol chl/m2)
#chl_f2=322.7893 #chlorophyll content (umol chl/m2)
#chl_f3=370.1116 #chlorophyll content (umol chl/m2)
#sla_f1=98.49966-10.6
#sla_f2=98.49966
#sla_f3=98.49966+10.6
#nm_f1=0.02 #(20.25+0.077*(sla_f1))/1000.0
#nm_f2=0.02 #(20.25+0.077*(sla_f2))/1000.0
#nm_f3=0.02 #(20.25+0.077*(sla_f3))/1000.0
#na_f1=2.0#nm_f1*(1/(sla_f1/10000))#(((chl_f1)/1000.0-0.022)/0.0036)/1000.0*14.0#nm_f1*(1/(sla_f1/10000))  #leaf N (gN/m2C)
#na_f2=2.0#nm_f2*(1/(sla_f2/10000))#(((chl_f2)/1000.0-0.022)/0.0036)/1000.0*14.0#nm_f2*(1/(sla_f2/10000))   #leaf N (gN/m2C)
#na_f3=2.0#nm_f3*(1/(sla_f3/10000)) #(((chl_f3)/1000.0-0.022)/0.0036)/1000.0*14.0#nm_f3*(1/(sla_f3/10000))  #leaf N (gN/m2C)
#ht_f1=5.423279#t+16 #temp of leaf (K)-->f(leaf area,angle,height)
#ht_f2=6.324087#t+16 #temp of leaf (K)-->f(leaf area,angle,height)
#ht_f3=7.224895#t+16 #temp of leaf (K)-->f(leaf area,angle,height)
#dia_f1=0.9141305/100. #leaf diameter (m)
#dia_f2=1.144449/100. #leaf diameter (m)
#dia_f3=1.3356245/100. #leaf diameter (m)

#dry meadow
chl_d1=371.3030 #chlorophyll content (umol chl/m2)
chl_d2=395.7132 #chlorophyll content (umol chl/m2)
chl_d3=420.1234 #chlorophyll content (umol chl/m2)
sla_d1=67.86595-7.6
sla_d2=67.86595
sla_d3=67.86595+7.6
#nm_d1=0.025#(20.25+0.077*(sla_d1))/1000.0
#nm_d2=0.025#(20.25+0.077*(sla_d2))/1000.0
#nm_d3=0.025#(20.25+0.077*(sla_d3))/1000.0
na_d1=2.45-0.6#nm_d1*(1/(sla_d1/10000))#(((chl_d1)/1000.0-0.022)/0.0036)/1000.0*14.0#nm_d1*(1/(sla_d1/10000)) #leaf N (gN/m2C)
na_d2=2.45#nm_d2*(1/(sla_d2/10000))#(((chl_d2)/1000.0-0.022)/0.0036)/1000.0*14.0#nm_d2*(1/(sla_d2/10000)) #leaf N (gN/m2C)
na_d3=2.45+0.6#nm_d3*(1/(sla_d3/10000))#(((chl_d3)/1000.0-0.022)/0.0036)/1000.0*14.0#nm_d3*(1/(sla_d3/10000)) #leaf N (gN/m2C)
ht_d1=7.688147#height of leaf (cm)
ht_d2=9.183549##height of leaf (cm)
ht_d3=10.678950#height of leaf (cm)
dia_d1=1.177997/100. #leaf diameter (m)
dia_d2=1.6/100. #leaf diameter (m)
dia_d3=1.642849/100. #leaf diameter (m)
t_d1=(15+2)*0.8
t_d2=15+2
t_d3=(15+2)*1.2
vwc_d1=0.12*0.8
vwc_d2=0.12
vwc_d3=0.12*1.2


#moist meadow
chl_m1=411.0208 #chlorophyll content (umol chl/m2)
chl_m2=464.606 #chlorophyll content (umol chl/m2)
chl_m3=518.1911 #chlorophyll content (umol chl/m2)
sla_m1=111.947-12.3
sla_m2=111.947
sla_m3=111.947+12.3
#nm_m1=0.03#(20.25+0.077*(sla_m1))/1000.0
#nm_m2=0.03#(20.25+0.077*(sla_m2))/1000.0
#nm_m3=0.03#(20.25+0.077*(sla_m3))/1000.0
na_m1=5.0-0.79#nm_m1*(1/(sla_m1/10000)) #nm_m1*(1/(sla_m1/10000))#(((chl_m1)/1000.0-0.022)/0.0036)/1000.0*14.0 #nm_m1*(1/(sla_m1/10000)) #leaf N (gN/m2C)
na_m2=5.0#nm_m2*(1/(sla_m2/10000)) #(((chl_m2)/1000.0-0.022)/0.0036)/1000.0*14.0 #nm_m2*(1/(sla_m2/10000)) #leaf N (gN/m2C)
na_m3=5.0+0.79#nm_m3*(1/(sla_m3/10000)) #(((chl_m3)/1000.0-0.022)/0.0036)/1000.0*14.0 #nm_m3*(1/(sla_m3/10000)) #leaf N (gN/m2C)
ht_m1=13.48271#height of leaf (cm)
ht_m2=19.19779#height of leaf (cm)
ht_m3=24.91287#height of leaf (cm)
dia_m1=2.042256/100. #leaf diameter (m)
dia_m2=2.6/100. #leaf diameter (m)
dia_m3=2.605534/100. #leaf diameter (m)
t_m1=(15)*0.8
t_m2=15
t_m3=(15)*1.2
vwc_m1=0.15*0.8
vwc_m2=0.15
vwc_m3=0.15*1.2


##wet meadow
chl_w1=446.7062 #chlorophyll content (umol chl/m2)
chl_w2=475.8913 #chlorophyll content (umol chl/m2)
chl_w3=505.0765 #chlorophyll content (umol chl/m2)
sla_w1=117.6779-5.7
sla_w2=117.6779
sla_w3=117.6779+5.7
#nm_w1=0.035#(20.25+0.077*(sla_w1))/1000.0
#nm_w2=0.035#(20.25+0.077*(sla_w2))/1000.0
#nm_w3=0.035#(20.25+0.077*(sla_w3))/1000.0
na_w1=6.25-0.61#nm_w1*(1/(sla_w1/10000))#(((chl_w1)/1000.0-0.022)/0.0036)/1000.0*14.0#nm_w1*(1/(sla_w1/10000)) #leaf N (gN/m2C)
na_w2=6.25#nm_w2*(1/(sla_w2/10000))#(((chl_w2)/1000.0-0.022)/0.0036)/1000.0*14.0#nm_w2*(1/(sla_w2/10000)) #leaf N (gN/m2C)
na_w3=6.25+0.61#nm_w3*(1/(sla_w3/10000))#(((chl_w3)/1000.0-0.022)/0.0036)/1000.0*14.0#nm_w3*(1/(sla_w3/10000)) #leaf N (gN/m2C)
ht_w1=16.86773#height of leaf (cm)
ht_w2=19.98519#height of leaf (cm)
ht_w3=23.10266#height of leaf (cm)
dia_w1=2.378119/100. #leaf diameter (m)
dia_w2=3.0/100. #leaf diameter (m)
dia_w3=2.846025/100. #leaf diameter (m)
t_w1=(15-2)*0.8
t_w2=15-2
t_w3=(15-2)*1.2
vwc_w1=0.17*0.8
vwc_w2=0.27
vwc_w3=0.27*1.2


##snowbed
#chl_s1=376.2863 #chlorophyll content (umol chl/m2)
#chl_s2=426.6331 #chlorophyll content (umol chl/m2)
#chl_s3=476.9798 #chlorophyll content (umol chl/m2)
#sla_s1=147.195-12.6
#sla_s2=147.195
#sla_s3=147.195+12.6
#nm_s1=0.021#(20.25+0.077*(sla_s1))/1000.0
#nm_s2=0.021#(20.25+0.077*(sla_s2))/1000.0
#nm_s3=0.021#(20.25+0.077*(sla_s3))/1000.0
#na_s1=#nm_s1*(1/(sla_s1/10000))#(((chl_s1)/1000.0-0.022)/0.0036)/1000.0*14.0#nm_s1*(1/(sla_s1/10000)) #leaf N (gN/m2C)
#na_s2=#nm_s2*(1/(sla_s2/10000))#(((chl_s2)/1000.0-0.022)/0.0036)/1000.0*14.0#nm_s2*(1/(sla_s2/10000)) #leaf N (gN/m2C)
#na_s3=#nm_s3*(1/(sla_s3/10000))#(((chl_s3)/1000.0-0.022)/0.0036)/1000.0*14.0#nm_s3*(1/(sla_s3/10000)) #leaf N (gN/m2C)
#ht_s1=9.055197#t #temp of leaf (K)-->f(leaf area,angle,height)
#ht_s2=12.84922#t #temp of leaf (K)-->f(leaf area,angle,height)
#ht_s3=16.643243#t #temp of leaf (K)-->f(leaf area,angle,height)
#dia_s1=1.979317/100. #leaf diameter (m)
#dia_s2=2.136071/100. #leaf diameter (m)
#dia_s3=2.282084/100. #leaf diameter (m)



##---In the Below Dictionary, sla, leaf n, and fraction of leaf n in rub are the only proportionally dependent variables---##



#traits vary for each community


        
        
leaf_params=[
      

       {'na':np.zeros(shape=1)+na_d2,
       'dia':np.zeros(shape=1)+dia_d2,
       'chl':np.zeros(shape=1)+chl_d2,
       'ht':np.zeros(shape=1)+ht_d2,
       'vwc':np.zeros(shape=1)+vwc_d2,
       'temp':np.zeros(shape=1)+t_d2},
       {'na':np.zeros(shape=1)+na_m2,
       'dia':np.zeros(shape=1)+dia_m2,
       'chl':np.zeros(shape=1)+chl_m2,
       'ht':np.zeros(shape=1)+ht_m2,
       'vwc':np.zeros(shape=1)+vwc_d2,
       'temp':np.zeros(shape=1)+t_d2},
       {'na':np.zeros(shape=1)+na_w2,
       'dia':np.zeros(shape=1)+dia_w2,
       'chl':np.zeros(shape=1)+chl_w2,
       'ht':np.zeros(shape=1)+ht_w2,
       'vwc':np.zeros(shape=1)+vwc_d2,
       'temp':np.zeros(shape=1)+t_d2},
       {'na':np.zeros(shape=1)+na_d2,
       'dia':np.zeros(shape=1)+dia_d2,
       'chl':np.zeros(shape=1)+chl_d2,
       'ht':np.zeros(shape=1)+ht_d2,
       'vwc':np.zeros(shape=1)+vwc_m2,
       'temp':np.zeros(shape=1)+t_m2},
       {'na':np.zeros(shape=1)+na_m2,
       'dia':np.zeros(shape=1)+dia_m2,
       'chl':np.zeros(shape=1)+chl_m2,
       'ht':np.zeros(shape=1)+ht_m2,
       'vwc':np.zeros(shape=1)+vwc_m2,
       'temp':np.zeros(shape=1)+t_m2},
       {'na':np.zeros(shape=1)+na_w2,
       'dia':np.zeros(shape=1)+dia_w2,
       'chl':np.zeros(shape=1)+chl_w2,
       'ht':np.zeros(shape=1)+ht_w2,
       'vwc':np.zeros(shape=1)+vwc_m2,
       'temp':np.zeros(shape=1)+t_m2},
       {'na':np.zeros(shape=1)+na_d2,
       'dia':np.zeros(shape=1)+dia_d2,
       'chl':np.zeros(shape=1)+chl_d2,
       'ht':np.zeros(shape=1)+ht_d2,
       'vwc':np.zeros(shape=1)+vwc_w2,
       'temp':np.zeros(shape=1)+t_w2},
       {'na':np.zeros(shape=1)+na_m2,
       'dia':np.zeros(shape=1)+dia_m2,
       'chl':np.zeros(shape=1)+chl_m2,
       'ht':np.zeros(shape=1)+ht_m2,
       'vwc':np.zeros(shape=1)+vwc_w2,
       'temp':np.zeros(shape=1)+t_w2},
       {'na':np.zeros(shape=1)+na_w2,
       'dia':np.zeros(shape=1)+dia_w2,
       'chl':np.zeros(shape=1)+chl_w2,
       'ht':np.zeros(shape=1)+ht_w2,
       'vwc':np.zeros(shape=1)+vwc_w2,
       'temp':np.zeros(shape=1)+t_w2},
        ]        
        
##---------------Define Number of Variable Parameters---------------#  
#
#num_params=4
#
#leaf_params=[]
#
##---------------Run For Loop of Different Dictionary Values---------------# 
#
#for i in range(len(dict)):
#
##---------------Make Dictionaries that Contain All Combinations of Parameter Inputs---------------#  
#
#    dict_combs_raw=[]
#    for x in it.combinations(dict[i],num_params):
#        dict_combs_raw+=[x]
#
#    #merge dictionaries 
#    dict_combs=[]
#    for x in range(len(dict_combs_raw)):
#        dict_combs+=[merge_dicts(dict_combs_raw[x])]
# 
#
##---------------Get Rid of Repeated Parameters (i.e. s_in & s) in Parameter Combination Dictionaries---------------#  
#                 
#    dict_combs_new=[]
#
#    for xx in range(len(dict_combs)):
#        keys=dict_combs[xx].keys()
#        key_param=[]
#        for ii in range(len(keys)):
#            key_param+=[keys[ii][:-3]]
#        if len(np.unique(key_param))<num_params:
#            continue
#        else:
#            dict_combs_new+=[dict_combs[xx]]
#
##---------------Rename All Parameters in Parameter Combination Dictionaries to Model Parameter Names---------------#  
#
#    params_int=[]
#
#
#    for xxx in range(len(dict_combs_new)):
#        params=dict_combs_new[xxx].copy()
#        
#             
#        
#        #rename leaf n
#        if 'na_f1' in params.keys():
#            params['na']=params.pop('na_f1')  
#        if 'na_f2' in params.keys():
#            params['na']=params.pop('na_f2')  
#        if 'na_f3' in params.keys():
#            params['na']=params.pop('na_f3')  
#            
#        if 'na_d1' in params.keys():
#            params['na']=params.pop('na_d1')
#        if 'na_d2' in params.keys():
#            params['na']=params.pop('na_d2')
#        if 'na_d3' in params.keys():
#            params['na']=params.pop('na_d3')
#            
#        if 'na_m1' in params.keys():
#            params['na']=params.pop('na_m1')
#        if 'na_m2' in params.keys():
#            params['na']=params.pop('na_m2')
#        if 'na_m3' in params.keys():
#            params['na']=params.pop('na_m3')
#            
#        if 'na_w1' in params.keys():
#            params['na']=params.pop('na_w1')
#        if 'na_w2' in params.keys():
#            params['na']=params.pop('na_w2')
#        if 'na_w3' in params.keys():
#            params['na']=params.pop('na_w3')
#            
#        if 'na_s1' in params.keys():
#            params['na']=params.pop('na_s1') 
#        if 'na_s2' in params.keys():
#            params['na']=params.pop('na_s2')  
#        if 'na_s3' in params.keys():
#            params['na']=params.pop('na_s3')  
#        
#            
#        #rename diameter
#        if 'dia_f1' in params.keys():
#            params['dia']=params.pop('dia_f1')
#        if 'dia_d1' in params.keys():
#            params['dia']=params.pop('dia_d1')
#        if 'dia_m1' in params.keys():
#            params['dia']=params.pop('dia_m1')
#        if 'dia_w1' in params.keys():
#            params['dia']=params.pop('dia_w1')            
#        if 'dia_s1' in params.keys():
#            params['dia']=params.pop('dia_s1')  
#
#        if 'dia_f2' in params.keys():
#            params['dia']=params.pop('dia_f2')
#        if 'dia_d2' in params.keys():
#            params['dia']=params.pop('dia_d2')
#        if 'dia_m2' in params.keys():
#            params['dia']=params.pop('dia_m2')
#        if 'dia_w2' in params.keys():
#            params['dia']=params.pop('dia_w2')            
#        if 'dia_s2' in params.keys():
#            params['dia']=params.pop('dia_s2')         
#
#        if 'dia_f3' in params.keys():
#            params['dia']=params.pop('dia_f3')
#        if 'dia_d3' in params.keys():
#            params['dia']=params.pop('dia_d3')
#        if 'dia_m3' in params.keys():
#            params['dia']=params.pop('dia_m3')
#        if 'dia_w3' in params.keys():
#            params['dia']=params.pop('dia_w3')            
#        if 'dia_s3' in params.keys():
#            params['dia']=params.pop('dia_s3')                  
#
#        
#        #rename chlorophyll content
#        if 'chl_f1' in params.keys():
#            params['chl']=params.pop('chl_f1')        
#        if 'chl_d1' in params.keys():
#            params['chl']=params.pop('chl_d1')
#        if 'chl_m1' in params.keys():
#            params['chl']=params.pop('chl_m1')
#        if 'chl_w1' in params.keys():
#            params['chl']=params.pop('chl_w1')
#        if 'chl_s1' in params.keys():
#            params['chl']=params.pop('chl_s1')            
#
#        if 'chl_f2' in params.keys():
#            params['chl']=params.pop('chl_f2')            
#        if 'chl_d2' in params.keys():
#            params['chl']=params.pop('chl_d2')
#        if 'chl_m2' in params.keys():
#            params['chl']=params.pop('chl_m2')
#        if 'chl_w2' in params.keys():
#            params['chl']=params.pop('chl_w2')
#        if 'chl_s2' in params.keys():
#            params['chl']=params.pop('chl_s2')              
# 
#        if 'chl_f3' in params.keys():
#            params['chl']=params.pop('chl_f3')            
#        if 'chl_d3' in params.keys():
#            params['chl']=params.pop('chl_d3')
#        if 'chl_m3' in params.keys():
#            params['chl']=params.pop('chl_m3')
#        if 'chl_w3' in params.keys():
#            params['chl']=params.pop('chl_w3')
#        if 'chl_s3' in params.keys():
#            params['chl']=params.pop('chl_s3')  
#       
#        #rename temp of leaf
#        if 'ht_f1' in params.keys():
#            params['ht']=params.pop('ht_f1')        
#        if 'ht_d1' in params.keys():
#            params['ht']=params.pop('ht_d1')
#        if 'ht_m1' in params.keys():
#            params['ht']=params.pop('ht_m1')
#        if 'ht_w1' in params.keys():
#            params['ht']=params.pop('ht_w1')
#        if 'ht_s1' in params.keys():
#            params['ht']=params.pop('ht_s1')            
#
#        if 'ht_f2' in params.keys():
#            params['ht']=params.pop('ht_f2')            
#        if 'ht_d2' in params.keys():
#            params['ht']=params.pop('ht_d2')
#        if 'ht_m2' in params.keys():
#            params['ht']=params.pop('ht_m2')
#        if 'ht_w2' in params.keys():
#            params['ht']=params.pop('ht_w2')
#        if 'ht_s2' in params.keys():
#            params['ht']=params.pop('ht_s2')              
#                
#        if 'ht_f3' in params.keys():
#            params['ht']=params.pop('ht_f3')            
#        if 'ht_d3' in params.keys():
#            params['ht']=params.pop('ht_d3')
#        if 'ht_m3' in params.keys():
#            params['ht']=params.pop('ht_m3')
#        if 'ht_w3' in params.keys():
#            params['ht']=params.pop('ht_w3')
#        if 'ht_s3' in params.keys():
#            params['ht']=params.pop('ht_s3')   
#            
#       #rename leaf sla
#        if 'sla_f1' in params.keys():
#            params['sla']=params.pop('sla_f1')  
#        if 'sla_f2' in params.keys():
#            params['sla']=params.pop('sla_f2')  
#        if 'sla_f3' in params.keys():
#            params['sla']=params.pop('sla_f3')  
#            
#        if 'sla_d1' in params.keys():
#            params['sla']=params.pop('sla_d1')
#        if 'sla_d2' in params.keys():
#            params['sla']=params.pop('sla_d2')
#        if 'sla_d3' in params.keys():
#            params['sla']=params.pop('sla_d3')
#            
#        if 'sla_m1' in params.keys():
#            params['sla']=params.pop('sla_m1')
#        if 'sla_m2' in params.keys():
#            params['sla']=params.pop('sla_m2')
#        if 'sla_m3' in params.keys():
#            params['sla']=params.pop('sla_m3')
#            
#        if 'sla_w1' in params.keys():
#            params['sla']=params.pop('sla_w1')
#        if 'sla_w2' in params.keys():
#            params['sla']=params.pop('sla_w2')
#        if 'sla_w3' in params.keys():
#            params['sla']=params.pop('sla_w3')
#            
#        if 'sla_s1' in params.keys():
#            params['sla']=params.pop('sla_s1') 
#        if 'sla_s2' in params.keys():
#            params['sla']=params.pop('sla_s2')  
#        if 'sla_s3' in params.keys():
#            params['sla']=params.pop('sla_s3')  
#
#
#
#       #rename leaf temp
#            
#        if 't_d1' in params.keys():
#            params['t']=params.pop('t_d1')
#        if 't_d2' in params.keys():
#            params['t']=params.pop('t_d2')
#        if 't_d3' in params.keys():
#            params['t']=params.pop('t_d3')
#            
#        if 't_m1' in params.keys():
#            params['t']=params.pop('t_m1')
#        if 't_m2' in params.keys():
#            params['t']=params.pop('t_m2')
#        if 't_m3' in params.keys():
#            params['t']=params.pop('t_m3')
#            
#        if 't_w1' in params.keys():
#            params['t']=params.pop('t_w1')
#        if 't_w2' in params.keys():
#            params['t']=params.pop('t_w2')
#        if 't_w3' in params.keys():
#            params['t']=params.pop('t_w3')
#  
#
#       #rename vwc
#            
#        if 'vwc_d1' in params.keys():
#            params['vwc']=params.pop('vwc_d1')
#        if 'vwc_d2' in params.keys():
#            params['vwc']=params.pop('vwc_d2')
#        if 'vwc_d3' in params.keys():
#            params['vwc']=params.pop('vwc_d3')
#            
#        if 'vwc_m1' in params.keys():
#            params['vwc']=params.pop('vwc_m1')
#        if 'vwc_m2' in params.keys():
#            params['vwc']=params.pop('vwc_m2')
#        if 'vwc_m3' in params.keys():
#            params['vwc']=params.pop('vwc_m3')
#            
#        if 'vwc_w1' in params.keys():
#            params['vwc']=params.pop('vwc_w1')
#        if 'vwc_w2' in params.keys():
#            params['vwc']=params.pop('vwc_w2')
#        if 'vwc_w3' in params.keys():
#            params['vwc']=params.pop('vwc_w3')        
#            
#            
#
#
#        params_int+=[params]       
#    
# 
#    leaf_params+=params_int