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
     
#---------------Create Dictionary For All Types of Parameter Inputs---------------#  

dict=[
      {'s': np.zeros(shape=100)+0.02},
      {'s_in':np.linspace(0.01,0.06,100)},
#      {'s_de':np.linspace(0.06,0.01,100)},
      {'ra':np.zeros(shape=100)+20.7}, 
      {'ra_in':np.linspace(17,23,100)},
#      {'ra_de':np.linspace(23,17,100)},
      {'nm':np.zeros(shape=100)+0.045},
      {'nm_in':np.linspace(0.04,0.06,100)},
#      {'nm_de':np.linspace(0.06,0.04,100)},
      {'flnr':np.zeros(shape=100)+0.7},
      {'flnr_in':np.linspace(0.5,0.9,100)},
#      {'flnr_de':np.linspace(0.9,0.5,100)},
#      {'m':np.zeros(shape=100)+9},
#      {'m_in':np.linspace(7,11,100)},
#      {'m_de':np.linspace(11,7,100)},
#      {'b':np.zeros(shape=100)+0.01},
#      {'b_in':np.linspace(0,0.02,100)},
#      {'b_de':np.linspace(0.02,0,100)}
        ]

#---------------Define Number of Variable Parameters---------------#  

num_params=4


#---------------Make Dictionaries that Contain All Combinations of Parameter Inputs---------------#  

dict_combs_raw=[]
for x in it.combinations(dict,num_params):
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


leaf_params=[]

for i in range(len(dict_combs_new)):
    params=dict_combs_new[i].copy()

    if 's_in' in params.keys():
        params['s']=params.pop('s_in')
    if 's_de' in params.keys():
        params['s']=params.pop('s_de')
    
    if 'ra_in' in params.keys():
        params['ra']=params.pop('ra_in')
    if 'ra_de' in params.keys():
        params['ra']=params.pop('ra_de')
        
    if 'nm_in' in params.keys():
        params['nm']=params.pop('nm_in')
    if 'nm_de' in params.keys():
        params['nm']=params.pop('nm_de')
    
    if 'flnr_in' in params.keys():
        params['flnr']=params.pop('flnr_in')
    if 'flnr_de' in params.keys():
        params['flnr']=params.pop('flnr_de')
        
    if 'm_in' in params.keys():
        params['m']=params.pop('m_in')
    if 'm_de' in params.keys():
        params['m']=params.pop('m_de')
    
        
    if 'b_in' in params.keys():
        params['b']=params.pop('b_in')
    if 'b_de' in params.keys():
        params['b']=params.pop('b_de')
    
    leaf_params+=[params]       