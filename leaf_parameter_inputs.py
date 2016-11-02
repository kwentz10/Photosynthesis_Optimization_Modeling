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

##---In the Below Dictionary, Stomatal Regulation is NOT a Function of Carboxylation Efficiency---##
dict=[
      [{'s_h':np.linspace(0.03,0.0309,3)},
       {'ra_h':np.linspace(20.7,21.0,3)},
       {'nm_h':np.linspace(0.04,0.05,3)},
       {'flnr_h':np.linspace(0.625,0.675,3)},
       {'m_h1':np.zeros(shape=3)+10.25},
       {'m_h2':np.zeros(shape=3)+10.75},
       {'m_l1':np.zeros(shape=3)+7.25},
       {'m_l2':np.zeros(shape=3)+7.75}],
      [{'s_l':np.linspace(0.01,0.0109,3)},
       {'ra_l':np.linspace(19.7,20,3)},
       {'nm_l':np.linspace(0.02,0.03,3)},
       {'flnr_l':np.linspace(0.425,0.475,3)},
       {'m_h1':np.zeros(shape=3)+10.25},
       {'m_h2':np.zeros(shape=3)+10.75},
       {'m_l1':np.zeros(shape=3)+7.25},
       {'m_l2':np.zeros(shape=3)+7.75}]
        ]

##---In the Below Dictionary, Stomatal Regulation is a Function of Carboxylation Efficiency---##

#dict=[
#      [{'s_h':np.zeros(shape=4)+0.06},
#       {'ra_h':np.zeros(shape=4)+20.7},
#       {'nm_h':np.linspace(0.06,0.08,4)},
#       {'flnr_h':np.linspace(0.65,0.7,4)},
#       {'m_h':np.linspace(10.25,10.75,4)},
#        {'m_l':np.linspace(7.25,7.75,4)}],
#      [{'s_l':np.zeros(shape=4)+0.025},
#       {'ra_l':np.zeros(shape=4)+19.0},
#       {'nm_l':np.linspace(0.02,0.04,4)},
#       {'flnr_l':np.linspace(0.45,0.5,4)},
#        {'m_h':np.linspace(10.25,10.75,4)},
#        {'m_l':np.linspace(7.25,7.75,4)}]
#   ]        
        

#---------------Define Number of Variable Parameters---------------#  

num_params=5

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

        if 's_l' in params.keys():
            params['s']=params.pop('s_l')
        if 's_h' in params.keys():
            params['s']=params.pop('s_h') 
        
        if 'ra_l' in params.keys():
            params['ra']=params.pop('ra_l')
        if 'ra_h' in params.keys():
            params['ra']=params.pop('ra_h')
        
        if 'nm_l' in params.keys():
            params['nm']=params.pop('nm_l')
        if 'nm_h' in params.keys():
            params['nm']=params.pop('nm_h')
    
        if 'flnr_l' in params.keys():
            params['flnr']=params.pop('flnr_l')
        if 'flnr_h' in params.keys():
            params['flnr']=params.pop('flnr_h')
        
        if 'm_l1' in params.keys():
            params['m']=params.pop('m_l1')
        if 'm_h1' in params.keys():
            params['m']=params.pop('m_h1')
        
        if 'm_l2' in params.keys():
            params['m']=params.pop('m_l2')
        if 'm_h2' in params.keys():
            params['m']=params.pop('m_h2')
        
        if 'm_l3' in params.keys():
            params['m']=params.pop('m_l3')
        if 'm_h3' in params.keys():
            params['m']=params.pop('m_h3')            
        
        params_int+=[params]       
    
    leaf_params+=params_int