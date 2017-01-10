#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 13:02:35 2016

@author: Katherine

Functions for Photosynthesis Model
"""

import numpy as np

#convert microbars to umol/mol
def bar_con_atmfrac(bars):
    
    atm=bars*0.987
    umolmol=atm*1000000
    
    return umolmol

#convert pascals to umol/mol
def pa_con_atmfrac(pa,atm_air):
    
    atm=pa*0.00000987 #this relationship holds for all altitudes
    umolmol=(atm/atm_air)*1000000 #altered this from atm*1000000 to (atm/atm_air)*1000000
    
    return umolmol
    
#Arhennius temperature function
def arr_temp(rate,ea,tref,tleaf):
    
    r=8.3144 #gas constant
    new_rate=rate*np.exp(((tleaf-tref)*ea)/(tref*r*tleaf))
    
    return new_rate
    
#Boltzmann temperature distrubtion function
def bol_temp(rate,ea,topt,tleaf):
    
    r=8.3144 #gas constant
    h=200000.0 #enthalpy term
        
    dtopt=tleaf-topt
    prodt=r*topt*tleaf
    
    numm=rate*h*np.exp(ea*dtopt/prodt)
    denom=h-ea*(1.0-np.exp(h*dtopt/prodt))
    new_rate=numm/denom
    
    return new_rate
    