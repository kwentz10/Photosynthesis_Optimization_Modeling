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
def pa_con_atmfrac(pa_v,z):
    Tb= 288.15 #standard temperature (K)
    Pb=101325.00 #static pressure (Pa)
    Lb=-0.0065 #standard temperature lapse rate (K/m) in ISA
    h=z
    hb=0.0 #height at bottom of layer b (m)
    g0=9.80665 #gravitational acceleration (m/s2)
    M=0.0289644 #molar mass of Earth's air (kg/mol)
    R=8.3144598 #universal gas constant (J/molK)
    
    #atm_v=pa_v*0.00000987 #this relationship holds for all altitudes
    pa_air=(Pb*(Tb/(Tb+Lb*(h-hb)))**((g0*M)/(R*Lb))) #Pa (look for equation)
    #atm_air=pa_air*0.00000987
    umolmol=(pa_v/pa_air)*1000000 #altered this from atm*1000000 to (atm/atm_air)*1000000
    
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
    