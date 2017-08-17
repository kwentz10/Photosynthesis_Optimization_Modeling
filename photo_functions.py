#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 13:02:35 2016

@author: Katherine

Functions for Photosynthesis Model
"""

import numpy as np

from scipy.interpolate import interp1d
from scipy.signal import savgol_filter
import math

from mpl_toolkits.mplot3d import Axes3D


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


def round_up_to_odd(f):
    return np.ceil(f) // 2 * 2 + 1


#Smoothing Function
def gs_smooth(yy,start,end,poly=3):
    
    i=0
    mylist=[]
    start_list=np.linspace(start,end,(end-start+1))
    
    for x in start_list:
        mylist+=[int(x)]
    
    for x in mylist:
        if yy[x]==0 or yy[x+1]==0 or yy[x+2]==0 or yy[x+3]==0 or yy[x+4]==0 or yy[x+5]==0:
            i+=1
        else:
            break
    
    start=start+i

    x=np.linspace(start,end,(end-start+1))
    y=yy[start:(end+1)]
    
    for i in range(len(y)):
        if math.isnan(y[i]):
            y[i]=np.mean([y[i-1],y[i+1]])
    
    xx = np.linspace(x.min(),x.max(),len(yy[start:end+1]))
    
    # interpolate + smooth
    itp = interp1d(x,y, kind='linear')
    window_size, poly_order = round_up_to_odd((end-start+1)/2.0), poly
    yy_sg = savgol_filter(itp(xx), window_size, poly_order)
    



    return np.zeros(len(yy[0:start])).tolist()+yy_sg.tolist()+np.zeros(len(yy[(end+1):])).tolist()
    



##Smoothing Function
#def gs_smooth_dy(yy):
#    x=np.linspace(122,264,264-122+1)
#    y=yy[122:265]
#    
#    for i in range(len(y)):
#        if math.isnan(y[i]):
#            y[i]=np.mean([y[i-1],y[i+1]])
#    
#    xx = np.linspace(x.min(),x.max(), 143)
#    
#    # interpolate + smooth
#    itp = interp1d(x,y, kind='linear')
#    window_size, poly_order = 59, 2
#    yy_sg = savgol_filter(itp(xx), window_size, poly_order)
#
#
#    return yy[0:122]+yy_sg.tolist()+yy[265:]
#
##Smoothing Function
#def gs_smooth_all(yy,poly=2):
#    x=np.linspace(1,365,365)
#    y=yy
#    
#    for i in range(len(y)):
#        if math.isnan(y[i]):
#            y[i]=np.mean([y[i-1],y[i+1]])
#    
#    xx = np.linspace(x.min(),x.max(), 365)
#    
#    # interpolate + smooth
#    itp = interp1d(x,y, kind='linear')
#    window_size, poly_order = 101, poly
#    yy_sg = savgol_filter(itp(xx), window_size, poly_order)
#
#
#    return yy_sg.tolist()
#
