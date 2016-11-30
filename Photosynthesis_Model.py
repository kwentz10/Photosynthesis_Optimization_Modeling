#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 26 16:37:08 2016

@author: Katherine
"""
#The line of code below imports the temperature functions used in this model
from photo_functions import arr_temp, bol_temp, pa_con_atmfrac
import numpy as np

def photo(s,nm,tl,ea,chl,crc,rub_max,ij,vwc,kc25,ko25,o,tau25,ca,rh,m,a,frnr,flnr,ra,j_b,j_m_max,q,vwc_min,vwc_max,b):

    ##---Calculated Parameter Arrays for Model(Constant+Variable Plant Trait(s))---##
    es_str=pa_con_atmfrac(611*np.exp(17.27*(tl-273.15)/((tl-273.15)+273.3))) #calculate saturation vapor pressure of surface (Pa)
    d=es_str-ea #calculate vapor pressure deficit (umol H2O/mol air)
    
    l=1/s #leaf mass per unit area (g C/m2 C)
    na=nm*l #leaf nitrogen (g N/ m2 C)
    
    #below is commented out because I am no longer using a variable lambda parameter
    #m=ca/(rh*d*lamb) ##Ball-Berry stomatal conductance slope parameter (unitless)

    rub=(chl*crc)/1000 # find ribulose bisphosphate content (umol RuBP/m2)
        
    if all(rub<rub_max):
        j_m=j_m_max*(rub/rub_max)*ij #find j_m slope based on ribulose bisphosphate content & leaf area/angle index
    else:
        j_m=j_m_max*ij
            
    vopt=frnr*flnr*ra*na #optimal carboxylation rate, limited by CO2 (umol CO2/m2s)
    jopt=vopt*j_m+j_b #optimal carboxylation rate, limited by RuBP (umol CO2/m2s)

    ##---Temperature Effects on Parameters---##
    
    #parameters
    tk_25=298.16; #absolute temperature at 25 C
    ekc=80500.0 #Activation energy for K of CO2 (J mol-1)
    eko=14500.0 #Activation energy for K of O2 (J mol-1)
    etau=-29000.0  #Activation energy for tau (???) (J mol-1)
    ev=55000.0 #Activation energy for carboxylation (J mol-1)
    ej=55000.0 #Activation energy for electron transport (J mol-1)
    toptv=298.0 #Optimum temperature for maximum carboxylation (K)
    toptj=298.0 #Optimum temperature for maximum electron transport (K)
        
    #calculated parameters due to temperature
    kc=arr_temp(kc25,ekc,tk_25,tl) #Michaelis-Menten kinetic coefficient for carbon dioxide at leaf temperature (umol/mol)
    ko=arr_temp(ko25,eko,tk_25,tl) #Michaelis-Menten kinetic coefficient for oxygen at leaf temperature (umol/mol) 
    tau=arr_temp(tau25,etau,tk_25,tl) #specifity coefficient of tau at leaf temperature (unitless) 
    gamma=o/(2*tau) #carbon dioxide compensation point (umol/mol)
    vmax1=bol_temp(vopt,ev,toptv,tl) #carboxylation rate at leaf temperature, limited by CO2 (umol CO2/m2s)
    jmax1=bol_temp(jopt,ej,toptj,tl) #carboxylation rate at leaf temperature, limited by RuBP (umol CO2/m2s)

    ##---Soil Moisture Effect on Parameters---##
    
    #below I removed the vwc constraint on photosynthesis because it is not a leaf trait
    Wfac=1
#    if all(vwc>=vwc_max):
#        Wfac=1
#    elif all(vwc<vwc_max):
#        Wfac=((vwc-vwc_min)/(vwc_max-vwc_min))**q
    
    vmax=Wfac*vmax1
    jmax=Wfac*jmax1
    
    ##---Define a1 and a2 depending on whether plant is rubisco limited or light limited---##

    #rubisco limited
    a1_r=vmax
    a2_r=kc*(1+(o/ko))
    #light limited
    a1_l=jmax/4
    a2_l=2*gamma

        
    ##---(1)Photosynthesis and Stomatal Conductance Models (b is not taken into account)---##

    if any(b==0.0):
    
        #In order to generate this model I combined the following equations:
        #A=gsc*(ca-ci)
        #gsc=gsw/a
        #gsw=mArh/ca
        #solve for ci
        #plug into A=a1(ci-gamma)/ci+a2
        #Rubisco Limiting: a1=vcmax; a2=kc(1+o/ko)
        #Light Limiting: a1=2.2*vcmax/4; a2=2*gamma

        #Solve for Assimilation
        
        ci=ca-((a*ca)/(m*rh)) #internal carbon dioxide (umol/mol)
        
        ##---Rubisco-Limited Assimilation---##
        A_r=(a1_r*(ci-gamma))/(ci+a2_r) #rubisco limited assimilation rate (umol CO2/m2s)
        
        ##---Light-Limited Assimilation---##
        A_l=(a1_l*(ci-gamma))/(ci+a2_l) #light limited assimilation rate (umol CO2/m2s)
        
        ##---Determine Rate-Limiting Assimilation---##
        A=[]
        for xx in range(len(A_r)):
            if A_r[xx]<A_l[xx]:
                A+=[A_r[xx]] #rubisco limited
            elif A_l[xx]<A_r[xx]:
                A+=[A_l[xx]] #light limited
            else: 
                A+=[A_l[xx]] #both light and rubisco limited
                
        ##---Solve for Stomatal Conductance to Water---##
        gsw=m*A*rh/ca #stomatal conductance to water (mol air/m2s)
        
        ##---Solve for Evapotranspiration---##
        E=gsw*d #(umol H2O/m2s)
    

    ##---(2)Photosynthesis and Stomatal Conductance Models (with b)---##

    elif any(b>0.0):
    
        #In order to generate this model I combined the following equations:
        #A=gsc*(ca-ci)
        #gsc=gsw/a
        #gsw=mArh/ca+b
        #solve for ci
        #plug into A=a1(ci-gamma)/ci+a2
        #Rubisco Limiting: a1=vcmax; a2=kc(1+o/ko)
        #Light Limiting: a1=2.2*vcmax/4; a2=2*gamma

        #Solve for Assimilation Using Quadratic Equation
        
        ##---Rubisco-Limited Assimilation---##
        aa_r=m*rh*ca-a*ca+m*rh*a2_r
        bb_r=b*(ca**2)+b*ca*a2_r-a1_r*m*rh*ca+a*ca*a1_r+a1_r*m*rh*gamma
        cc_r=a1_r*b*(ca**2)+gamma*b*ca*a1_r

        A1_r=(-bb_r+np.sqrt(bb_r**2-4*aa_r*cc_r))/(2*aa_r)
        A2_r=(-bb_r-np.sqrt(bb_r**2-4*aa_r*cc_r))/(2*aa_r)
                    
        #Choose Highest Values for Assimilation and Conductance
        A_r=[]
        for j in range(len(A1_r)):
            if A1_r[j]>A2_r[j]:
                A_r+=[A1_r[j]]
            elif A2_r[j]>A1_r[j]:
                A_r+=[A2_r[j]]
            else:
                A_r+=[A1_r[j]]

        ##---Light-Limited Assimilation---##
        aa_l=m*rh*ca-a*ca+m*rh*a2_l
        bb_l=b*(ca**2)+b*ca*a2_l-a1_l*m*rh*ca+a*ca*a1_l+a1_l*m*rh*gamma
        cc_l=a1_l*b*(ca**2)+gamma*b*ca*a1_l

        A1_l=(-bb_l+np.sqrt(bb_l**2-4*aa_l*cc_l))/(2*aa_l)
        A2_l=(-bb_l-np.sqrt(bb_l**2-4*aa_l*cc_l))/(2*aa_l)
            
        #Choose Highest Values for Assimilation and Conductance
        A_l=[]
        for j in range(len(A1_l)):
            if A1_l[j]>A2_l[j]:
                A_l+=[A1_l[j]]
            elif A2_l[j]>A1_l[j]:
                A_l+=[A2_l[j]]
            else:
                A_l+=[A1_l[j]]
    
        ##---Determine Rate-Limiting Assimilation---##
        A=[]
        for xx in range(len(A_r)):
            if A_r[xx]<A_l[xx]:
                A+=[A_r[xx]] #rubisco limited
            elif A_l[xx]<A_r[xx]:
                A+=[A_l[xx]] #light limited
            else: 
                A+=[A_l[xx]] #both light and rubisco limited         
        
        ##---Solve for Stomatal Conductance to Water---##
        gsw=m*A*rh/ca #stomatal conductance to water (mol H2O/m2s) #make array from list
            
        ##---Solve for Evapotranspiration---##
        E=gsw*d #(umol H2O/m2s)

        
        #---------------Test for Nan or Negative Values---------------#       
            
    for xxx in range(len(A)):
        if np.isnan(A[xxx]):
            print "A array contains nan values"
            return -999,-999
        if A[xxx]<0.0:
            print "A array contains negative values"
            return -999,-999
        if np.isnan(gsw[xxx]):
            print "gsw array contains nan values"
            return -999,-999
        if gsw[xxx]<0.0:
            print "gsw array contains negative values"
            return -999,-999
        if np.isnan(E[xxx]):
            print "E array contains nan values"
            return -999,-999
        if E[xxx]<0.0:
            print "E array contains negative values"
            return -999,-999

        
    #---------------WUE vs. NUE---------------#    
    
    
    wue=np.diff(A)/np.diff(E)*1000.0 #multiply by 1000 to get from umol CO2/umol H20 to umol CO2/mmol H20
    nue=np.diff(A)/np.diff(na)
    
    return wue, nue