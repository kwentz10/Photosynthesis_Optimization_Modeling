#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 26 16:37:08 2016

In this code I solve for assimilation by three approaches. All
three approaches use a combination of Ohm's Law and Ball-Berry's 
Model of stomatal conductances. 

The first photo model solves for Assimilation without factoring in boundary 
layer or mesophyll conductance.

The second photo model solves for Assimilation using a cubic equation which 
only includes boundary layer conductance.

The third photo model solves for Assimilation using a cubic equation which
includes boundary layer and mesophyll conductance. 

The fourth photo model solves for Assimilation using a cubic equation which
includes boundary layer and mesophyll conductance. In this model I set 
mesophyll conductance equal to stomatal conductance.

See Supplementary Eq. Worksheet for more information.


@author: Katherine
"""

#############------Import Modules------#############

from photo_functions import arr_temp, bol_temp, pa_con_atmfrac
import numpy as np
import math


#############------Photosynthesis Model 1------#############

def photo_no_bound_meso(tk_25,ekc,eko,etau,ev,ej,toptv,toptj,na,qeff, PAR, tl,ea,chl,ij,kc25,ko25,o,ca,rh,m,a,frnr,flnr,ra,j_m,g0,q,vwc,vwc_max,vwc_min):
    
    if all(ij>1.0):
        ij=np.zeros(shape=2)+1.0
    else:
        ij=ij
    
    ##---Calculated Parameter Arrays for Model(Constant+Variable Plant Trait(s))---##
    pa_v=611*np.exp((17.27*tl)/(tl+237.3)) #calculate saturation vapor pressure of surface (Pa)
    es_str=pa_con_atmfrac(pa_v,3528) #calculate saturation vapor pressure of surface (Pa-->umol/mol)
    d=es_str-ea #calculate vapor pressure deficit (umol H2O/mol air)
    
             
    vmaxopt=frnr*flnr*ra*na #optimal carboxylation rate, limited by CO2 (umol CO2/m2s)
    jmaxopt=vmaxopt*j_m #optimal carboxylation rate, limited by RuBP (umol CO2/m2s)
    
    ##---Temperature Effects on Parameters---##
    if tl<=0.0:
        nue=[0.0]
        wue=[0.0]
        A=[0.0]
        E=[0.0]
        cs=[0.0]
        ci=[0.0]
        gsw=[0.0]
        gs=[0.0]
        gbw=[0.0]
        gb=[0.0]
        gm=[0.0]
        cc=[0.0]
        dd=[0.0]
        return wue, nue, A, E, cs, ci, gsw, gs, gbw,gb,gm,cc,dd
    
        
    #calculated parameters due to temperature
    kc=arr_temp(pa_con_atmfrac(kc25,3528),ekc,tk_25,tl+273.15) #Michaelis-Menten kinetic coefficient for carbon dioxide at leaf temperature (umol/mol)
    ko=arr_temp(pa_con_atmfrac(ko25,3528),eko,tk_25,tl+273.15) #Michaelis-Menten kinetic coefficient for oxygen at leaf temperature (umol/mol) 
    gamma=(kc*o)/(2.0*ko)*0.21 #carbon dioxide compensation point (umol/mol)
    vmax1=bol_temp(vmaxopt,ev,toptv,tl+273.15) #carboxylation rate at leaf temperature, limited by CO2 (umol CO2/m2s)
    jmax1=bol_temp(jmaxopt,ej,toptj,tl+273.15) #carboxylation rate at leaf temperature, limited by RuBP (umol CO2/m2s)

    ##---Soil Moisture Effect on Parameters---##

   #below I removed the vwc constraint on photosynthesis because it is not a leaf trait

    if math.isnan(vwc):
        nue=[np.nan]
        wue=[np.nan]
        A=[np.nan]
        E=[np.nan]
        cs=[np.nan]
        ci=[np.nan]
        gsw=[np.nan]
        gs=[np.nan]
        gbw=[np.nan]
        gb=[np.nan]
        gm=[np.nan]
        cc=[np.nan]
        dd=[np.nan]
        return wue, nue, A, E, cs, ci, gsw, gs, gbw,gb,gm,cc,dd   
    
    if vwc<=vwc_min:
        nue=[0.0]
        wue=[0.0]
        A=[0.0]
        E=[0.0]
        cs=[0.0]
        ci=[0.0]
        gsw=[0.0]
        gs=[0.0]
        gbw=[0.0]
        gb=[0.0]
        gm=[0.0]
        cc=[0.0]
        dd=[0.0]
        return wue, nue, A, E, cs, ci, gsw, gs, gbw,gb,gm,cc,dd
   
    
    #below I removed the vwc constraint on photosynthesis because it is not a leaf trait
    Wfac=1.0
    if all(vwc>=vwc_max):
        Wfac=1.0
    elif all(vwc<vwc_max):
        Wfac=((vwc-vwc_min)/(vwc_max-vwc_min))**q
    

    
    vmax=Wfac*vmax1
    jmax=Wfac*jmax1

    if vmax<=0:
        nue=[0.0]
        wue=[0.0]
        A=[0.0]
        E=[0.0]
        cs=[0.0]
        ci=[0.0]
        gsw=[0.0]
        gs=[0.0]
        gbw=[0.0]
        gb=[0.0]
        gm=[0.0]
        cc=[0.0]
        dd=[0.0]
        return wue, nue, A, E, cs, ci, gsw, gs, gbw,gb,gm,cc,dd
    
    ##---Determine J---## 
    alpha=(chl/1000.)/((chl/1000.)+0.076) #from Developmental Constratins on Photosynthesis: Effects of Light and Nutrition (Evans)
    qalpha=alpha*qeff
    Iphoton=PAR*ij
    
    
    if all(jmax>0.0):
        jj=(qalpha*Iphoton)/(np.sqrt(1.+(((qalpha**2)*(Iphoton**2))/(jmax**2))))
    
    
    ##---Define a1 and a2 depending on whether plant is rubisco limited or light limited---##

    #rubisco limited
    a1_r=vmax
    a2_r=kc*(1+(o/ko))
    #light limited
    a1_l=jj/4
    a2_l=2*gamma
    
        
    ##---(1)Photosynthesis and Stomatal Conductance Models (b is not taken into account)---##

    if all(g0==0.0):
    
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
        
        ##---Solve for Stomatal Conductance to CO2---##
        gs=gsw*a
        
        ##---Solve for Evapotranspiration---##
        E=gsw*d #(umol H2O/m2s)
    

    ##---(2)Photosynthesis and Stomatal Conductance Models (with b)---##

    elif any(g0>0.0):
    
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
        bb_r=g0*(ca**2)+g0*ca*a2_r-a1_r*m*rh*ca+a*ca*a1_r+a1_r*m*rh*gamma
        cc_r=a1_r*g0*(ca**2)+gamma*g0*ca*a1_r

        A1_r=(-bb_r+np.sqrt(bb_r**2-4*aa_r*cc_r))/(2*aa_r) #(umol CO2/m2s)
        A2_r=(-bb_r-np.sqrt(bb_r**2-4*aa_r*cc_r))/(2*aa_r) #(umol CO2/m2s)
                    
        #Choose Highest Values for Assimilation and Conductance
        A_r=[]
        
        for j in range(len(A1_r)):
            
            if A1_r[j]<A2_r[j]:
                if A1_r[j]>0.0:
                    A_r+=[A1_r[j]]
                else:
                    A_r+=[A2_r[j]]
            
            elif A1_r[j]>A2_r[j]:
                if A2_r[j]>0.0:
                    A_r+=[A2_r[j]]
                else:
                    A_r+=[A1_r[j]]
            
            elif A1_r[j]==A2_r[j]:
                A_r+=[A2_r[j]]


        ##---Light-Limited Assimilation---##
        aa_l=m*rh*ca-a*ca+m*rh*a2_l
        bb_l=g0*(ca**2)+g0*ca*a2_l-a1_l*m*rh*ca+a*ca*a1_l+a1_l*m*rh*gamma
        cc_l=a1_l*g0*(ca**2)+gamma*g0*ca*a1_l

        A1_l=(-bb_l+np.sqrt(bb_l**2-4*aa_l*cc_l))/(2*aa_l) #(umol CO2/m2s)
        A2_l=(-bb_l-np.sqrt(bb_l**2-4*aa_l*cc_l))/(2*aa_l) #(umol CO2/m2s)
            
        #Choose Highest Values for Assimilation and Conductance (do minimum that is not negative!!)
        A_l=[]
        
        for j in range(len(A1_l)):
            
            if A1_l[j]<A2_l[j]:
                if A1_l[j]>0.0:
                    A_l+=[A1_l[j]]
                else:
                    A_l+=[A2_l[j]]
            
            elif A1_l[j]>A2_l[j]:
                if A2_l[j]>0.0:
                    A_l+=[A2_l[j]]
                else:
                    A_l+=[A1_l[j]]
            
            elif A1_l[j]==A2_l[j]:
                A_l+=[A2_l[j]]
    
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
        gsw=(m*A*rh/ca)+g0 #stomatal conductance to water (mol H2O/m2s) #make array from list
         
        ##---Solve for Stomatal Conductance to CO2---##
        gs=gsw*a
   
        ##---Solve for Evapotranspiration---##
        E=gsw*d #(umol H2O/m2s)
        
        
        
        #---------------Test for Nan or Negative Values---------------#       
  
    neg_vals=([np.zeros(1)+0]*13)
    
    for xxx in range(len(A)):
        if np.isnan(A[xxx]):
            print "A array contains nan values"
            return neg_vals[0],neg_vals[1],neg_vals[2],neg_vals[3],neg_vals[4],neg_vals[5],neg_vals[6],neg_vals[7],neg_vals[8],neg_vals[9],neg_vals[10],neg_vals[11],neg_vals[12]
        if A[xxx]<0.0:
            print "A array contains negative values"
            return neg_vals[0],neg_vals[1],neg_vals[2],neg_vals[3],neg_vals[4],neg_vals[5],neg_vals[6],neg_vals[7],neg_vals[8],neg_vals[9],neg_vals[10],neg_vals[11],neg_vals[12]
        if np.isnan(gsw[xxx]):
            print "gsw array contains nan values"
            return neg_vals[0],neg_vals[1],neg_vals[2],neg_vals[3],neg_vals[4],neg_vals[5],neg_vals[6],neg_vals[7],neg_vals[8],neg_vals[9],neg_vals[10],neg_vals[11],neg_vals[12]
        if gsw[xxx]<0.0:
            print "gsw array contains negative values"
            return neg_vals[0],neg_vals[1],neg_vals[2],neg_vals[3],neg_vals[4],neg_vals[5],neg_vals[6],neg_vals[7],neg_vals[8],neg_vals[9],neg_vals[10],neg_vals[11],neg_vals[12]
        if np.isnan(E[xxx]):
            print "E array contains nan values"
            return neg_vals[0],neg_vals[1],neg_vals[2],neg_vals[3],neg_vals[4],neg_vals[5],neg_vals[6],neg_vals[7],neg_vals[8],neg_vals[9],neg_vals[10],neg_vals[11],neg_vals[12]
        if E[xxx]<0.0:
            print "E array contains negative values"
            return neg_vals[0],neg_vals[1],neg_vals[2],neg_vals[3],neg_vals[4],neg_vals[5],neg_vals[6],neg_vals[7],neg_vals[8],neg_vals[9],neg_vals[10],neg_vals[11],neg_vals[12]

        
    #---------------WUE vs. NUE---------------#    

    wue=A/E*1000.0 #multiply by 1000 to get from umol CO2/umol H20 to umol CO2/mmol H20
    nue=A/na   
    
    
    return wue, nue, A, gsw, E, na, gs
  
    
    
    
    

    
#############------Photosynthesis Model 2------#############

def photo_bound(tk_25,ekc,eko,etau,ev,ej,toptv,toptj,na, qeff, PAR, tl,ea,chl,ij,kc25,ko25,o,ca,rh,m,a,frnr,flnr,ra,j_m,g0,b,dia,u,q,vwc,vwc_max,vwc_min):  
    
    if all(ij>1.0):
        ij=np.zeros(shape=2)+1.0
    else:
        ij=ij
    
    ##---Calculated Parameter Arrays for Model(Constant+Variable Plant Trait(s))---##
    pa_v=611*np.exp((17.27*tl)/(tl+237.3)) #calculate saturation vapor pressure of surface (Pa)
    es_str=pa_con_atmfrac(pa_v,3528) #calculate saturation vapor pressure of surface (Pa-->umol/mol)
    d=es_str-ea #calculate vapor pressure deficit (umol H2O/mol air)
    
          
    vmaxopt=frnr*flnr*ra*na #optimal carboxylation rate, limited by CO2 (umol CO2/m2s)
    jmaxopt=vmaxopt*j_m #optimal carboxylation rate, limited by RuBP (umol CO2/m2s)
    
    ##---Temperature Effects on Parameters---##
    if tl<=0.0:
        nue=[0.0]
        wue=[0.0]
        A=[0.0]
        E=[0.0]
        cs=[0.0]
        ci=[0.0]
        gsw=[0.0]
        gs=[0.0]
        gbw=[0.0]
        gb=[0.0]
        gm=[0.0]
        cc=[0.0]
        dd=[0.0]
        return wue, nue, A, E, cs, ci, gsw, gs, gbw,gb,gm,cc,dd
        
    #calculated parameters due to temperature
    kc=arr_temp(pa_con_atmfrac(kc25,3528),ekc,tk_25,tl+273.15) #Michaelis-Menten kinetic coefficient for carbon dioxide at leaf temperature (umol/mol)
    ko=arr_temp(pa_con_atmfrac(ko25,3528),eko,tk_25,tl+273.15) #Michaelis-Menten kinetic coefficient for oxygen at leaf temperature (umol/mol) 
    gamma=(kc*o)/(2.0*ko)*0.21 #carbon dioxide compensation point (umol/mol)
    vmax1=bol_temp(vmaxopt,ev,toptv,tl+273.15) #carboxylation rate at leaf temperature, limited by CO2 (umol CO2/m2s)
    jmax1=bol_temp(jmaxopt,ej,toptj,tl+273.15) #carboxylation rate at leaf temperature, limited by RuBP (umol CO2/m2s)

    ##---Soil Moisture Effect on Parameters---##

   #below I removed the vwc constraint on photosynthesis because it is not a leaf trait

    if math.isnan(vwc):
        nue=[np.nan]
        wue=[np.nan]
        A=[np.nan]
        E=[np.nan]
        cs=[np.nan]
        ci=[np.nan]
        gsw=[np.nan]
        gs=[np.nan]
        gbw=[np.nan]
        gb=[np.nan]
        gm=[np.nan]
        cc=[np.nan]
        dd=[np.nan]
        return wue, nue, A, E, cs, ci, gsw, gs, gbw,gb,gm,cc,dd   
    
    if vwc<=vwc_min:
        nue=[0.0]
        wue=[0.0]
        A=[0.0]
        E=[0.0]
        cs=[0.0]
        ci=[0.0]
        gsw=[0.0]
        gs=[0.0]
        gbw=[0.0]
        gb=[0.0]
        gm=[0.0]
        cc=[0.0]
        dd=[0.0]
        return wue, nue, A, E, cs, ci, gsw, gs, gbw,gb,gm,cc,dd
   
    
    #below I removed the vwc constraint on photosynthesis because it is not a leaf trait
    Wfac=1.0
    if all(vwc>=vwc_max):
        Wfac=1.0
    elif all(vwc<vwc_max):
        Wfac=((vwc-vwc_min)/(vwc_max-vwc_min))**q
    
    vmax=Wfac*vmax1
    jmax=Wfac*jmax1
    
    if vmax<=0:
        nue=[0.0]
        wue=[0.0]
        A=[0.0]
        E=[0.0]
        cs=[0.0]
        ci=[0.0]
        gsw=[0.0]
        gs=[0.0]
        gbw=[0.0]
        gb=[0.0]
        gm=[0.0]
        cc=[0.0]
        dd=[0.0]
        return wue, nue, A, E, cs, ci, gsw, gs, gbw,gb,gm,cc,dd
    
    ##---Determine J---## 
    alpha=(chl/1000.)/((chl/1000.)+0.076) #from Developmental Constratins on Photosynthesis: Effects of Light and Nutrition (Evans)
    qalpha=alpha*qeff
    Iphoton=PAR*ij
    
    if all(jmax>0.0):
        jj=(qalpha*Iphoton)/(np.sqrt(1.+(((qalpha**2)*(Iphoton**2))/(jmax**2))))
    
    
    ##---Define a1 and a2 depending on whether plant is rubisco limited or light limited---##

    #rubisco limited
    a1_r=vmax
    a2_r=kc*(1+(o/ko))
    #light limited
    a1_l=jj/4
    a2_l=2*gamma
    
    
     ##---Photosynthesis Model---##
    

    #Solve for Assimilation Using Cubic Equation
    
    #constants
    gbw=27.0/(200*np.sqrt(dia/u))
    gb=b*gbw
    
    C1=(-a*m*rh*gb)+(a*g0)+gb
    C2=(ca*(gb**2)*a*m*rh)-(ca*gb*a*g0)-(a*g0*ca*gb)-(ca*(gb**2))
    C3=(ca**2)*(gb**2)*a*g0
    C4=(a*m*rh*(gb**2))-(a*g0*gb)
    C5=a*g0*ca*(gb**2)
    
    ##---Rubisco-Limited Assimilation---##
    
    #make list of coefficients for cubic eq.
    coef=[]
    for i in range(len(C1)):
        coef+=[[C1[i],C2[i]+C4[i]*a2_r[i]-a1_r[i]*C1[i],C3[i]+C5[i]*a2_r[i]-a1_r[i]*C2[i]+C4[i]*a1_r[i]*gamma[i],-a1_r[i]*C3[i]+C5[i]*a1_r[i]*gamma[i]]]
   

    #Solve Cubic Eq. for Roots
    roots=[]
    for i in range(len(coef)):
        roots+=[np.roots(coef[i])]
    
    #make list of real roots
    roots_real=[]
    for i in range(len(roots)):
        roots_real+=[roots[i].real[abs(roots[i].imag)<1e-5]]
    
    #make list of roots greater than zero
    roots_pos=[]
    for i in range(len(roots_real)):  
        roots_pos_0=[]
        for ii in range(len(roots_real[i])):
            if roots_real[i][ii]>0.0:
                roots_pos_0+=[roots_real[i][ii]]
            else:
                continue
        roots_pos+=[roots_pos_0]
  
    #assimilation is the minimum of the roots that are greater than zero
    A_r=[]
    for i in range(len(roots_pos)):
        A_r+=[np.min(roots_pos[i])]


    ##---Light-Limited Assimilation---##
    #make list of coefficients for cubic eq.
    coef=[]
    for i in range(len(C1)):
        coef+=[[C1[i],C2[i]+C4[i]*a2_l[i]-a1_l[i]*C1[i],C3[i]+C5[i]*a2_l[i]-a1_l[i]*C2[i]+C4[i]*a1_l[i]*gamma[i],-a1_l[i]*C3[i]+C5[i]*a1_l[i]*gamma[i]]]
    
    #Solve Cubic Eq. for Roots
    roots=[]
    for i in range(len(coef)):
        roots+=[np.roots(coef[i])]
    
    #make list of real roots
    roots_real=[]
    for i in range(len(roots)):
        roots_real+=[roots[i].real[abs(roots[i].imag)<1e-5]] #or do: roots[np.isreal(roots)]
    
    #make list of roots greater than zero
    roots_pos=[]
    for i in range(len(roots_real)):  
        roots_pos_0=[]
        for ii in range(len(roots_real[i])):
            if roots_real[i][ii]>0.0:
                roots_pos_0+=[roots_real[i][ii]]
            else:
                continue
        roots_pos+=[roots_pos_0]
        
    #assimilation is the minimum of the roots that are greater than zero
    A_l=[]
    for i in range(len(roots_pos)):
        A_l+=[np.min(roots_pos[i])]
                         
    ##---Determine Rate-Limiting Assimilation---##
    A=[]
    for i in range(len(A_r)):
        if A_r[i]<A_l[i]:
            A+=[A_r[i]] #rubisco limited
        elif A_l[i]<A_r[i]:
            A+=[A_l[i]] #light limited
        else: 
            A+=[A_l[i]] #both light and rubisco limited         
    
    
    

    ##---Solve for Stomatal Conductance to Water---##
    cs=ca-(A/gb) #stomatal co2 (umol CO2/mol)
    gsw=(m*A*rh/cs)+g0 #stomatal conductance to water (mol H2O/m2s) #make array from list

    ##---Solve for Stomatal Conductance to CO2---##
    gs=gsw*a
    
    ##---Solve for Internal CO2 Concentration (in Mesophyll)---##
    ci=cs-(A/gs)
        
    ##---Solve for Evapotranspiration---##
    E=gsw*d #(umol H2O/m2s)

    
    #---------------Test for Nan or Negative Values---------------#       
  
    neg_vals=([np.zeros(1)+0]*13)
    
    for xxx in range(len(A)):
        if np.isnan(A[xxx]):
            print "A array contains nan values"
            return neg_vals[0],neg_vals[1],neg_vals[2],neg_vals[3],neg_vals[4],neg_vals[5],neg_vals[6],neg_vals[7],neg_vals[8],neg_vals[9],neg_vals[10],neg_vals[11],neg_vals[12]
        if A[xxx]<0.0:
            print "A array contains negative values"
            return neg_vals[0],neg_vals[1],neg_vals[2],neg_vals[3],neg_vals[4],neg_vals[5],neg_vals[6],neg_vals[7],neg_vals[8],neg_vals[9],neg_vals[10],neg_vals[11],neg_vals[12]
        if np.isnan(gsw[xxx]):
            print "gsw array contains nan values"
            return neg_vals[0],neg_vals[1],neg_vals[2],neg_vals[3],neg_vals[4],neg_vals[5],neg_vals[6],neg_vals[7],neg_vals[8],neg_vals[9],neg_vals[10],neg_vals[11],neg_vals[12]
        if gsw[xxx]<0.0:
            print "gsw array contains negative values"
            return neg_vals[0],neg_vals[1],neg_vals[2],neg_vals[3],neg_vals[4],neg_vals[5],neg_vals[6],neg_vals[7],neg_vals[8],neg_vals[9],neg_vals[10],neg_vals[11],neg_vals[12]
        if np.isnan(E[xxx]):
            print "E array contains nan values"
            return neg_vals[0],neg_vals[1],neg_vals[2],neg_vals[3],neg_vals[4],neg_vals[5],neg_vals[6],neg_vals[7],neg_vals[8],neg_vals[9],neg_vals[10],neg_vals[11],neg_vals[12]
        if E[xxx]<0.0:
            print "E array contains negative values"
            return neg_vals[0],neg_vals[1],neg_vals[2],neg_vals[3],neg_vals[4],neg_vals[5],neg_vals[6],neg_vals[7],neg_vals[8],neg_vals[9],neg_vals[10],neg_vals[11],neg_vals[12]

        
    #---------------WUE vs. NUE---------------#    
    
    wue=A/E*1000.0 #multiply by 1000 to get from umol CO2/umol H20 to umol CO2/mmol H20
    nue=A/na   
    
    return wue, nue, A, E, na, cs, ci, gsw, gs, gbw, gb        

    
    
    
    
    

    
#############------Photosynthesis Model 3------#############
    
def photo_bound_meso(tk_25,ekc,eko,etau,ev,ej,toptv,toptj,na, qeff, PAR, s,tl,ea,chl,ij,kc25,ko25,o,ca,rh,m,a,frnr,flnr,ra,j_m,g0,b,dia,u,gm,q,vwc_min,vwc_max,vwc):
  
    
    if all(ij>1.0):
        ij=np.zeros(shape=2)+1.0
    else:
        ij=ij
    
    ##---Calculated Parameter Arrays for Model(Constant+Variable Plant Trait(s))---##
    pa_v=611*np.exp((17.27*tl)/(tl+237.3)) #calculate saturation vapor pressure of surface (Pa)
    es_str=pa_con_atmfrac(pa_v,3528) #calculate saturation vapor pressure of surface (Pa-->umol/mol)
    d=es_str-ea #calculate vapor pressure deficit (umol H2O/mol air)

        
    vmaxopt=frnr*flnr*ra*na #optimal carboxylation rate, limited by CO2 (umol CO2/m2s)
    jmaxopt=vmaxopt*j_m #optimal carboxylation rate, limited by RuBP (umol CO2/m2s)
    
    ##---Temperature Effects on Parameters---##
    if tl<=0.0:
        nue=[0.0]
        wue=[0.0]
        A=[0.0]
        E=[0.0]
        cs=[0.0]
        ci=[0.0]
        gsw=[0.0]
        gs=[0.0]
        gbw=[0.0]
        gb=[0.0]
        gm=[0.0]
        cc=[0.0]
        dd=[0.0]
        return wue, nue, A, E, cs, ci, gsw, gs, gbw,gb,gm,cc,dd

        
    #calculated parameters due to temperature
    kc=arr_temp(pa_con_atmfrac(kc25,3528),ekc,tk_25,tl+273.15) #Michaelis-Menten kinetic coefficient for carbon dioxide at leaf temperature (umol/mol)
    ko=arr_temp(pa_con_atmfrac(ko25,3528),eko,tk_25,tl+273.15) #Michaelis-Menten kinetic coefficient for oxygen at leaf temperature (umol/mol) 
    gamma=(kc*o)/(2.0*ko)*0.21 #carbon dioxide compensation point (umol/mol)
    vmax1=bol_temp(vmaxopt,ev,toptv,tl+273.15) #carboxylation rate at leaf temperature, limited by CO2 (umol CO2/m2s)
    jmax1=bol_temp(jmaxopt,ej,toptj,tl+273.15) #carboxylation rate at leaf temperature, limited by RuBP (umol CO2/m2s)

    ##---Soil Moisture Effect on Parameters---##
    
   #below I removed the vwc constraint on photosynthesis because it is not a leaf trait

    if math.isnan(vwc):
        nue=[np.nan]
        wue=[np.nan]
        A=[np.nan]
        E=[np.nan]
        cs=[np.nan]
        ci=[np.nan]
        gsw=[np.nan]
        gs=[np.nan]
        gbw=[np.nan]
        gb=[np.nan]
        gm=[np.nan]
        cc=[np.nan]
        dd=[np.nan]
        return wue, nue, A, E, cs, ci, gsw, gs, gbw,gb,gm,cc,dd   
    
    if vwc<=vwc_min:
        nue=[0.0]
        wue=[0.0]
        A=[0.0]
        E=[0.0]
        cs=[0.0]
        ci=[0.0]
        gsw=[0.0]
        gs=[0.0]
        gbw=[0.0]
        gb=[0.0]
        gm=[0.0]
        cc=[0.0]
        dd=[0.0]
        return wue, nue, A, E, cs, ci, gsw, gs, gbw,gb,gm,cc,dd
   

    #below I removed the vwc constraint on photosynthesis because it is not a leaf trait
    Wfac=1.0
    if all(vwc>=vwc_max):
        Wfac=1.0
    elif all(vwc<vwc_max):
        Wfac=((vwc-vwc_min)/(vwc_max-vwc_min))**q
    
    vmax=Wfac*vmax1
    jmax=Wfac*jmax1

    if vmax<=0:
        nue=[0.0]
        wue=[0.0]
        A=[0.0]
        E=[0.0]
        cs=[0.0]
        ci=[0.0]
        gsw=[0.0]
        gs=[0.0]
        gbw=[0.0]
        gb=[0.0]
        gm=[0.0]
        cc=[0.0]
        dd=[0.0]
        return wue, nue, A, E, cs, ci, gsw, gs, gbw,gb,gm,cc,dd
    
    
    ##---Determine J---## 
    alpha=(chl/1000.)/((chl/1000.)+0.076) #from (Photosynthetic acclimation of plants to growth irradiance)Developmental Constratins on Photosynthesis: Effects of Light and Nutrition (Evans)
    qalpha=alpha*qeff
    Iphoton=PAR*ij
    
    if all(jmax>0.0):
        jj=(qalpha*Iphoton)/(np.sqrt(1.+(((qalpha**2)*(Iphoton**2))/(jmax**2)))) #see baldocchi lecture 10 notes for this eq.
    
    
    ##---Define a1 and a2 depending on whether plant is rubisco limited or light limited---##

    #rubisco limited
    a1_r=vmax
    a2_r=kc*(1+(o/ko))
    #light limited
    a1_l=jj/4
    a2_l=2*gamma
    
    
    ##---Photosynthesis Model---##
    

    #Solve for Assimilation Using Cubic Equation
    
    #constants
    gbw=27./(200*np.sqrt(dia/u))
    gb=b*gbw
    
    C1=(-a*m*rh*gb)+(a*g0)+gb
    C2=(ca*(gb**2)*a*m*rh)-(ca*gb*a*g0)-(a*g0*ca*gb)-(ca*(gb**2))
    C3=(ca**2)*(gb**2)*a*g0
    C4=(a*m*rh*(gb**2))-(a*g0*gb)
    C5=a*g0*ca*(gb**2)
    C6=C1*gm-C4
    C7=C2*gm-C5
    C8=C3*gm
    C9=C4*gm
    C10=C5*gm
    
    ##---Rubisco-Limited Assimilation---##
    
    #make list of coefficients for cubic eq.
    coef=[]
    for i in range(len(C1)):
        coef+=[[C6[i],C7[i]+C9[i]*a2_r[i]-a1_r[i]*C6[i],C8[i]+C10[i]*a2_r[i]-a1_r[i]*C7[i]+C9[i]*a1_r[i]*gamma[i],-a1_r[i]*C8[i]+C10[i]*a1_r[i]*gamma[i]]]
    
    #Solve Cubic Eq. for Roots
    roots=[]
    for i in range(len(coef)):
        roots+=[np.roots(coef[i])]
    
    #make list of real roots
    roots_real=[]
    for i in range(len(roots)):
        roots_real+=[roots[i].real[abs(roots[i].imag)<1e-5]]
    
    #make list of roots greater than zero
    roots_pos=[]
    for i in range(len(roots_real)):  
        roots_pos_0=[]
        for ii in range(len(roots_real[i])):
            if roots_real[i][ii]>0.0:
                roots_pos_0+=[roots_real[i][ii]]
            else:
                continue
        roots_pos+=[roots_pos_0]
    
    #assimilation is the minimum of the roots that are greater than zero
    A_r=[]
    for i in range(len(roots_pos)):
        A_r+=[np.min(roots_pos[i])]


    ##---Light-Limited Assimilation---##
    #make list of coefficients for cubic eq.
    coef=[]
    for i in range(len(C1)):
        coef+=[[C6[i],C7[i]+C9[i]*a2_l[i]-a1_l[i]*C6[i],C8[i]+C10[i]*a2_l[i]-a1_l[i]*C7[i]+C9[i]*a1_l[i]*gamma[i],-a1_l[i]*C8[i]+C10[i]*a1_l[i]*gamma[i]]]
    
    #Solve Cubic Eq. for Roots
    roots=[]
    for i in range(len(coef)):
        roots+=[np.roots(coef[i])]
    
    #make list of real roots
    roots_real=[]
    for i in range(len(roots)):
        roots_real+=[roots[i].real[abs(roots[i].imag)<1e-5]] #or do: roots[np.isreal(roots)]
    
    #make list of roots greater than zero
    roots_pos=[]
    for i in range(len(roots_real)):  
        roots_pos_0=[]
        for ii in range(len(roots_real[i])):
            if roots_real[i][ii]>0.0:
                roots_pos_0+=[roots_real[i][ii]]
            else:
                continue
        roots_pos+=[roots_pos_0]
    
    #assimilation is the minimum of the roots that are greater than zero
    A_l=[]
    for i in range(len(roots_pos)):
        A_l+=[np.min(roots_pos[i])]
                         
    ##---Determine Rate-Limiting Assimilation---##
    A=[]
    for i in range(len(A_r)):
        if A_r[i]<A_l[i]:
            A+=[A_r[i]] #rubisco limited
        elif A_l[i]<A_r[i]:
            A+=[A_l[i]] #light limited
        else: 
            A+=[A_l[i]] #both light and rubisco limited         
    
    
    

    ##---Solve for Stomatal Conductance to Water---##
    cs=ca-(A/gb) #stomatal co2 (umol CO2/mol)
    gsw=(m*A*rh/cs)+g0 #stomatal conductance to water (mol H2O/m2s) #make array from list

    ##---Solve for Stomatal Conductance to CO2---##
    gs=gsw*a
    
    ##---Solve for Internal CO2 Concentration (in Mesophyll)---##
    ci=cs-(A/gs)
    
    ##---Solve for Internal CO2 Concentration (in Chloroplast)---##    
    cc=ci-(A/gm)    
    
    ##---Solve for Evapotranspiration---##
    E=gsw*d #(umol H2O/m2s)

    
    #---------------Test for Nan or Negative Values---------------#       
  
    neg_vals=([np.zeros(1)+0]*13)
    
    for xxx in range(len(A)):
        if np.isnan(A[xxx]):
            print "A array contains nan values"
            return neg_vals[0],neg_vals[1],neg_vals[2],neg_vals[3],neg_vals[4],neg_vals[5],neg_vals[6],neg_vals[7],neg_vals[8],neg_vals[9],neg_vals[10],neg_vals[11],neg_vals[12]
        if A[xxx]<0.0:
            print "A array contains negative values"
            return neg_vals[0],neg_vals[1],neg_vals[2],neg_vals[3],neg_vals[4],neg_vals[5],neg_vals[6],neg_vals[7],neg_vals[8],neg_vals[9],neg_vals[10],neg_vals[11],neg_vals[12]
        if np.isnan(gsw[xxx]):
            print "gsw array contains nan values"
            return neg_vals[0],neg_vals[1],neg_vals[2],neg_vals[3],neg_vals[4],neg_vals[5],neg_vals[6],neg_vals[7],neg_vals[8],neg_vals[9],neg_vals[10],neg_vals[11],neg_vals[12]
        if gsw[xxx]<0.0:
            print "gsw array contains negative values"
            return neg_vals[0],neg_vals[1],neg_vals[2],neg_vals[3],neg_vals[4],neg_vals[5],neg_vals[6],neg_vals[7],neg_vals[8],neg_vals[9],neg_vals[10],neg_vals[11],neg_vals[12]
        if np.isnan(E[xxx]):
            print "E array contains nan values"
            return neg_vals[0],neg_vals[1],neg_vals[2],neg_vals[3],neg_vals[4],neg_vals[5],neg_vals[6],neg_vals[7],neg_vals[8],neg_vals[9],neg_vals[10],neg_vals[11],neg_vals[12]
        if E[xxx]<0.0:
            print "E array contains negative values"
            return neg_vals[0],neg_vals[1],neg_vals[2],neg_vals[3],neg_vals[4],neg_vals[5],neg_vals[6],neg_vals[7],neg_vals[8],neg_vals[9],neg_vals[10],neg_vals[11],neg_vals[12]

        
    #---------------WUE vs. NUE---------------#    


    wue=A/E*1000.0 #multiply by 1000 to get from umol CO2/umol H20 to umol CO2/mmol H20
    nue=A/na   

    
    return wue, nue, A, E, na, cs, ci, gsw, gs, gbw, gb, cc    
    
  
    
    
    
    
    
#############------Photosynthesis Model 4------#############

def photo_bound_meso_eqstom(tk_25,ekc,eko,etau,ev,ej,toptv,toptj,na, qeff, PAR,tl,ea,chl,ij,kc25,ko25,o,ca,rh,m,a,frnr,flnr,ra,jm,g0,b,dia,u,q,vwc_min,vwc_max,vwc):

    
    if all(ij>1.0):
        ij=np.zeros(shape=2)+1.0
    else:
        ij=ij
    
    ##---Calculated Parameter Arrays for Model(Constant+Variable Plant Trait(s))---##
    pa_v=611*np.exp((17.27*tl)/(tl+237.3)) #calculate saturation vapor pressure of surface (Pa)
    es_str=pa_con_atmfrac(pa_v,3528) #calculate saturation vapor pressure of surface (Pa-->umol/mol)
    dd=es_str-ea #calculate vapor pressure deficit (umol H2O/mol air)

    
    vmaxopt=frnr*flnr*ra*na #optimal carboxylation rate, limited by CO2 (umol CO2/m2s)
    jmaxopt=vmaxopt*jm #optimal electron transport rate, limited by RuBP (umol electrons/m2s)
    
    
    ##---Temperature Effects on Parameters---##

    if tl<=0.0:
        nue=[0.0]
        wue=[0.0]
        A=[0.0]
        E=[0.0]
        cs=[0.0]
        ci=[0.0]
        gsw=[0.0]
        gs=[0.0]
        gbw=[0.0]
        gb=[0.0]
        gm=[0.0]
        cc=[0.0]
        dd=[0.0]
        return wue, nue, A, E, cs, ci, gsw, gs, gbw,gb,gm,cc,dd

    if math.isnan(tl):
        nue=[np.nan]
        wue=[np.nan]
        A=[np.nan]
        E=[np.nan]
        cs=[np.nan]
        ci=[np.nan]
        gsw=[np.nan]
        gs=[np.nan]
        gbw=[np.nan]
        gb=[np.nan]
        gm=[np.nan]
        cc=[np.nan]
        dd=[np.nan]
        return wue, nue, A, E, cs, ci, gsw, gs, gbw,gb,gm,cc,dd 
        
    #calculated parameters due to temperature
    kc=arr_temp(pa_con_atmfrac(kc25,3528),ekc,tk_25,tl+273.15) #Michaelis-Menten kinetic coefficient for carbon dioxide at leaf temperature (umol/mol)
    ko=arr_temp(pa_con_atmfrac(ko25,3528),eko,tk_25,tl+273.15) #Michaelis-Menten kinetic coefficient for oxygen at leaf temperature (umol/mol) 
    gamma=(kc*o)/(2.0*ko)*0.21 #carbon dioxide compensation point (umol/mol)
    vmax1=bol_temp(vmaxopt,ev,toptv,tl+273.15) #carboxylation rate at leaf temperature, limited by CO2 (umol CO2/m2s)
    jmax1=bol_temp(jmaxopt,ej,toptj,tl+273.15) #rate of electron transport at leaf temperature, limited by RuBP (umol CO2/m2s)

    ##---Soil Moisture Effect on Parameters---##
    
    #below I removed the vwc constraint on photosynthesis because it is not a leaf trait
 
    if math.isnan(vwc):
        nue=[np.nan]
        wue=[np.nan]
        A=[np.nan]
        E=[np.nan]
        cs=[np.nan]
        ci=[np.nan]
        gsw=[np.nan]
        gs=[np.nan]
        gbw=[np.nan]
        gb=[np.nan]
        gm=[np.nan]
        cc=[np.nan]
        dd=[np.nan]
        return wue, nue, A, E, cs, ci, gsw, gs, gbw,gb,gm,cc,dd    
    
    if vwc<=vwc_min:
        nue=[0.0]
        wue=[0.0]
        A=[0.0]
        E=[0.0]
        cs=[0.0]
        ci=[0.0]
        gsw=[0.0]
        gs=[0.0]
        gbw=[0.0]
        gb=[0.0]
        gm=[0.0]
        cc=[0.0]
        dd=[0.0]
        return wue, nue, A, E, cs, ci, gsw, gs, gbw,gb,gm,cc,dd
    
    if vwc>=vwc_max:
        Wfac=1
    elif vwc<vwc_max:
        Wfac=((vwc-vwc_min)/(vwc_max-vwc_min))**q
    
    
    vmax=Wfac*vmax1
    jmax=Wfac*jmax1

    if vmax<=0:
        nue=[0.0]
        wue=[0.0]
        A=[0.0]
        E=[0.0]
        cs=[0.0]
        ci=[0.0]
        gsw=[0.0]
        gs=[0.0]
        gbw=[0.0]
        gb=[0.0]
        gm=[0.0]
        cc=[0.0]
        dd=[0.0]
        return wue, nue, A, E, cs, ci, gsw, gs, gbw,gb,gm,cc,dd
        
    
    ##---Determine J---## 
    alpha=(chl/1000.)/((chl/1000.)+0.076) #from Developmental Constratins on Photosynthesis: Effects of Light and Nutrition (Evans)
    qalpha=alpha*qeff
    Iphoton=PAR*ij
    
    if all(jmax>0.0):
        jj=(qalpha*Iphoton)/(np.sqrt(1.+(((qalpha**2)*(Iphoton**2))/(jmax**2))))

    
    ##---Define a1 and a2 depending on whether plant is rubisco limited or light limited---##

    #rubisco limited
    a1_r=vmax
    a2_r=kc*(1+(o/ko))
    #light limited
    a1_l=jj/4
    a2_l=2*gamma
    
    
     ##---Photosynthesis Model---##
    

    #Solve for Assimilation Using Cubic Equation
    
    #constants
    gbw=(27.0)/(200.0*np.sqrt(dia/u)) #26 is from: PV=nRT; n/V=P/RT; g(mol/m2s)=g(m/s)*P/RT (mol/m3); @ atm pressure at 25C, g(mol/m2s)=0.04g(mm/s)
    gb=b*gbw
    
    C1=(-a*m*rh*gb)+(a*g0)+gb
    C2=(ca*(gb**2)*a*m*rh)-(ca*gb*a*g0)-(a*g0*ca*gb)-(ca*(gb**2))
    C3=(ca**2)*(gb**2)*a*g0
    C4=(a*m*rh*(gb**2))-(a*g0*gb)
    C5=a*g0*ca*(gb**2)
    C11=C1+gb
    C12=C2-(ca*(gb**2))
    
    
    
    ##---Rubisco-Limited Assimilation---##
    
    
    #make list of coefficients for cubic eq.
    coef=[]
    for i in range(len(C1)):
        coef+=[[C11[i],C12[i]+C4[i]*a2_r[i]-a1_r[i]*C11[i],C3[i]+C5[i]*a2_r[i]-a1_r[i]*C12[i]+C4[i]*a1_r[i]*gamma[i],-a1_r[i]*C3[i]+C5[i]*a1_r[i]*gamma[i]]]
    
    #Solve Cubic Eq. for Roots
    roots=[]
    for i in range(len(coef)):
        roots+=[np.roots(coef[i])]

    
    #make list of real roots
    roots_real=[]
    for i in range(len(roots)):
        roots_real+=[roots[i].real[abs(roots[i].imag)<1e-5]]
    
    
    #make list of roots greater than zero
    roots_pos=[]
    for i in range(len(roots_real)):  
        roots_pos_0=[]
        for ii in range(len(roots_real[i])):
            if roots_real[i][ii]>0.0:
                roots_pos_0+=[roots_real[i][ii]]
            else:
                continue
        roots_pos+=[roots_pos_0]
    
   
    
    #assimilation is the minimum of the roots that are greater than zero
    A_r=[]
    for i in range(len(roots_pos)):
        A_r+=[np.min(roots_pos[i])]
    
    
    
    
    ##---Light-Limited Assimilation---##
    #make list of coefficients for cubic eq.
    coef=[]
    for i in range(len(C1)):
        coef+=[[C11[i],C12[i]+C4[i]*a2_l[i]-a1_l[i]*C11[i],C3[i]+C5[i]*a2_l[i]-a1_l[i]*C12[i]+C4[i]*a1_l[i]*gamma[i],-a1_l[i]*C3[i]+C5[i]*a1_l[i]*gamma[i]]]
    
    #Solve Cubic Eq. for Roots
    roots=[]
    for i in range(len(coef)):
        roots+=[np.roots(coef[i])]
    
    #make list of real roots
    roots_real=[]
    for i in range(len(roots)):
        roots_real+=[roots[i].real[abs(roots[i].imag)<1e-5]] #or do: roots[np.isreal(roots)]
    
    #make list of roots greater than zero
    roots_pos=[]
    for i in range(len(roots_real)):  
        roots_pos_0=[]
        for ii in range(len(roots_real[i])):
            if roots_real[i][ii]>0.0:
                roots_pos_0+=[roots_real[i][ii]]
            else:
                continue
        roots_pos+=[roots_pos_0]
    
    #assimilation is the minimum of the roots that are greater than zero
    A_l=[]
    for i in range(len(roots_pos)):
        A_l+=[np.min(roots_pos[i])]
                         
    ##---Determine Rate-Limiting Assimilation---##
    A=[]
    for i in range(len(A_r)):
        if A_r[i]<A_l[i]:
            A+=[A_r[i]] #rubisco limited
        elif A_l[i]<A_r[i]:
            A+=[A_l[i]] #light limited
        else: 
            A+=[A_l[i]] #both light and rubisco limited         
    
    
    

    ##---Solve for Stomatal Conductance to Water---##
    cs=ca-(A/gb) #stomatal co2 (umol CO2/mol)
    gsw=(m*A*rh/cs)+g0 #stomatal conductance to water (mol H2O/m2s) #make array from list
    
    
    ##---Solve for Stomatal Conductance to CO2---##
    gs=gsw*a
    
    ##---Solve for Internal CO2 Concentration (in Mesophyll)---##
    ci=cs-(A/gs)
    
    ##---Solve for Mesophyll Conductance to CO2---##
    gm=((a*m*A*rh)+(a*g0*cs))/cs
    
    ##---Solve for Internal CO2 Concentration (in Chloroplast)---##    
    cc=ci-(A/gm)    
    
    
    ##---Solve for Evapotranspiration---##
    E=gsw*dd #(umol H2O/m2s)

    
    #---------------Test for Nan or Negative Values---------------#       
    
    neg_vals=([np.zeros(1)+0]*13)
    
    for xxx in range(len(A)):
        if np.isnan(A[xxx]):
            print "A array contains nan values"
            return neg_vals[0],neg_vals[1],neg_vals[2],neg_vals[3],neg_vals[4],neg_vals[5],neg_vals[6],neg_vals[7],neg_vals[8],neg_vals[9],neg_vals[10],neg_vals[11],neg_vals[12]
        if A[xxx]<0.0:
            print "A array contains negative values"
            return neg_vals[0],neg_vals[1],neg_vals[2],neg_vals[3],neg_vals[4],neg_vals[5],neg_vals[6],neg_vals[7],neg_vals[8],neg_vals[9],neg_vals[10],neg_vals[11],neg_vals[12]
        if np.isnan(gsw[xxx]):
            print "gsw array contains nan values"
            return neg_vals[0],neg_vals[1],neg_vals[2],neg_vals[3],neg_vals[4],neg_vals[5],neg_vals[6],neg_vals[7],neg_vals[8],neg_vals[9],neg_vals[10],neg_vals[11],neg_vals[12]
        if gsw[xxx]<0.0:
            print "gsw array contains negative values"
            return neg_vals[0],neg_vals[1],neg_vals[2],neg_vals[3],neg_vals[4],neg_vals[5],neg_vals[6],neg_vals[7],neg_vals[8],neg_vals[9],neg_vals[10],neg_vals[11],neg_vals[12]
        if np.isnan(E[xxx]):
            print "E array contains nan values"
            return neg_vals[0],neg_vals[1],neg_vals[2],neg_vals[3],neg_vals[4],neg_vals[5],neg_vals[6],neg_vals[7],neg_vals[8],neg_vals[9],neg_vals[10],neg_vals[11],neg_vals[12]
        if E[xxx]<0.0:
            print "E array contains negative values"
            return neg_vals[0],neg_vals[1],neg_vals[2],neg_vals[3],neg_vals[4],neg_vals[5],neg_vals[6],neg_vals[7],neg_vals[8],neg_vals[9],neg_vals[10],neg_vals[11],neg_vals[12]

        
    #---------------WUE vs. NUE---------------#    

    wue=A/E*1000.0 #multiply by 1000 to get from umol CO2/umol H20 to umol CO2/mmol H20
    nue=A/na   
    
    return wue, nue, A, E, cs, ci, gsw, gs, gbw, gb, gm, cc, dd    