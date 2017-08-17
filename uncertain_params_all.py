#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 20 12:08:50 2017

@author: Katherine
"""

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May 24 09:04:31 2017

@author: Katherine

Draw parameter values from uniform distribution
"""

import numpy as np

def monte_carlo(chl_mean, chl_sd, dia_mean, dia_sd, na_min, na_max, ht_mean, ht_sd):

    #meadow specific parameter min and max
    

    temp_mean=[15.+2.,15.,15.-2.]
    temp_sd=[2.5,2.5,2.5]
    vwc_mean=[0.12,0.15,0.27]
    vwc_sd=[0.008,0.013,0.022]
    
    
    #parameter min and max for all meadows
    ekc_lw=80500.0*0.8 #Activation energy for K of CO2 (J mol-1) VARIABLE
    ekc_hg=80500.0*1.2 #Activation energy for K of CO2 (J mol-1) VARIABLE
    
    eko_lw=14500.0*0.8 #Activation energy for K of O2 (J mol-1) VARIABLE
    eko_hg=14500.0*1.2 #Activation energy for K of O2 (J mol-1) VARIABLE
    
    
    etau_lw=-29000.0*0.8  #Activation energy for tau (???) (J mol-1) VARIABLE
    etau_hg=-29000.0*1.2  #Activation energy for tau (???) (J mol-1) VARIABLE
    
    
    ev_lw=55000.0*0.8 #Activation energy for carboxylation (J mol-1) VARIABLE
    ev_hg=55000.0*1.2 #Activation energy for carboxylation (J mol-1) VARIABLE
    
    ej_lw=55000.0*0.8 #Activation energy for electron transport (J mol-1) VARIABLE
    ej_hg=55000.0*1.2 #Activation energy for electron transport (J mol-1) VARIABLE
    
    ra_lw=20.7*0.8 #specific rubisco activity (umol CO2/g Rub s) VARIABLE
    ra_hg=20.7*1.2 #specific rubisco activity (umol CO2/g Rub s) VARIABLE
    
    flnr_lw=0.1*0.8 #fraction of leaf nitrogen in rubisco (g N Rub/g N leaf) VARIABLE
    flnr_hg=0.1*1.2 #fraction of leaf nitrogen in rubisco (g N Rub/g N leaf) VARIABLE
    
    rh_lw=0.5*0.8 #relative humidity (kPa/kPa) VARIABLE
    rh_hg=0.5*1.2 #relative humidity (kPa/kPa) VARIABLE
    
    
    ko25_lw=30000*0.8 #Michaelis-Menten kinetic coefficient for oxygen at 25 C(Pa) VARIABLE
    ko25_hg=30000*1.2 #Michaelis-Menten kinetic coefficient for oxygen at 25 C(Pa) VARIABLE
    
    
    kc25_lw=30*0.8 #Michaelis-Menten kinetic coefficient for carbon dioxide at 25 C (Pa) VARIABLE
    kc25_hg=30*1.2 #Michaelis-Menten kinetic coefficient for carbon dioxide at 25 C (Pa) VARIABLE
    
    
    g0_lw=0.002*0.8 #Ball-Berry stomatal conductance intercept parameter (mol H2O/m2s) VARIABLE
    g0_hg=0.002*1.2 #Ball-Berry stomatal conductance intercept parameter (mol H2O/m2s) VARIABLE
    
    
    m_lw=9.0*0.8 #ball-berry parameter (unitless) VARIABLE
    m_hg=9.0*1.2 #ball-berry parameter (unitless) VARIABLE
    
    
    u_lw=5.0*0.8 #windspeed (m/s) VARIABLE
    u_hg=5.0*1.2 #windspeed (m/s) VARIABLE
    
    qeff_lw=0.32*0.8 #leaf quantum yield, electrons VARIABLE
    qeff_hg=0.32*1.2 #leaf quantum yield, electrons VARIABLE
    
    PAR_lw=2000*0.8 #photosynthetic active radiation (umol/m2s) VARIABLE
    PAR_hg=2000*1.2 #photosynthetic active radiation (umol/m2s) VARIABLE
    
    
    jm_lw=2.68*0.8 #slope coefficient  VARIABLE
    jm_hg=2.68*1.2 #slope coefficient  VARIABLE
    
    
    vwc_min_lw=0.08*0.8 #minimum soil water content for photosynthesis to occur (permanent wilting point) (cm3/cm3)  VARIABLE
    vwc_min_hg=0.08*1.2 #minimum soil water content for photosynthesis to occur (permanent wilting point) (cm3/cm3)  VARIABLE
    
    vwc_max_lw=0.66*0.8 #maximum soil water content where increases in soil water do not affect photosynthesis (field capacity?) (cm3/cm3) VARIABLE
    vwc_max_hg=0.66*1.2 #maximum soil water content where increases in soil water do not affect photosynthesis (field capacity?) (cm3/cm3) VARIABLE
    
    
    q_lw=0.2*0.8 #parameter for soil water affect on photosynthesis (unitless) VARIABLE
    q_hg=0.2*1.2 #parameter for soil water affect on photosynthesis (unitless) VARIABLE
    
    #make dictionary of parameters with uncertainty
    
    params=[]
    
    for i in range(len(chl_mean)):
        meadow_params=[]
        for ii in range(30):
            d={} #create dictionary
            
            #parameters unique to each meadow
    
            d["chl"]=np.random.normal(chl_mean[i],chl_sd[i],1)

            d["na"]=np.random.uniform(na_min[i],na_max[i],1)
            d["dia"]=np.random.normal(dia_mean[i],dia_sd[i],1)
            if d["dia"]<0.0:
                d["dia"]=[dia_mean[i]-dia_sd[i]]
            d["ht"]=np.random.normal(ht_mean[i],ht_sd[i],1)
            d["temp"]=np.random.normal(temp_mean[i],temp_sd[i],1)
            d["vwc"]=np.random.normal(vwc_mean[i],vwc_sd[i],1)
    
    
            #parameters for all meadows
            d["ekc"]=np.random.uniform(ekc_lw,ekc_hg,1)
            d["eko"]=np.random.uniform(eko_lw,eko_hg,1)
            d["etau"]=np.random.uniform(etau_lw,etau_hg,1)
            d["ev"]=np.random.uniform(ev_lw,ev_hg,1)
            d["ej"]=np.random.uniform(ej_lw,ej_hg,1)
            d["ra"]=np.random.uniform(ra_lw,ra_hg,1)
            d["flnr"]=np.random.uniform(flnr_lw,flnr_hg,1)
            d["rh"]=np.random.uniform(rh_lw,rh_hg,1)
            d["ko25"]=np.random.uniform(ko25_lw,ko25_hg,1)
            d["kc25"]=np.random.uniform(kc25_lw,kc25_hg,1)
            d["g0"]=np.random.uniform(g0_lw,g0_hg,1)
            d["m"]=np.random.uniform(m_lw,m_hg,1)
            d["u"]=np.random.uniform(u_lw,u_hg,1)
            d["qeff"]=np.random.uniform(qeff_lw,qeff_hg,1)
            d["PAR"]=np.random.uniform(PAR_lw,PAR_hg,1)
            d["jm"]=np.random.uniform(jm_lw,jm_hg,1)
            d["vwc_min"]=np.random.uniform(vwc_min_lw,vwc_min_hg,1)
            d["vwc_max"]=np.random.uniform(vwc_max_lw,vwc_max_hg,1)
            d["q"]=np.random.uniform(q_lw,q_hg,1)
            
            meadow_params+=[d]
        
        params+=[meadow_params]
