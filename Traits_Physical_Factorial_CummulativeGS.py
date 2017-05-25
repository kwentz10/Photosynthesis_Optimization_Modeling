#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 25 10:10:28 2017

@author: Katherine
"""

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue May 23 16:04:21 2017

@author: Katherine
"""

# -*- coding: utf-8 -*-
"""
Photosynthesis and Stomatal Conductance Model 
Created 9/27/2016
Katherine Wentz

This is a program that runs photosynthesis and
stomatal conductance models given changes in leaf-
level traits. 

The end product is graphs of NUE vs. WUE.


Update: I am going to run the model for plants with 
traits that are distinctive of the meadow moisture 
gradient in the alpine tundra.

Fix: correct for atmospheric pressure differences in co2, o2, and vapor pressure

Fix: vcmax temp dependence (pg 63 in plant physiological ecology book)

Fix: NEW VARIBALE TRAIT-->make the fraction of leaf N in rubisco go down with increasing SLA,
chlorophyll content, and decreasing light (wet meadow)--more N is allocated
to thylakoids. The only way for chl/m2 to increase even when g N/m2 goes down
or is constant is for the leaf to allocate more of leaf N to chl...also, note
that there is more organic N designated to photo in leaf when SLA goes up
because less N is used in structure. see "Photosynthesis or persistence: N allocation
in leaves of evergreen and deciduous... by Takashima et al. 2004. Also see Photosynthetic
nitrogen-use efficiency of species...by Poorter and Evans 1998

Note to self: NUE and WUE relationship flipflops with change in air temperature;
NUE makes sense because C:N decreases from dry to wet meadows; WUE increasing
in snowbed does not necessarilly make sense--look in the literature for this

herbs have a higher NUE

"""

#---------------Import Modules---------------#

import itertools as it
import numpy as np
from matplotlib import pyplot as plt

#Import combinations of variable parameters 
from Traits_Physical_Factorial_Inputs import leaf_params

#Import photosynthesis model
from Photosynthesis_Model import photo_bound_meso_eqstom as photo

#Import functions to switch between Pa and umol/mol at sea level
from photo_functions import pa_con_atmfrac


#import timeseries of vwc and temp
from time_dep_params import surtemp_dm, surtemp_mm, surtemp_wm, vwc_dm, vwc_mm, vwc_wm, na_dm, na_mm, na_wm


#---------------Determine if I Want to Keep Any of the Variable Parameters Constant---------------#

const_params=[]
for xxx in it.combinations(['ht'],0): #keep ht and t constant for constant vpd
    const_params+=[xxx]

#do this when I do not put any of the variable parameters as constant. instead I 
#vary each parameter one at a time while keeping the other parameters constant.
if const_params==[()]:
    const_params=[[-999999]]   

#---------------Begin Looping Through Photosynthesis Model---------------#

#each loop is for a constant value, or combinatin of constant values, of variable parameter as determined above
for ii in range(len(const_params)):




    #---------------Run through time series---------------#
    
    days=np.linspace(1,365,365)
    
    #dry meadow
    tot_nue_dm_avg=[]
    tot_wue_dm_avg=[]
    tot_nue_dm_min=[]
    tot_wue_dm_min=[]
    tot_nue_dm_max=[]
    tot_wue_dm_max=[]
    tot_A_dm_avg=[]
    tot_gs_dm_avg=[]
    
    #moist meadow
    tot_nue_mm_avg=[]
    tot_wue_mm_avg=[]
    tot_nue_mm_min=[]
    tot_wue_mm_min=[]
    tot_nue_mm_max=[]
    tot_wue_mm_max=[] 
    tot_A_mm_avg=[]
    tot_gs_mm_avg=[]

    #wet meadow
    tot_nue_wm_avg=[]
    tot_wue_wm_avg=[]
    tot_nue_wm_min=[]
    tot_wue_wm_min=[]
    tot_nue_wm_max=[]
    tot_wue_wm_max=[]
    tot_A_wm_avg=[]
    tot_gs_wm_avg=[]


        

    
    



    #---------------Photosynthesis + Stomatal Conductance Model---------------#

    
    ##---Constant Parameter Arrays for Model---##

    #----Params Used in Model Currently----#
      
    tk_25=298.16; #absolute temperature at 25 C
    ekc=80500.0 #Activation energy for K of CO2 (J mol-1)
    eko=14500.0 #Activation energy for K of O2 (J mol-1)
    etau=-29000.0  #Activation energy for tau (???) (J mol-1)
    ev=55000.0 #Activation energy for carboxylation (J mol-1)
    ej=55000.0 #Activation energy for electron transport (J mol-1)
    toptv=303.0 #Optimum temperature for maximum carboxylation (K)
    toptj=303.0 #Optimum temperature for maximum electron transport (K)
    ra=np.zeros(shape=1)+20.7 #specific rubisco activity (umol CO2/g Rub s)
    flnr=np.zeros(shape=1)+0.1 #fraction of leaf nitrogen in rubisco (g N Rub/g N leaf)
    frnr=np.zeros(shape=1)+6.25 #weight fraction of nitrogen in rubisco molecule (g Rub/g N Rub) 
    rh=np.zeros(shape=1)+0.5 #relative humidity (kPa/kPa)
    ca=np.zeros(shape=1)+405 #ambient carbon dioxide (umol CO2/mol air)
    ko25=np.zeros(shape=1)+30000 #Michaelis-Menten kinetic coefficient for oxygen at 25 C(Pa) 
    kc25=np.zeros(shape=1)+30 #Michaelis-Menten kinetic coefficient for carbon dioxide at 25 C (Pa)
    o=np.zeros(shape=1)+210000 #concentration of ambient oxygen (umol/mol)
    g0=np.zeros(shape=1)+0.002 #Ball-Berry stomatal conductance intercept parameter (mol H2O/m2s)
    a=np.zeros(shape=1)+1.6 #Conversion Coefficient between stomatal conductance to water and carbon dioxide (unitless)
    ij=np.zeros(shape=1)+1.0 #leaf angle index--downregulates jmax
    m=np.zeros(shape=1)+9.0 #ball-berry parameter (unitless)
    b=1.37 #Conversion Coefficient between boundary layer conductance to water and carbon dioxide 
    u=5.0 #windspeed (m/s)
    qeff=0.32 #leaf quantum yield, electrons
    PAR=2000 #photosynthetic active radiation (umol/m2s)
    jm=2.68 #slope coefficient 
    vwc_min=0.08 #minimum soil water content for photosynthesis to occur (permanent wilting point) (cm3/cm3) 
    vwc_max=0.68 #maximum soil water content where increases in soil water do not affect photosynthesis (field capacity?) (cm3/cm3)
    q=0.2 #parameter for soil water affect on photosynthesis (unitless)
   
    
    #------constant variable params for sensitivty analysis-----#
    
    chl_c=np.zeros(shape=1)+(np.mean([396,465,476])) #Chlorophyll Content of the Leaf (umol chl/m2)
    ht_c=np.zeros(shape=1)+(np.mean([9.2,19.5,20.0])) #Temperature of the Leaf (K)
    dia_c=np.zeros(shape=1)+(np.mean([1.4,2.3,2.6])/100.) #Mean diameter or size of leaf (m)
    na_c=np.zeros(shape=1)+(np.mean([2.5,5.6,6.3])) #leaf nitrogen (g N/ m2)
    t_c=np.zeros(shape=1)+15.0 #temp (C)


    #-----which timeseries should I use--based on factorial meadow type---#
    
    na_type=[na_dm, na_mm, na_wm,na_dm, na_mm, na_wm,na_dm, na_mm, na_wm]
    vwc_type=[vwc_dm,vwc_dm,vwc_dm,vwc_mm,vwc_mm,vwc_mm,vwc_wm,vwc_wm,vwc_wm]
    temp_type=[surtemp_dm,surtemp_dm,surtemp_dm,surtemp_mm,surtemp_mm,surtemp_mm,surtemp_wm,surtemp_wm,surtemp_wm]


    A_tot_all=[]

#---------------Import Variable Parameter Arrays from Leaf Parameter File---------------#
    
    for xx in range(len(leaf_params)):
        for key,val in leaf_params[xx].items():
            exec(key + '=val')

            
                             
        #set variable parameters constant if I specify this above
        if 'na' in const_params[ii]:
            na=na_c
        if 'dia' in const_params[ii]:
            dia=dia_c
        if 'chl' in const_params[ii]:
            chl=chl_c
        if 'ht' in const_params[ii]:
            ht=ht_c

        
        A_tot=0
        
        for time in range(len(surtemp_dm)):


      
            #------calculate vapor pressure-----#
            pa_v=611*np.exp((17.27*temp_type[xx][time])/(temp_type[xx][time]+237.3)) #saturation vapor pressure of air (Pa)
            ea_str=pa_con_atmfrac(pa_v,3528) #saturation vapor pressure of air (Pa-->umol h20/mol air)
            ea=rh*ea_str #vapor pressure (umol h2O/mol air)                


            #correct for leaf temperatures using leaf height
   
            t_diff=18-0.4*ht
        
            tl=temp_type[xx][time]+t_diff                
            
            #---------------Photosynthesis Function---------------#
        
            #alter this line of code for when implementing different photosynthesis functions
            wue, nue, A, E, cs, ci, gsw, gs, gbw, gb, gm, cc,dd =photo(tk_25,ekc,eko,etau,ev,ej,toptv,toptj,na, qeff, PAR,tl,ea,chl,ij,kc25,ko25,o,ca,rh,m,a,frnr,flnr,ra,jm,g0,b,dia,u,q,vwc_min,vwc_max,vwc_type[xx][time])
 
       
            #test to make sure wue and nue are positive at not 'nan'
            if wue[0]==-999 and nue[0]==-999:
           
                continue                                
            
            if np.isnan(A[0]):
                A[0]=0.0

            A_tot+=(A[0]*3600*6)/1000000.*44.
                
       
        A_tot_all+=[A_tot]

        
        
                    
    

