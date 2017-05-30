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
from uncertain_params import params

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

 
    #---------------Initialize Plots---------------#
    
    #--figure 1--#
         
    #put in correct ax value (e.g. axA, axB)
    fig1,axA = plt.subplots(figsize=(30,15))
    
    #twin axis
#    axA2=axA.twinx()

    axA.set_xlabel('Time (days)',fontsize=28, fontname='Times New Roman')
    axA.set_ylabel('NUE (umol CO2/g N s)',fontsize=28, fontname='Times New Roman')
#    axA2.set_ylabel('Temperature ($^\circ$C)', fontsize=28, fontname='Times New Roman')  
#    axA.set_ylim([0,3])
#    axA2.set_ylim([0,25])
    axA.set_xlim([0,365])
    axA.set_title('NUE in the Dry Meadow', fontname='Times New Roman',fontsize=36,fontweight='bold')

     #-----figure 2----#

        
    #put in correct ax value (e.g. axA, axB)
    fig2,axB = plt.subplots(figsize=(30,15))
    
    #twin axis
#    axB2=axB.twinx()

    axB.set_xlabel('Time (days)',fontsize=28, fontname='Times New Roman')
    axB.set_ylabel('NUE (umol CO2/g N s)',fontsize=28, fontname='Times New Roman')
#    axB2.set_ylabel('Temperature ($^\circ$C)', fontsize=28, fontname='Times New Roman')  
#    axB.set_ylim([0,5])
#    axB2.set_ylim([0,20])
    axB.set_xlim([0,365])
    axB.set_title('NUE in the Moist Meadow', fontname='Times New Roman',fontsize=36,fontweight='bold')

     #-----figure 3----#
          
    #put in correct ax value (e.g. axA, axB)
    fig3,axC = plt.subplots(figsize=(30,15))
    
    #twin axis
#    axC2=axC.twinx()

    axC.set_xlabel('Time (days)',fontsize=28, fontname='Times New Roman')
    axC.set_ylabel('NUE (umol CO2/g N s)',fontsize=28, fontname='Times New Roman')
#    axC2.set_ylabel('Temperature ($^\circ$C)', fontsize=28, fontname='Times New Roman')  
#    axC.set_ylim([0,8])
#    axC2.set_ylim([0,20])
    axC.set_xlim([0,365])
    axC.set_title('NUE in the Wet Meadow', fontname='Times New Roman',fontsize=36,fontweight='bold')


    #--figure 4--#
         
    #put in correct ax value (e.g. axA, axB)
    fig4,axD = plt.subplots(figsize=(30,15))

    axD.set_xlabel('Time (days)',fontsize=28, fontname='Times New Roman')
    axD.set_ylabel('WUE (umol CO2/mmol H2O)', fontsize=28, fontname='Times New Roman')  
    axD.set_title('WUE in the Dry Meadow', fontname='Times New Roman',fontsize=36,fontweight='bold')
    axD.set_xlim([0,365])

     #-----figure 5----#

        
    #put in correct ax value (e.g. axA, axB)
    fig5,axE = plt.subplots(figsize=(30,15))
    
    axE.set_xlabel('Time (days)',fontsize=28, fontname='Times New Roman')
    axE.set_ylabel('WUE (umol CO2/mmol H2O)', fontsize=28, fontname='Times New Roman')  
    axE.set_title('WUE in the Moist Meadow', fontname='Times New Roman',fontsize=36,fontweight='bold')
    axE.set_xlim([0,365])

     #-----figure 6----#
          
    #put in correct ax value (e.g. axA, axB)
    fig6,axF = plt.subplots(figsize=(30,15))

    axF.set_xlabel('Time (days)',fontsize=28, fontname='Times New Roman')
    axF.set_ylabel('WUE (umol CO2/mmol H2O)', fontsize=28, fontname='Times New Roman')  
    axF.set_title('WUE in the Wet Meadow', fontname='Times New Roman',fontsize=36,fontweight='bold')
    axF.set_xlim([0,365])

    #--figure 7--#
         
    #put in correct ax value (e.g. axA, axB)
    fig7,axG = plt.subplots(figsize=(30,15))

    axG.set_xlabel('Time (days)',fontsize=28, fontname='Times New Roman')
    axG.set_ylabel('Assimilation (umol CO2/m2s)', fontsize=28, fontname='Times New Roman')  
    axG.set_title('Assimilation in the Dry Meadow', fontname='Times New Roman',fontsize=36,fontweight='bold')
    axG.set_xlim([0,365])
    
     #-----figure 8----#

        
    #put in correct ax value (e.g. axA, axB)
    fig8,axH = plt.subplots(figsize=(30,15))
    
    axH.set_xlabel('Time (days)',fontsize=28, fontname='Times New Roman')
    axH.set_ylabel('Assimilation (umol CO2/m2s)', fontsize=28, fontname='Times New Roman')  
    axH.set_title('Assimilation in the Moist Meadow', fontname='Times New Roman',fontsize=36,fontweight='bold')
    axH.set_xlim([0,365])
     #-----figure 9----#
          
    #put in correct ax value (e.g. axA, axB)
    fig9,axI = plt.subplots(figsize=(30,15))

    axI.set_xlabel('Time (days)',fontsize=28, fontname='Times New Roman')
    axI.set_ylabel('Assimilation (umol CO2/mmol H2O)', fontsize=28, fontname='Times New Roman')  
    axI.set_title('Assimilation in the Wet Meadow', fontname='Times New Roman',fontsize=36,fontweight='bold')
    axI.set_xlim([0,365])

    #--figure 10--#
         
    #put in correct ax value (e.g. axA, axB)
    fig10,axJ = plt.subplots(figsize=(30,15))

    axJ.set_xlabel('Time (days)',fontsize=28, fontname='Times New Roman')
    axJ.set_ylabel('Temperature ($^\circ$C)', fontsize=28, fontname='Times New Roman')  
    axJ.set_title('Surface Temperature in Three Plant Communities', fontname='Times New Roman',fontsize=36,fontweight='bold')
    axJ.set_xlim([0,365])

     #-----figure 11----#

        
    #put in correct ax value (e.g. axA, axB)
    fig11,axK = plt.subplots(figsize=(30,15))
    
    axK.set_xlabel('Time (days)',fontsize=28, fontname='Times New Roman')
    axK.set_ylabel('Volumetric Water Content (m3/m3)', fontsize=28, fontname='Times New Roman')  
    axK.set_title('Volumetric Water Content in Three Plant Communities', fontname='Times New Roman',fontsize=36,fontweight='bold')
    axK.set_xlim([0,365])

     #-----figure 12----#
          
    #put in correct ax value (e.g. axA, axB)
    fig12,axL = plt.subplots(figsize=(30,15))

    axL.set_xlabel('Time (days)',fontsize=28, fontname='Times New Roman')
    axL.set_ylabel('Leaf Nitrogen Content (g N/m2)', fontsize=28, fontname='Times New Roman')  
    axL.set_title('Leaf Nitrogen Content in Three Plant Communities', fontname='Times New Roman',fontsize=36,fontweight='bold')
    axL.set_xlim([0,365])

    #---------------Run through time series---------------#
    
    days=np.linspace(1,365,365)
    
    #dry meadow
    tot_nue_dm_avg=[]
    tot_wue_dm_avg=[]
    tot_A_dm_avg=[]    
    tot_nue_dm_min=[]
    tot_wue_dm_min=[]
    tot_A_dm_min=[]    
    tot_nue_dm_max=[]
    tot_wue_dm_max=[]
    tot_A_dm_max=[]

    
    #moist meadow
    tot_nue_mm_avg=[]
    tot_wue_mm_avg=[]
    tot_A_mm_avg=[]    
    tot_nue_mm_min=[]
    tot_wue_mm_min=[]
    tot_A_mm_min=[]    
    tot_nue_mm_max=[]
    tot_wue_mm_max=[] 
    tot_A_mm_max=[]


    #wet meadow
    tot_nue_wm_avg=[]
    tot_wue_wm_avg=[]
    tot_A_wm_avg=[]    
    tot_nue_wm_min=[]
    tot_wue_wm_min=[]
    tot_A_wm_min=[]    
    tot_nue_wm_max=[]
    tot_wue_wm_max=[]
    tot_A_wm_max=[]


    for time in range(len(surtemp_dm)):
        
        #total nue and wue
        nue_tot=[]
        wue_tot=[]
       
        #wue and nue arrays
        wue_d=[]
        nue_d=[]
    
        wue_m=[]
        nue_m=[]
           
        wue_w=[]
        nue_w=[]    
         
    
        #gsw arrays
    
        gsw_d=[]
         
        gsw_m=[]
    
        gsw_w=[]
    
    
        #assimilation arrays
    
        A_d=[]
    
        A_m=[]
    
        A_w=[]
    
    
        #evapo arrays
    
        E_d=[]
        
        E_m=[]
    
        E_w=[]
    
        
        #vapor pressure deficit arrays
        vpd_d=[]
    
        vpd_m=[]
    
        vpd_w=[]

    
    



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

    #---------------Import Variable Parameter Arrays from Leaf Parameter File---------------#
        
        for xx in range(len(params)):
            for yy in range(len(params[xx])):
                for key,val in params[xx][yy].items():
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
    
                
                #number of simulations per meadow type:
                m_sim=len(params) #meadow simulations 
                
                if xx==0: 
        
                    #------calculate vapor pressure-----#
                    pa_v=611*np.exp((17.27*surtemp_dm[time])/(surtemp_dm[time]+237.3)) #saturation vapor pressure of air (Pa)
                    ea_str=pa_con_atmfrac(pa_v,3528) #saturation vapor pressure of air (Pa-->umol h20/mol air)
                    ea=rh*ea_str #vapor pressure (umol h2O/mol air)                
        
        
                    #correct for leaf temperatures using leaf height
           
                    t_diff=18-0.4*ht
                
                    tl=surtemp_dm[time]+t_diff
    
                    
                    #---------------Photosynthesis Function---------------#
                
                    #alter this line of code for when implementing different photosynthesis functions
                    wue, nue, A, E, cs, ci, gsw, gs, gbw, gb, gm, cc,dd =photo(tk_25,ekc,eko,etau,ev,ej,toptv,toptj,na_dm[time], qeff, PAR,tl,ea,chl,ij,kc25,ko25,o,ca,rh,m,a,frnr,flnr,ra,jm,g0,b,dia,u,q,vwc_min,vwc_max,vwc_dm[time])
     
               
                    #test to make sure wue and nue are positive at not 'nan'
                    if wue[0]==-999 and nue[0]==-999:
                   
                        continue
                    
                    wue_d+=[wue[0]]
                    nue_d+=[nue[0]]
                    gsw_d+=[gsw[0]]
                    A_d+=[A[0]]
                    E_d+=[E[0]]
                    vpd_d+=[dd[0]]
    
                    
                elif xx==1:
                    
                    #------calculate vapor pressure-----#
                    pa_v=611*np.exp((17.27*surtemp_mm[time])/(surtemp_mm[time]+237.3)) #saturation vapor pressure of air (Pa)
                    ea_str=pa_con_atmfrac(pa_v,3528) #saturation vapor pressure of air (Pa-->umol h20/mol air)
                    ea=rh*ea_str #vapor pressure (umol h2O/mol air)                
        
        
                    #correct for leaf temperatures using leaf height
           
                    t_diff=18-0.4*ht
                
                    tl=surtemp_mm[time]+t_diff
                    
                    #---------------Photosynthesis Function---------------#
                
                    #alter this line of code for when implementing different photosynthesis functions
                    wue, nue, A, E, cs, ci, gsw, gs, gbw, gb, gm, cc,dd =photo(tk_25,ekc,eko,etau,ev,ej,toptv,toptj,na_mm[time], qeff, PAR,tl,ea,chl,ij,kc25,ko25,o,ca,rh,m,a,frnr,flnr,ra,jm,g0,b,dia,u,q,vwc_min,vwc_max,vwc_mm[time])
     
               
                    #test to make sure wue and nue are positive at not 'nan'
                    if wue[0]==-999 and nue[0]==-999:
                   
                        continue                
    
                    
                    wue_m+=[wue[0]]
                    nue_m+=[nue[0]]
                    gsw_m+=[gsw[0]]
                    A_m+=[A[0]]
                    E_m+=[E[0]]
                    vpd_m+=[dd[0]]
                    
                    if vwc_mm[time]==-999:
                        print "ahhhh"
                    
    
    
           
            
                elif xx==2:
    
                    #------calculate vapor pressure-----#
                    pa_v=611*np.exp((17.27*surtemp_wm[time])/(surtemp_wm[time]+237.3)) #saturation vapor pressure of air (Pa)
                    ea_str=pa_con_atmfrac(pa_v,3528) #saturation vapor pressure of air (Pa-->umol h20/mol air)
                    ea=rh*ea_str #vapor pressure (umol h2O/mol air)                
        
        
                    #correct for leaf temperatures using leaf height
           
                    t_diff=18-0.4*ht
                
                    tl=surtemp_wm[time]+t_diff                
                    
                    #---------------Photosynthesis Function---------------#
                
                    #alter this line of code for when implementing different photosynthesis functions
                    wue, nue, A, E, cs, ci, gsw, gs, gbw, gb, gm, cc,dd =photo(tk_25,ekc,eko,etau,ev,ej,toptv,toptj,na_wm[time], qeff, PAR,tl,ea,chl,ij,kc25,ko25,o,ca,rh,m,a,frnr,flnr,ra,jm,g0,b,dia,u,q,vwc_min,vwc_max,vwc_wm[time])
     
               
                    #test to make sure wue and nue are positive at not 'nan'
                    if wue[0]==-999 and nue[0]==-999:
                   
                        continue                                
    
                        
                    wue_w+=[wue[0]]
                    nue_w+=[nue[0]]
                    gsw_w+=[gsw[0]]
                    A_w+=[A[0]]
                    E_w+=[E[0]]
                    vpd_w+=[dd[0]]
    
                    
    
        #dry meadow
        tot_nue_dm_avg+=[np.mean(nue_d)]
        tot_wue_dm_avg+=[np.mean(wue_d)]
        tot_A_dm_avg+=[np.mean(A_d)]
        tot_nue_dm_min+=[np.min(nue_d)]
        tot_wue_dm_min+=[np.min(wue_d)]
        tot_A_dm_min+=[np.min(A_d)]
        tot_nue_dm_max+=[np.max(nue_d)]
        tot_wue_dm_max+=[np.max(wue_d)]
        tot_A_dm_max+=[np.max(A_d)]

        
        #moist meadow
        tot_nue_mm_avg+=[np.mean(nue_m)]
        tot_wue_mm_avg+=[np.mean(wue_m)]
        tot_A_mm_avg+=[np.mean(A_m)]
        tot_nue_mm_min+=[np.min(nue_m)]
        tot_wue_mm_min+=[np.min(wue_m)]
        tot_A_mm_min+=[np.min(A_m)]
        tot_nue_mm_max+=[np.max(nue_m)]
        tot_wue_mm_max+=[np.max(wue_m)]   
        tot_A_mm_max+=[np.max(A_m)]
  
    
        #wet meadow
        tot_nue_wm_avg+=[np.mean(nue_w)]
        tot_wue_wm_avg+=[np.mean(wue_w)]
        tot_A_wm_avg+=[np.mean(A_w)]
        tot_nue_wm_min+=[np.min(nue_w)]
        tot_wue_wm_min+=[np.min(wue_w)]
        tot_A_wm_min+=[np.min(A_w)]
        tot_nue_wm_max+=[np.max(nue_w)]
        tot_wue_wm_max+=[np.max(wue_w)]
        tot_A_wm_max+=[np.max(A_w)]        

    
#---------------Plot Plant Communities vs. NUE & WUE---------------#      



    #dry meadow nue 
    axA.plot(days, tot_nue_dm_avg, 'r-')
#    axA.plot(days,na_dm,'y',linewidth=5)
#    axA.plot(days,np.array(vwc_dm)*10.,'c-',linewidth=5)
#    axA2.plot(days, tot_wue_dm_avg, 'b-')
#    axA2.plot(days,tot_A_dm_avg,'k-',linewidth=5)
#    axA2.plot(days,np.array(surtemp_dm)/1.0,'g-',linewidth=5)
    axA.fill_between(days, tot_nue_dm_min, tot_nue_dm_max,alpha=0.3,color='red')
#    axA2.fill_between(days, tot_wue_dm_min, tot_wue_dm_max,alpha=0.3,color='blue')
 
    #moist meadow nue 
    axB.plot(days, tot_nue_mm_avg, 'r-')
#    axB.plot(days,na_mm,'y',linewidth=5)
#    axB.plot(days,np.array(vwc_mm)*10.,'c-',linewidth=5)   
#    axB2.plot(days, tot_wue_mm_avg, 'b-')
#    axB2.plot(days,tot_A_mm_avg,'k-',linewidth=5)
#    axB2.plot(days,np.array(surtemp_mm)/1.0,'g-',linewidth=5)
    axB.fill_between(days, tot_nue_mm_min, tot_nue_mm_max,alpha=0.3,color='red')
#    axB2.fill_between(days, tot_wue_mm_min, tot_wue_mm_max,alpha=0.3,color='blue')    

    #wet meadow nue 
    axC.plot(days, tot_nue_wm_avg, 'r-')
#    axC.plot(days,na_wm,'y',linewidth=5)    
#    axC.plot(days,np.array(vwc_wm)*10.,'c-',linewidth=5)    
#    axC2.plot(days, tot_wue_wm_avg, 'b-')
#    axC2.plot(days,tot_A_wm_avg,'k-',linewidth=5)    
#    axC2.plot(days,np.array(surtemp_wm)/1.0,'g-',linewidth=5)    
    axC.fill_between(days, tot_nue_wm_min, tot_nue_wm_max,alpha=0.3, color='red')
#    axC2.fill_between(days, tot_wue_wm_min, tot_wue_wm_max,alpha=0.3,color='blue')
 



###------------------------------------------------------######    
    
    #dry meadow wue 
#    axA.plot(days,na_dm,'y',linewidth=5)
#    axA.plot(days,np.array(vwc_dm)*10.,'c-',linewidth=5)
    axD.plot(days, tot_wue_dm_avg, 'b-')
#    axA2.plot(days,tot_A_dm_avg,'k-',linewidth=5)
#    axA2.plot(days,np.array(surtemp_dm)/1.0,'g-',linewidth=5)
    axD.fill_between(days, tot_wue_dm_min, tot_wue_dm_max,alpha=0.3,color='blue')
 
    #moist meadow wue 
#    axB.plot(days,na_mm,'y',linewidth=5)
#    axB.plot(days,np.array(vwc_mm)*10.,'c-',linewidth=5)   
    axE.plot(days, tot_wue_mm_avg, 'b-')
#    axB2.plot(days,tot_A_mm_avg,'k-',linewidth=5)
#    axB2.plot(days,np.array(surtemp_mm)/1.0,'g-',linewidth=5)
    axE.fill_between(days, tot_wue_mm_min, tot_wue_mm_max,alpha=0.3,color='blue')    

    #wet meadow wue 

#    axC.plot(days,na_wm,'y',linewidth=5)    
#    axC.plot(days,np.array(vwc_wm)*10.,'c-',linewidth=5)    
    axF.plot(days, tot_wue_wm_avg, 'b-')
#    axC2.plot(days,tot_A_wm_avg,'k-',linewidth=5)    
#    axC2.plot(days,np.array(surtemp_wm)/1.0,'g-',linewidth=5)    
    axF.fill_between(days, tot_wue_wm_min, tot_wue_wm_max,alpha=0.3,color='blue')


###------------------------------------------------------######    
    

    #dry meadow A
#    axA.plot(days,na_dm,'y',linewidth=5)
#    axA.plot(days,np.array(vwc_dm)*10.,'c-',linewidth=5)
    axG.plot(days, tot_A_dm_avg, 'k-')
#    axA2.plot(days,tot_A_dm_avg,'k-',linewidth=5)
#    axA2.plot(days,np.array(surtemp_dm)/1.0,'g-',linewidth=5)
    axG.fill_between(days, tot_A_dm_min, tot_A_dm_max,alpha=0.3,color='black')
 
    #moist meadow A 
#    axB.plot(days,na_mm,'y',linewidth=5)
#    axB.plot(days,np.array(vwc_mm)*10.,'c-',linewidth=5)   
    axH.plot(days, tot_A_mm_avg, 'k-')
#    axB2.plot(days,tot_A_mm_avg,'k-',linewidth=5)
#    axB2.plot(days,np.array(surtemp_mm)/1.0,'g-',linewidth=5)
    axH.fill_between(days, tot_A_mm_min, tot_A_mm_max,alpha=0.3,color='black')    

    #wet meadow A 

#    axC.plot(days,na_wm,'y',linewidth=5)    
#    axC.plot(days,np.array(vwc_wm)*10.,'c-',linewidth=5)    
    axI.plot(days, tot_A_wm_avg, 'k-')
#    axC2.plot(days,tot_A_wm_avg,'k-',linewidth=5)    
#    axC2.plot(days,np.array(surtemp_wm)/1.0,'g-',linewidth=5)    
    axI.fill_between(days, tot_A_wm_min, tot_A_wm_max,alpha=0.3,color='black')


###------------------------------------------------------######   

    axJ.plot(days,surtemp_dm,'k-',linewidth=5,label="dry meadow")
    axJ.plot(days,surtemp_mm,'g-',linewidth=5,label="moist meadow")
    axJ.plot(days,surtemp_wm,'b-',linewidth=5,label="wet meadow")
    axJ.legend()
    
    axK.plot(days,vwc_dm,'k-',linewidth=5,label="dry meadow")
    axK.plot(days,vwc_mm,'g-',linewidth=5,label="moist meadow")
    axK.plot(days,vwc_wm,'b-',linewidth=5,label="wet meadow")
    axK.legend()
    
    axL.plot(days,na_dm,'k-',linewidth=5,label="dry meadow")
    axL.plot(days,na_mm,'g-',linewidth=5,label="moist meadow")
    axL.plot(days,na_wm,'b-',linewidth=5,label="wet meadow")
    axL.legend()
    
    #color axes
    axA.tick_params(axis='y', labelsize=20)
#    axA.yaxis.label.set_color('red')
    axA.tick_params(axis='x', pad=15,labelsize=20)


    axB.tick_params(axis='y', labelsize=20)
#    axB.yaxis.label.set_color('red')  
    axB.tick_params(axis='x', pad=15,labelsize= 20)

    axC.tick_params(axis='y', labelsize=20)
#    axC.yaxis.label.set_color('red')  
    axC.tick_params(axis='x', pad=15,labelsize= 20)

#    axA2.tick_params(axis='y', labelsize=20,colors='blue')
#    axA2.yaxis.label.set_color('blue')
#
#    axB2.tick_params(axis='y', labelsize=20,colors='blue')
#    axB2.yaxis.label.set_color('blue')    
#    
#    axC2.tick_params(axis='y', labelsize=20,colors='blue')
#    axC2.yaxis.label.set_color('blue')   

    axD.tick_params(axis='y', labelsize=20)
#    axC.yaxis.label.set_color('red')
    axD.tick_params(axis='x', pad=15,labelsize=20)


    axE.tick_params(axis='y', labelsize=20)
#    axB.yaxis.label.set_color('red')  
    axE.tick_params(axis='x', pad=15,labelsize= 20)

    axF.tick_params(axis='y', labelsize=20)
#    axC.yaxis.label.set_color('red')  
    axF.tick_params(axis='x', pad=15,labelsize= 20)

    axG.tick_params(axis='y', labelsize=20)
#    axC.yaxis.label.set_color('red')
    axG.tick_params(axis='x', pad=15,labelsize=20)


    axH.tick_params(axis='y', labelsize=20)
#    axB.yaxis.label.set_color('red')  
    axH.tick_params(axis='x', pad=15,labelsize= 20)

    axI.tick_params(axis='y', labelsize=20)
#    axC.yaxis.label.set_color('red')  
    axI.tick_params(axis='x', pad=15,labelsize= 20)


    #tight layout
    fig1.tight_layout()
    fig2.tight_layout()
    fig3.tight_layout()
    fig4.tight_layout()
    fig5.tight_layout()
    fig6.tight_layout() 
    fig7.tight_layout()
    fig8.tight_layout()
    fig9.tight_layout()   
    fig10.tight_layout()
    fig11.tight_layout()
    fig12.tight_layout()

    


    
#---------------Finalize Figure---------------#    

    ##---Save Figure--##
    fig1.savefig('NUE_TimeSeries_DM.png') 
    fig2.savefig('NUE_TimeSeries_MM.png') 
    fig3.savefig('NUE_TimeSeries_WM.png') 
    fig4.savefig('WUE_TimeSeries_DM.png') 
    fig5.savefig('WUE_TimeSeries_MM.png') 
    fig6.savefig('WUE_TimeSeries_WM.png')     
    fig7.savefig('A_TimeSeries_DM.png') 
    fig8.savefig('A_TimeSeries_MM.png') 
    fig9.savefig('A_TimeSeries_WM.png')     
    fig10.savefig('Temp_TimeSeries_DM.png') 
    fig11.savefig('vwc_TimeSeries_MM.png') 
    fig12.savefig('na_TimeSeries_WM.png')        