#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May  3 13:33:48 2017

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
from mpl_toolkits.mplot3d import Axes3D
#from matplotlib import rcParams

#Import combinations of variable parameters 
from leaf_params_sensitivity_analysis import leaf_params

#Import photosynthesis model
from Photosynthesis_Model import photo_bound_meso_eqstom as photo

#Import functions to switch between Pa and umol/mol at sea level
from photo_functions import pa_con_atmfrac


#---------------Determine if I Want to Keep Any of the Variable Parameters Constant---------------#

var_params=[]
for xy in it.combinations(['na','chl','dia','sla','ht'],5):
    var_params+=[xy]

#do this when I do not put any of the variable parameters as constant. instead I 
#vary each parameter one at a time while keeping the other parameters constant.
if var_params==[()]:
    var_params=[[-999999]]   

#---------------Begin Looping Through Photosynthesis Model---------------#

#each loop is for a constant value, or combinatin of constant values, of variable parameter as determined above
for ii in range(len(var_params)):
    
    #---------------Initialize Plots---------------#

   
    #--figure 1--#
    
    #put in correct ax value (e.g. axA, axB)
    fig1,axA = plt.subplots(figsize=(15,15))

#    axA.set_xlabel('Volumetric Water Content (m3/m3)',fontsize=36, fontname='Times New Roman')
#    axA.set_xlabel('Surface Temperature (C)',fontsize=36, fontname='Times New Roman')
#    axA.set_xlabel('Total Climate Variation (%)',fontsize=36, fontname='Times New Roman')
#    axA.set_xlabel('Leaf Chlorophyll (umol chl/m2)',fontsize=36, fontname='Times New Roman')
#    axA.set_xlabel('Leaf Diameter (m2)',fontsize=36, fontname='Times New Roman')
#    axA.set_xlabel('Leaf Height (cm)',fontsize=36, fontname='Times New Roman')
#    axA.set_xlabel('Specific Leaf Area (cm2/g)',fontsize=36, fontname='Times New Roman')
#    axA.set_xlabel('Leaf Nitrogen (gN/m2)',fontsize=36, fontname='Times New Roman')
#    axA.set_xlabel('Total Physiological Variation (%)',fontsize=36, fontname='Times New Roman')
    axA.set_xlabel('Climate and Physiological Variation (%)',fontsize=36, fontname='Times New Roman')


    axA.set_ylabel('NUE (umol CO2/g N s)',fontsize=36, fontname='Times New Roman')
    axA.set_ylim([0,4])
    axA.set_title('NUE vs. WUE for Dry Meadow', fontname='Times New Roman',fontsize=36,fontweight='bold')
    #twin axis
    axA2=axA.twinx()
    axA2.set_ylabel('WUE (umol CO2/mmol H2O)', fontsize=36, fontname='Times New Roman')  
    axA2.set_ylim([0,4])


    #--figure 2--#
    fig2,axB = plt.subplots(figsize=(15,15))

#    axB.set_xlabel('Volumetric Water Content (m3/m3)',fontsize=36, fontname='Times New Roman')
#    axB.set_xlabel('Surface Temperature (C)',fontsize=36, fontname='Times New Roman')
#    axB.set_xlabel('Total Climate Variation (%)',fontsize=36, fontname='Times New Roman')
#    axB.set_xlabel('Leaf Chlorophyll (umol chl/m2)',fontsize=36, fontname='Times New Roman')
#    axB.set_xlabel('Leaf Diameter (m2)',fontsize=36, fontname='Times New Roman')
#    axB.set_xlabel('Leaf Height (cm)',fontsize=36, fontname='Times New Roman')
#    axB.set_xlabel('Specific Leaf Area (cm2/g)',fontsize=36, fontname='Times New Roman')
#    axB.set_xlabel('Leaf Nitrogen (gN/m2)',fontsize=36, fontname='Times New Roman')
#    axB.set_xlabel('Total Physiological Variation (%)',fontsize=36, fontname='Times New Roman')
    axB.set_xlabel('Climate and Physiological Variation (%)',fontsize=36, fontname='Times New Roman')


    axB.set_ylabel('NUE (umol CO2/g N s)',fontsize=36, fontname='Times New Roman')
    axB.set_ylim([0,4])
    axB.set_title('NUE vs. WUE for Moist Meadow', fontname='Times New Roman',fontsize=36,fontweight='bold')
    #twin axis
    axB2=axB.twinx()
    axB2.set_ylabel('WUE (umol CO2/mmol H2O)', fontsize=36, fontname='Times New Roman')  
    axB2.set_ylim([0,4])

    #--figure 3--#

    fig3,axC = plt.subplots(figsize=(15,15))

#    axC.set_xlabel('Volumetric Water Content (m3/m3)',fontsize=36, fontname='Times New Roman')
#    axC.set_xlabel('Surface Temperature (C)',fontsize=36, fontname='Times New Roman')
#    axC.set_xlabel('Total Climate Variation (%)',fontsize=36, fontname='Times New Roman')
#    axC.set_xlabel('Leaf Chlorophyll (umol chl/m2)',fontsize=36, fontname='Times New Roman')
#    axC.set_xlabel('Leaf Diameter (m)',fontsize=36, fontname='Times New Roman')
#    axC.set_xlabel('Leaf Height (cm)',fontsize=36, fontname='Times New Roman')
#    axC.set_xlabel('Specific Leaf Area (cm2/g)',fontsize=36, fontname='Times New Roman')
#    axC.set_xlabel('Leaf Nitrogen (gN/m2)',fontsize=36, fontname='Times New Roman')
#    axC.set_xlabel('Total Physiological Variation (%)',fontsize=36, fontname='Times New Roman')
    axC.set_xlabel('Climate and Physiological Variation (%)',fontsize=36, fontname='Times New Roman')
    
    axC.set_ylabel('NUE (umol CO2/g N s)',fontsize=36, fontname='Times New Roman')
    axC.set_ylim([0,4])
    axC.set_title('NUE vs. WUE for Wet Meadow', fontname='Times New Roman',fontsize=36,fontweight='bold')
    #twin axis
    axC2=axC.twinx()
    axC2.set_ylabel('WUE (umol CO2/mmol H2O)', fontsize=36, fontname='Times New Roman')  
    axC2.set_ylim([0,4])    

#3D Plots
   
    fig4 = plt.figure()
    axD = fig4.add_subplot(111, projection='3d') 
    axD.set_title('Dry Meadow NUE As a Function of Climate and Plant Physiology', fontname='Times New Roman',fontsize=16,fontweight='bold')
    axD.set_xlabel('Total Climate Variability (%)',fontsize=12, fontname='Times New Roman')
    axD.set_ylabel('Total Physiological Variability (%)',fontsize=12, fontname='Times New Roman')
    axD.set_zlabel('NUE (umol CO2/g N s)',fontsize=12, fontname='Times New Roman')

    fig5 = plt.figure()
    axE = fig5.add_subplot(111, projection='3d') 
    axE.set_title('Dry Meadow WUE As a Function of Climate and Plant Physiology', fontname='Times New Roman',fontsize=16,fontweight='bold')
    axE.set_xlabel('Total Climate Variability (%)',fontsize=12, fontname='Times New Roman')
    axE.set_ylabel('Total Physiological Variability (%)',fontsize=12, fontname='Times New Roman')
    axE.set_zlabel('WUE (umol CO2/mmol H2O)',fontsize=12, fontname='Times New Roman')

    fig6 = plt.figure()
    axF = fig6.add_subplot(111, projection='3d') 
    axF.set_title('Moist Meadow NUE As a Function of Climate and Plant Physiology', fontname='Times New Roman',fontsize=16,fontweight='bold')
    axF.set_xlabel('Total Climate Variability (%)',fontsize=12, fontname='Times New Roman')
    axF.set_ylabel('Total Physiological Variability (%)',fontsize=12, fontname='Times New Roman')
    axF.set_zlabel('NUE (umol CO2/g N s)',fontsize=12, fontname='Times New Roman')

    fig7 = plt.figure()
    axG = fig7.add_subplot(111, projection='3d') 
    axG.set_title('Moist Meadow WUE As a Function of Climate and Plant Physiology', fontname='Times New Roman',fontsize=16,fontweight='bold')
    axG.set_xlabel('Total Climate Variability (%)',fontsize=12, fontname='Times New Roman')
    axG.set_ylabel('Total Physiological Variability (%)',fontsize=12, fontname='Times New Roman')
    axG.set_zlabel('WUE (umol CO2/mmol H2O)',fontsize=12, fontname='Times New Roman')

    fig8 = plt.figure()
    axH = fig8.add_subplot(111, projection='3d') 
    axH.set_title('Wet Meadow NUE As a Function of Climate and Plant Physiology', fontname='Times New Roman',fontsize=16,fontweight='bold')
    axH.set_xlabel('Total Climate Variability (%)',fontsize=12, fontname='Times New Roman')
    axH.set_ylabel('Total Physiological Variability (%)',fontsize=12, fontname='Times New Roman')
    axH.set_zlabel('NUE (umol CO2/g N s)',fontsize=12, fontname='Times New Roman')

    fig9 = plt.figure()
    axI = fig9.add_subplot(111, projection='3d') 
    axI.set_title('Wet Meadow WUE As a Function of Climate and Plant Physiology', fontname='Times New Roman',fontsize=16,fontweight='bold')
    axI.set_xlabel('Total Climate Variability (%)',fontsize=12, fontname='Times New Roman')
    axI.set_ylabel('Total Physiological Variability (%)',fontsize=12, fontname='Times New Roman')
    axI.set_zlabel('WUE (umol CO2/mmol H2O)',fontsize=12, fontname='Times New Roman')
    
    ####-----------run through various temperatures---------######

   
    #make 30% changes
    
    #abiotic  variability: temp and moisture
#    temps=np.linspace(0.7,1.3,50)
#    vwcs=np.linspace(1.3,0.7,50)

#    #abiotic  variability: temp 
#    temps=np.linspace(0.7,1.3,50)
#    vwcs=np.linspace(1.0,1.0,50)

#    #abiotic  variability: moisture
#    temps=np.linspace(1.0,1.0,50)
#    vwcs=np.linspace(1.3,0.7,50)


#    abiotic and physiological variability (set all traits to varying)
    temps=np.linspace(0.7,1.3,10)
    vwcs=np.linspace(1.3,0.7,10)
    
#    #physiological variability --remember to set plant traits as varying (for total--set all plant traits to varying)
#    temps=np.linspace(1.0,1.0,1)
#    vwcs=np.linspace(1.0,1.0,1)

    abiotic_index=np.linspace(-30,30,len(temps))
#---------------Initialize Arrays for Each Meadow---------------#
    
    #total nue and wue
    nue_tot=[]
    wue_tot=[]
   
    #wue and nue arrays
    wue_f=[]
    nue_f=[]
    wue_d=[]
    nue_d=[]
    wue_m=[]
    nue_m=[]
    wue_w=[]
    nue_w=[]    
    wue_s=[]
    nue_s=[]
    
    #na arrays
    na_f=[]
    na_d=[]
    na_m=[]
    na_w=[]
    na_s=[]
    
    #gsw arrays
    gsw_f=[]
    gsw_d=[]
    gsw_m=[]
    gsw_w=[]
    gsw_s=[]

    #assimilation arrays
    A_f=[]
    A_d=[]
    A_m=[]
    A_w=[]
    A_s=[]    

    #evapo arrays
    E_f=[]
    E_d=[]
    E_m=[]
    E_w=[]
    E_s=[]   
    
    #vapor pressure deficit arrays
    vpd_d=[]
    vpd_m=[]
    vpd_w=[]


    #arrays for x axis
    t_d_x=[]
    vwc_d_x=[]
    ab_i_d_x=[]
 
    t_m_x=[]
    vwc_m_x=[]
    ab_i_m_x=[]  
 
    t_w_x=[]
    vwc_w_x=[]
    ab_i_w_x=[]
    
    #x arrays
    chl_d_x=[]
    dia_d_x=[]
    ht_d_x=[]
    sla_d_x=[]
    na_d_x=[]
    b_i_d_x=[]

    chl_m_x=[]
    dia_m_x=[]
    ht_m_x=[]
    sla_m_x=[]
    na_m_x=[]
    b_i_m_x=[]

    chl_w_x=[]
    dia_w_x=[]
    ht_w_x=[]
    sla_w_x=[]
    na_w_x=[]
    b_i_w_x=[]

 
    color=[]
    for xx in range(len(temps)):
        color+=[[float(xx)/float(len(temps)),0.0,float(len(temps)-xx)/float(len(temps))]]
        temp_per=temps[xx]
        vwc_per=vwcs[xx]
        ab_i=abiotic_index[xx] #abiotic index

        #---------------Photosynthesis + Stomatal Conductance Model---------------#

        
        ##---Constant Parameter Arrays for Model---##
        
        #I have commented out parameters that I am assuming are variable (for the time being)
        
        #-------Params Not Currently Used-----#
        gm=np.zeros(shape=1)+0.84 #(mesophyll conductance)
        vwc=np.zeros(shape=1)+0.15 #Soil Volumetric Water Content (cm3/cm3)
        vwc_min=0.08 #minimum soil water content for photosynthesis to occur (permanent wilting point) (cm3/cm3) 
        vwc_max=0.68 #maximum soil water content where increases in soil water do not affect photosynthesis (field capacity?) (cm3/cm3)
        q=0.2 #parameter for soil water affect on photosynthesis (unitless)
       
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
        flnr=np.zeros(shape=1)+0.2 #fraction of leaf nitrogen in rubisco (g N Rub/g N leaf)
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
        j_m=2.68 #slope coefficient 
        

        #----abiotic parameters that vary for each meadow type----#
        t=[(15+2)*temp_per,(15+0)*temp_per,(15-2)*temp_per] #surface temperature for dry, moist, and wet meadows 
        vwc=[0.12*vwc_per,0.15*vwc_per,0.27*vwc_per] #volumetric water content for dry, moist, and wet meadows (June through August)
        

    
    #---------------Import Variable Parameter Arrays from Leaf Parameter File---------------#
        
        for i in range(len(leaf_params)):
            for key,val in leaf_params[i].items():
                exec(key + '=val')
            
 
            #set variable parameters constant if I specify this above
            if 'na' in var_params[ii]:
                na=np.linspace(na*0.7,na*1.3,10)
            else:
                na=np.linspace(na,na,10)
            
            if 'dia' in var_params[ii]:
                dia=np.linspace(dia*0.7,dia*1.3,10)
            else:
                dia=np.linspace(dia,dia,10)
            
            if 'chl' in var_params[ii]:
                chl=np.linspace(chl*0.7,chl*1.3,10)
            else:
                chl=np.linspace(chl,chl,10)
           
            if 'ht' in var_params[ii]:
                ht=np.linspace(ht*0.7,ht*1.3,10)
            else:
                ht=np.linspace(ht,ht,10)
           
            if 'sla' in var_params[ii]:
                sla=np.linspace(sla*0.7,sla*1.3,10)
            else:
                sla=np.linspace(sla,sla,10)
            


            #---------------Make Array of Values for Each Meadow---------------#  
            
            #biotic index
            biotic_index=np.linspace(-30,30,len(na))
            
            #number of simulations per meadow type:
            m_sim=len(leaf_params)/3.0 #meadow simulations 

            #dry meadow            
            if i<(m_sim):     
                
                for iii in range(len(na)):
                
                    b_i=biotic_index[iii]
                    
                    #correct for microclimate abiotic conditions
                    t_dm=t[0]
    
                    #correct for volumetric moisture content
                    vwc_dm=vwc[0]
                    
                    #------calculate vapor pressure-----#
                    pa_v=611*np.exp((17.27*t_dm)/(t_dm+237.3)) #saturation vapor pressure of air (Pa)
                    ea_str=pa_con_atmfrac(pa_v,3528) #saturation vapor pressure of air (Pa-->umol h20/mol air)
                    ea=rh*ea_str #vapor pressure (umol h2O/mol air)                
    
                    #correct for leaf temperatures using leaf height
           
                    t_diff=18-0.4*ht[iii]
                
                    tl=t_dm+t_diff
                
       
                    #solve for flnr
      
                    sla_corr=sla[iii]/10.0
                    flnr=sla_corr*0.0040+0.0703
    
                
        
                    #---------------Photosynthesis Function---------------#
                
                    #alter this line of code for when implementing different photosynthesis functions
                    wue, nue, A, E, cs, ci, gsw, gs, gbw, gb, gm, cc,dd =photo(tk_25,ekc,eko,etau,ev,ej,toptv,toptj,na[iii], qeff, PAR,tl,ea,chl[iii],ij,kc25,ko25,o,ca,rh,m,a,frnr,flnr,ra,j_m,g0,b,dia[iii],u,q,vwc_min,vwc_max,vwc_dm)     
                   
                    #test to make sure wue and nue are positive at not 'nan'
                    if wue[0]==-999 and nue[0]==-999:
                   
                        continue
            
 
                    wue_d+=[wue[0]]
                    nue_d+=[nue[0]]
                    na_d+=[na[0]]
                    gsw_d+=[gsw[0]]
                    A_d+=[A[0]]
                    E_d+=[E[0]]
                    vpd_d+=[dd[0]]
                    t_d_x+=[t_dm]
                    vwc_d_x+=[vwc_dm]
                    chl_d_x+=[chl[iii]]
                    dia_d_x+=[dia[iii]]
                    ht_d_x+=[ht[iii]]
                    sla_d_x+=[sla[iii]]
                    na_d_x+=[na[iii]]
                    ab_i_d_x+=[ab_i]
                    b_i_d_x+=[b_i]
                
                
                
            #moist meadow       
            elif i>=(m_sim) and i<(m_sim*2):
                
                for iii in range(len(na)):

                    b_i=biotic_index[iii]
                    
                    #correct for microclimate abiotic conditions
                    t_mm=t[1]
    
                    #correct for volumetric moisture content
                    vwc_mm=vwc[1]
                    
                    #------calculate vapor pressure-----#
                    pa_v=611*np.exp((17.27*t_mm)/(t_mm+237.3)) #saturation vapor pressure of air (Pa)
                    ea_str=pa_con_atmfrac(pa_v,3528) #saturation vapor pressure of air (Pa-->umol h20/mol air)
                    ea=rh*ea_str #vapor pressure (umol h2O/mol air)                
    
    
                    #correct for leaf temperatures using leaf height
           
                    t_diff=18-0.4*ht[iii]
                
                    tl=t_mm+t_diff
                
       
                    #solve for flnr
      
                    sla_corr=sla[iii]/10.0
                    flnr=sla_corr*0.0040+0.0703
    
                
        
                    #---------------Photosynthesis Function---------------#
                
                    #alter this line of code for when implementing different photosynthesis functions
                    wue, nue, A, E, cs, ci, gsw, gs, gbw, gb, gm, cc,dd =photo(tk_25,ekc,eko,etau,ev,ej,toptv,toptj,na[iii], qeff, PAR,tl,ea,chl[iii],ij,kc25,ko25,o,ca,rh,m,a,frnr,flnr,ra,j_m,g0,b,dia[iii],u,q,vwc_min,vwc_max,vwc_mm)
          
               
                    #test to make sure wue and nue are positive at not 'nan'
                    if wue[0]==-999 and nue[0]==-999:
                   
                        continue
                
                    wue_m+=[wue[0]]
                    nue_m+=[nue[0]]
                    na_m+=[na[0]]
                    gsw_m+=[gsw[0]]
                    A_m+=[A[0]]
                    E_m+=[E[0]]
                    vpd_m+=[dd[0]]
   
                    t_m_x+=[t_mm]
                    vwc_m_x+=[vwc_mm]
                    chl_m_x+=[chl[iii]]
                    dia_m_x+=[dia[iii]]
                    ht_m_x+=[ht[iii]]
                    sla_m_x+=[sla[iii]]
                    na_m_x+=[na[iii]]
                    ab_i_m_x+=[ab_i]
                    b_i_m_x+=[b_i]
        


            #wet meadow
            elif i>=(m_sim*2) and i<(m_sim*3):
                
                for iii in range(len(na)):

                    b_i=biotic_index[iii]

                    #correct for microclimate abiotic conditions
                    t_wm=t[2]
                    
                    #correct for volumetric moisture content
                    vwc_wm=vwc[2]
    
                    #------calculate vapor pressure-----#
                    pa_v=611*np.exp((17.27*t_wm)/(t_wm+237.3)) #saturation vapor pressure of air (Pa)
                    ea_str=pa_con_atmfrac(pa_v,3528) #saturation vapor pressure of air (Pa-->umol h20/mol air)
                    ea=rh*ea_str #vapor pressure (umol h2O/mol air)                
    
    
                    #correct for leaf temperatures using leaf height
           
                    t_diff=18-0.4*ht[iii]
                
                    tl=t_wm+t_diff
                
       
                    #solve for flnr
      
                    sla_corr=sla[iii]/10.0
                    flnr=sla_corr*0.0040+0.0703
    
                
        
                    #---------------Photosynthesis Function---------------#
                
                    #alter this line of code for when implementing different photosynthesis functions
                    wue, nue, A, E, cs, ci, gsw, gs, gbw, gb, gm, cc,dd =photo(tk_25,ekc,eko,etau,ev,ej,toptv,toptj,na[iii], qeff, PAR,tl,ea,chl[iii],ij,kc25,ko25,o,ca,rh,m,a,frnr,flnr,ra,j_m,g0,b,dia[iii],u,q,vwc_min,vwc_max,vwc_wm)
      
               
                    #test to make sure wue and nue are positive at not 'nan'
                    if wue[0]==-999 and nue[0]==-999:
                   
                        continue                
                
                    wue_w+=[wue[0]]
                    nue_w+=[nue[0]]
                    na_w+=[na[0]]
                    gsw_w+=[gsw[0]]
                    A_w+=[A[0]]
                    E_w+=[E[0]]
                    vpd_w+=[dd[0]]
                    t_w_x+=[t_wm]
                    vwc_w_x+=[vwc_wm]
                    chl_w_x+=[chl[iii]]
                    dia_w_x+=[dia[iii]]
                    ht_w_x+=[ht[iii]]
                    sla_w_x+=[sla[iii]]
                    na_w_x+=[na[iii]]
                    ab_i_w_x+=[ab_i]
                    b_i_w_x+=[b_i]
    
#---------------Regression Plot WUE vs. NUE---------------#   
   

    axA.plot(ab_i_d_x,nue_d,color='red',linewidth=10,alpha=0.5)
#    axA.scatter(t_d_x,nue_d,color='red',facecolors='red',marker='o',s=50,alpha=0.5)
    axA2.plot(ab_i_d_x,wue_d,color='blue',linewidth=10,alpha=0.5)
#    axA2.scatter(t_d_x,wue_d,color='blue',facecolors='blue',marker='o',s=50,alpha=0.5)

    axB.plot(ab_i_m_x,nue_m,color='red',linewidth=10,alpha=0.5)
#    axB.scatter(t_m_x,nue_m,color='red',facecolors='red',marker='o',s=50,alpha=0.5)
    axB2.plot(ab_i_m_x,wue_m,color='blue',linewidth=10,alpha=0.5)
#    axB2.scatter(t_m_x,wue_m,color='blue',facecolors='blue',marker='o',s=50,alpha=0.5)
 
    axC.plot(ab_i_w_x,nue_w,color='red',linewidth=10,alpha=0.5)
#    axC.scatter(t_w_x,nue_w,color='red',facecolors='red',marker='o',s=50,alpha=0.5)
    axC2.plot(ab_i_w_x,wue_w,color='blue',linewidth=10,alpha=0.5)
#    axC2.scatter(t_w_x,wue_w,color='blue',facecolors='blue',marker='o',s=50,alpha=0.5)
#   

#3D Plots

    axD.plot_trisurf(ab_i_d_x,b_i_d_x,nue_d,linewidth=0.2, antialiased=False,color="red",alpha=0.5)
    axE.plot_trisurf(ab_i_d_x,b_i_d_x,wue_d,linewidth=0.2, antialiased=False,color="blue",alpha=0.5)

    axF.plot_trisurf(ab_i_m_x,b_i_m_x,nue_m,linewidth=0.2, antialiased=False,color="red",alpha=0.5)
    axG.plot_trisurf(ab_i_m_x,b_i_m_x,wue_m,linewidth=0.2, antialiased=False,color="blue",alpha=0.5)

    axH.plot_trisurf(ab_i_w_x,b_i_w_x,nue_w,linewidth=0.2, antialiased=False,color="red",alpha=0.5)
    axI.plot_trisurf(ab_i_w_x,b_i_w_x,wue_w,linewidth=0.2, antialiased=False,color="blue",alpha=0.5)

    
    axA.tick_params(axis='y', labelsize=16)
    axB.tick_params(axis='y', labelsize=16)
    axC.tick_params(axis='y', labelsize=16)
    axA2.tick_params(axis='y', labelsize=16)
    axB2.tick_params(axis='y', labelsize=16)
    axC2.tick_params(axis='y', labelsize=16)
    axA.tick_params(axis='x', labelsize=16)
    axB.tick_params(axis='x', labelsize=16)
    axC.tick_params(axis='x', labelsize=16)
             
    axD.tick_params(axis='x', labelsize=8)
    axD.tick_params(axis='y', labelsize=8)
    axD.tick_params(axis='z', labelsize=8)
    
    axE.tick_params(axis='x', labelsize=8)
    axE.tick_params(axis='y', labelsize=8)
    axE.tick_params(axis='z', labelsize=8)
    
    axF.tick_params(axis='x', labelsize=8)
    axF.tick_params(axis='y', labelsize=8)
    axF.tick_params(axis='z', labelsize=8)
    
    axG.tick_params(axis='x', labelsize=8)
    axG.tick_params(axis='y', labelsize=8)
    axG.tick_params(axis='z', labelsize=8)
    
    axH.tick_params(axis='x', labelsize=8)
    axH.tick_params(axis='y', labelsize=8)
    axH.tick_params(axis='z', labelsize=8)
    
    axI.tick_params(axis='x', labelsize=8)
    axI.tick_params(axis='y', labelsize=8)
    axI.tick_params(axis='z', labelsize=8)
      

    
#---------------Finalize Figure---------------#    

    ##---Legend---##
    #axA.legend(bbox_to_anchor=(1, 1), loc='left', prop={'size':15})

    ##---Save Figure--##
#    fig1.savefig('NUE_vs_WUE_dm_a_biotic.png')
#    fig2.savefig('NUE_vs_WUE_mm_a_biotic.png')
#    fig3.savefig('NUE_vs_WUE_wm_a_biotic.png')
    fig4.savefig('NUE_dm_a_biotic')
    fig5.savefig('WUE_dm_a_biotic')
    fig6.savefig('NUE_mm_a_biotic')
    fig7.savefig('WUE_mm_a_biotic')
    fig8.savefig('NUE_wm_a_biotic')
    fig9.savefig('WUE_wm_a_biotic')    