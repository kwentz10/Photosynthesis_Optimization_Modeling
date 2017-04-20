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
#from matplotlib import rcParams

#Import combinations of variable parameters 
from leaf_parameter_inputs import leaf_params

#Import photosynthesis model
from Photosynthesis_Model import photo_bound_meso_eqstom as photo

#Import functions to switch between Pa and umol/mol at sea level
from photo_functions import pa_con_atmfrac


#---------------Determine if I Want to Keep Any of the Variable Parameters Constant---------------#

const_params=[]
for xx in it.combinations(['dia','na','ht'],0):
    const_params+=[xx]

#do this when I do not put any of the variable parameters as constant. instead I 
#vary each parameter one at a time while keeping the other parameters constant.
if const_params==[()]:
    const_params=[[-999999]]   

#---------------Begin Looping Through Photosynthesis Model---------------#

#each loop is for a constant value, or combinatin of constant values, of variable parameter as determined above
for ii in range(len(const_params)):
    
    #---------------Initialize Plots---------------#

    ##---Figure With Subplots Blueprint---##

    #fb1=plt.figure(1,figsize=(12,2)) 
    #axA = fb1.add_subplot(121)
    #axB = fb1.add_subplot(122)

    ##---Figures Without Subplots Blueprint---##
    
    #--figure 1--#
    
    #put in correct ax value (e.g. axA, axB)
    fig1,axA = plt.subplots(figsize=(30,15))
    
    #twin axis
    axA2=axA.twinx()

    #axA.set_xlabel('Plant Communities',fontsize=20, fontname='Times New Roman')
    axA.set_ylabel('NUE (umol CO2/g N s)',fontsize=36, fontname='Times New Roman')
    axA2.set_ylabel('WUE (umol CO2/mmol H2O)', fontsize=36, fontname='Times New Roman')  
#    axA.set_xlim([0,30])
    axA.set_ylim([0,7])
    axA2.set_ylim([0,1.5])
    axA.set_title('WUE and NUE for Alpine Tundra Plant Communities', fontname='Times New Roman',fontsize=36,fontweight='bold')

    
    #--figure 2--#
    
    #put in correct ax value (e.g. axA, axB)
    fig2,axB = plt.subplots(figsize=(15,15))

    #axA.set_xlabel('Plant Communities',fontsize=20, fontname='Times New Roman')
    axB.set_ylabel('gsw (mol H2O/m2s)',fontsize=36, fontname='Times New Roman')
    #axA.set_xlim([0,20])
    #axA.set_ylim([0,10])
    axB.set_title('Plant Communities vs. Stomatal Conductance', fontname='Times New Roman',fontsize=36,fontweight='bold')

    #--figure 3--#
    
    #put in correct ax value (e.g. axA, axB)
    fig3,axC = plt.subplots(figsize=(15,15))

    #axA.set_xlabel('Plant Communities',fontsize=20, fontname='Times New Roman')
    axC.set_ylabel('Leaf N (g N/m2)',fontsize=36, fontname='Times New Roman')
    #axA.set_xlim([0,20])
    #axA.set_ylim([0,10])
    axC.set_title('Plant Communities vs. Leaf Nitrogen', fontname='Times New Roman',fontsize=36,fontweight='bold')

    #--figure 4--#
    
    #put in correct ax value (e.g. axA, axB)
    fig4,axD = plt.subplots(figsize=(15,15))

    #axA.set_xlabel('Plant Communities',fontsize=20, fontname='Times New Roman')
    axD.set_ylabel('Assimilation (umol CO2/m2s)',fontsize=36, fontname='Times New Roman')
    #axA.set_xlim([0,20])
    #axA.set_ylim([0,10])
    axD.set_title('Plant Communities vs. Assimilation', fontname='Times New Roman',fontsize=36,fontweight='bold')

  
    #--figure 5--#
    
    #put in correct ax value (e.g. axA, axB)
    fig5,axE = plt.subplots(figsize=(15,15))

    #axA.set_xlabel('Plant Communities',fontsize=20, fontname='Times New Roman')
    axE.set_ylabel('Evapotranspiration (mmol H2O/m2s)',fontsize=36, fontname='Times New Roman')
    #axA.set_xlim([0,20])
    #axA.set_ylim([0,10])
    axE.set_title('Plant Communities vs. Evapotranspiration', fontname='Times New Roman',fontsize=36,fontweight='bold')

   
    #--figure 6--#
    
    #put in correct ax value (e.g. axA, axB)
    fig6,axF = plt.subplots(figsize=(15,15))

    axF.set_xlabel('NUE (umol CO2/g N s)',fontsize=36, fontname='Times New Roman')
    axF.set_ylabel('WUE (umol CO2/mmol H2O)',fontsize=36, fontname='Times New Roman')
    axF.set_xlim([0,7])
    axF.set_ylim([0,1.5])
    axF.set_title('NUE vs. WUE for all Plant Communities', fontname='Times New Roman',fontsize=36,fontweight='bold')

    
    ####-----------run through various temperatures---------######
    
    
#    temps=[8.0,10.0,12.0,14.0,16.0,18.0,20.0,22.0,24.0,26.0]
#    lw=[1,2,3,4,5,6,7,8,9,10]
    temps=[12.0]
    lw=[1]
#    temps=[8.0,26.0]
#    lw=[1,5]

    
    color=[]
    for xx in range(len(temps)):
        color+=[[float(xx)/float(len(temps)),0.0,float(len(temps)-xx)/float(len(temps))]]

    for iii in range(len(temps)):
        t=temps[iii]

        #---------------Photosynthesis + Stomatal Conductance Model---------------#

        
        ##---Constant Parameter Arrays for Model---##
        
        #I have commented out parameters that I am assuming are variable (for the time being)
        
        #-------Params Not Currently Used-----#
        gm=np.zeros(shape=1)+0.84 #(mesophyll conductance)
        vwc=np.zeros(shape=1)+0.15 #Soil Volumetric Water Content (cm3/cm3)
        vwc_min=0.08 #minimum soil water content for photosynthesis to occur (permanent wilting point) (cm3/cm3) 
        vwc_max=0.5 #maximum soil water content where increases in soil water do not affect photosynthesis (field capacity?) (cm3/cm3)
        q=0.1 #parameter for soil water affect on photosynthesis (unitless)
       
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
        flnr=np.zeros(shape=1)+0.13 #fraction of leaf nitrogen in rubisco (g N Rub/g N leaf)
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
        
        #------constant variable params for sensitivty analysis-----#
        
        chl_c=np.zeros(shape=1)+(np.mean([396,465,476])) #Chlorophyll Content of the Leaf (umol chl/m2)
        ht_c=np.zeros(shape=1)+(np.mean([9.2,19.5,20.0])) #Temperature of the Leaf (K)
        dia_c=np.zeros(shape=1)+(np.mean([1.4,2.3,2.6])/100.) #Mean diameter or size of leaf (m)
        na_c=np.zeros(shape=1)+(np.mean([2.5,5.6,6.3])) #leaf nitrogen (g N/ m2)
        
        #------preliminary calculations-----#
        pa_v=611*np.exp((17.27*t)/(t+237.3)) #saturation vapor pressure of air (Pa)
        ea_str=pa_con_atmfrac(pa_v,3528) #saturation vapor pressure of air (Pa-->umol h20/mol air)
        ea=rh*ea_str #vapor pressure deficit (umol h2O/mol air)
        
#don't initialize plots everytime I run this!!            
    
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
    
    #---------------Import Variable Parameter Arrays from Leaf Parameter File---------------#
        
        for i in range(len(leaf_params)):
            for key,val in leaf_params[i].items():
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

    
            
    #correct for leaf temperatures using leaf height
#            if ht[0]>=1.0 and ht[0]<=5.0:
#                tl=25+t
#            elif ht[0]>=6.0 and ht[0]<=10.0:
#                tl=20+t
#            elif ht[0]>=11.0 and ht[0]<=15.0:
#                tl=15+t
#            elif ht[0]>=16.0 and ht[0]<=20.0:
#                tl=10+t
#            elif ht[0]>=21.0 and ht[0]<=25.0:
#                tl=5+t
       
            t_diff=25.0-0.5*ht
            
            tl=t+t_diff
#    
   
 
            
            
                
    #---------------Photosynthesis Function---------------#
        
    
            
            #alter this line of code for when implementing different photosynthesis functions
            wue, nue, A, E, na, cs, ci, gsw, gs, gbw, gb, gm, cc =photo(tk_25,ekc,eko,etau,ev,ej,toptv,toptj,na, qeff, PAR,tl,ea,chl,ij,kc25,ko25,o,ca,rh,m,a,frnr,flnr,ra,j_m,g0,b,dia,u,q,vwc_min,vwc_max,vwc)
 
           
            #test to make sure wue and nue are positive at not 'nan'
            if wue[0]==-999 and nue[0]==-999:
               
                continue
        
    #---------------Make Array of Values for Each Meadow---------------#  
            
            #number of simulations per meadow type:
            m_sim=len(leaf_params)/3.0 #meadow simulations 
            
            if i<(m_sim):  
                wue_d+=[wue]
                nue_d+=[nue]
                na_d+=[na]
                gsw_d+=[gsw]
                A_d+=[A]
                E_d+=[E]
            if i>=(m_sim) and i<(m_sim*2):
                wue_m+=[wue]
                nue_m+=[nue]
                na_m+=[na]
                gsw_m+=[gsw]
                A_m+=[A]
                E_m+=[E]

            elif i>=(m_sim*2) and i<(m_sim*3):
                wue_w+=[wue]
                nue_w+=[nue]
                na_w+=[na]
                gsw_w+=[gsw]
                A_w+=[A]
                E_w+=[E]
              
            
           
                
#            elif i>=(m_sim*3) and i<(m_sim*4):
#                wue_w+=[wue]
#                nue_w+=[nue]
#                na_w+=[na]
#                gsw_w+=[gsw]
#                A_w+=[A]
#                E_w+=[E]
#            elif i>=(m_sim*4) and i<(m_sim*5):
#                wue_s+=[wue]
#                nue_s+=[nue]
#                na_s+=[na]
#                gsw_s+=[gsw]
#                A_s+=[A]
#                E_s+=[E]
            
            nue_tot+=[nue]
            wue_tot+=[wue]
    
    #---------------Plot Plant Communities vs. NUE & WUE---------------#          
        
    #    nue_bp=axA.boxplot([nue_f,nue_d,nue_m,nue_w,nue_s], positions=[0.875,1.875,2.875,3.875,4.875],widths=0.25, patch_artist=True, showmeans=True, showfliers=False)
    #    wue_bp=axA2.boxplot([wue_f,wue_d,wue_m,wue_w,wue_s], positions=[1.125,2.125,3.125,4.125,5.125],widths=0.25, patch_artist=True, showmeans=True, showfliers=False)
    #    axA.plot([0.875,1.875,2.875,3.875,4.875],[np.mean(nue_f),np.mean(nue_d),np.mean(nue_m),np.mean(nue_w),np.mean(nue_s)],'-r')
    #    axA2.plot([1.125,2.125,3.125,4.125,5.125],[np.mean(wue_f),np.mean(wue_d),np.mean(wue_m),np.mean(wue_w),np.mean(wue_s)],'-b')
    
        axA.bar([0.875,1.875,2.875],[np.mean(nue_d),np.mean(nue_m),np.mean(nue_w)], yerr=[np.std(nue_d),np.std(nue_m),np.std(nue_w)], edgecolor='black',linewidth=lw[iii], width=0.2, color='red',alpha=1.0,error_kw={'ecolor':'black', 'lw':2, 'capsize':5, 'capthick':2})
        axA2.bar([1.125,2.125,3.125],[np.mean(wue_d),np.mean(wue_m),np.mean(wue_w)], yerr=[np.std(wue_d),np.std(wue_m),np.std(wue_w)], edgecolor='black',linewidth=lw[iii], width=0.2, color='blue',alpha=1.0,error_kw={'ecolor':'black', 'lw':2, 'capsize':5, 'capthick':2})
        
        axA.set_xticks([1, 2, 3])
        axA.set_xticklabels(['Dry Meadow','Moist Meadow','Wet Meadow'],fontname='Times New Roman')
        axA.tick_params(axis='x', labelsize=28)
        axA.tick_params(axis='y', labelsize=18)
        axA2.tick_params(axis='y', labelsize=18)
        
    #    #nue boxplot specs
    #    for box in nue_bp['boxes']:
    #        #change outline color
    #        box.set(color='red',linewidth=2)
    #        #change fill color
    #        box.set(facecolor='red',alpha=0.2)
    #
    #    for whisker in nue_bp['whiskers']:
    #        whisker.set(color='red',linewidth=2)
    #    
    #    for cap in nue_bp['caps']:
    #        cap.set(color='red',linewidth=2)
    #
    #    for median in nue_bp['medians']:
    #        median.set(color='red', linewidth=2)
    #
    #    for flier in nue_bp['fliers']:
    #        flier.set(marker='*',color='red',alpha=0.5)
    #
    #    for means in nue_bp['means']:
    #        means.set(marker='o',markerfacecolor='red')    
    #
    #    
    #    #wue boxplot specs    
    #    for box in wue_bp['boxes']:
    #        #change outline color
    #        box.set(color='blue',linewidth=2)
    #        #change fill color
    #        box.set(facecolor='blue',alpha=0.2)
    #
    #    for whisker in wue_bp['whiskers']:
    #        whisker.set(color='blue',linewidth=2)
    #    
    #    for cap in wue_bp['caps']:
    #        cap.set(color='blue',linewidth=2)
    #
    #    for median in wue_bp['medians']:
    #        median.set(color='blue', linewidth=2)
    #
    #    for flier in wue_bp['fliers']:
    #        flier.set(marker='*',color='blue',alpha=0.5)
    #
    #    for means in wue_bp['means']:
    #        means.set(marker='o',markerfacecolor='blue')
              
    
    #---------------Box Plot Plant Communities vs. Stomatal Condcutance--------------- #     
    
    #    gsw_bp=axB.boxplot([gsw_f,gsw_d,gsw_m,gsw_w,gsw_s], patch_artist=True, showmeans=True, showfliers=False)
    #    axB.plot([1,2,3,4,5],[np.mean(gsw_f),np.mean(gsw_d),np.mean(gsw_m),np.mean(gsw_w),np.mean(gsw_s)],'-c')
        
        
        axB.bar([1,2,3],[np.mean(gsw_d),np.mean(gsw_m),np.mean(gsw_w)], yerr=[np.std(gsw_d),np.std(gsw_m),np.std(gsw_w)],edgecolor='black',linewidth=lw[iii], width=0.2, color='cyan',alpha=1.0,error_kw={'ecolor':'black', 'lw':2, 'capsize':5, 'capthick':2})
        axB.set_xticks([1, 2, 3])
        axB.set_xticklabels(['Dry Meadow','Moist Meadow','Wet Meadow'],fontname='Times New Roman')
        axB.tick_params(axis='x', labelsize=25)
        axB.tick_params(axis='y', labelsize=16)
        
        #gsw boxplot specs
    #    for box in gsw_bp['boxes']:
    #        #change outline color
    #        box.set(color='cyan',linewidth=2)
    #        #change fill color
    #        box.set(facecolor='cyan',alpha=0.5)
    #
    #    for whisker in gsw_bp['whiskers']:
    #        whisker.set(color='cyan',linewidth=2)
    #    
    #    for cap in gsw_bp['caps']:
    #        cap.set(color='cyan',linewidth=2)
    #
    #    for median in gsw_bp['medians']:
    #        median.set(color='cyan', linewidth=2)
    #
    #    for flier in gsw_bp['fliers']:
    #        flier.set(marker='*',color='cyan',alpha=0.5)
    #
    #    for means in gsw_bp['means']:
    #        means.set(marker='o',markerfacecolor='cyan')    
    
    #---------------Box Plot Plant Communities vs. Leaf N--------------- #     
    
    #    na_bp=axC.boxplot([na_f,na_d,na_m,na_w,na_s], patch_artist=True, showmeans=True, showfliers=False)
    #    axC.plot([1,2,3,4,5],[np.mean(na_f),np.mean(na_d),np.mean(na_m),np.mean(na_w),np.mean(na_s)],'-',color='orange')
        
        axC.bar([1,2,3],[np.mean(na_d),np.mean(na_m),np.mean(na_w)], yerr=[np.std(na_d),np.std(na_m),np.std(na_w)],edgecolor='black',linewidth=lw[iii], width=0.2, color='orange',alpha=1.0,error_kw={'ecolor':'black', 'lw':2, 'capsize':5, 'capthick':2})
        axC.set_xticks([1, 2, 3])
        axC.set_xticklabels(['Dry Meadow','Moist Meadow','Wet Meadow'],fontname='Times New Roman')
        axC.tick_params(axis='x', labelsize=25)
        axC.tick_params(axis='y', labelsize=16)
        
    #    #na boxplot specs
    #    for box in na_bp['boxes']:
    #        #change outline color
    #        box.set(color='orange',linewidth=2)
    #        #change fill color
    #        box.set(facecolor='orange',alpha=0.5)
    #
    #    for whisker in na_bp['whiskers']:
    #        whisker.set(color='orange',linewidth=2)
    #    
    #    for cap in na_bp['caps']:
    #        cap.set(color='orange',linewidth=2)
    #
    #    for median in na_bp['medians']:
    #        median.set(color='orange', linewidth=2)
    #
    #    for flier in na_bp['fliers']:
    #        flier.set(marker='*',color='orange',alpha=0.5)
    #
    #    for means in na_bp['means']:
    #        means.set(marker='o',markerfacecolor='orange')    
    
            
    #---------------Box Plot Plant Communities vs. Assimilation--------------- #     
    
    #    A_bp=axD.boxplot([A_f,A_d,A_m,A_w,A_s], patch_artist=True, showmeans=True, showfliers=False)
    #    axD.plot([1,2,3,4,5],[np.mean(A_f),np.mean(A_d),np.mean(A_m),np.mean(A_w),np.mean(A_s)],'-',color='purple')
        
        axD.bar([1,2,3],[np.mean(A_d),np.mean(A_m),np.mean(A_w)],yerr=[np.std(A_d),np.std(A_m),np.std(A_w)],edgecolor='black',linewidth=lw[iii], width=0.2,color='purple',alpha=1.0,error_kw={'ecolor':'black', 'lw':2, 'capsize':5, 'capthick':2})
        axD.set_xticks([1, 2, 3])
        axD.set_xticklabels(['Dry Meadow','Moist Meadow','Wet Meadow'],fontname='Times New Roman')
        axD.tick_params(axis='x', labelsize=25)
        axD.tick_params(axis='y', labelsize=16)
        
    #    #A boxplot specs
    #    for box in A_bp['boxes']:
    #        #change outline color
    #        box.set(color='purple',linewidth=2)
    #        #change fill color
    #        box.set(facecolor='purple',alpha=0.5)
    #
    #    for whisker in A_bp['whiskers']:
    #        whisker.set(color='purple',linewidth=2)
    #    
    #    for cap in A_bp['caps']:
    #        cap.set(color='purple',linewidth=2)
    #
    #    for median in A_bp['medians']:
    #        median.set(color='purple', linewidth=2)
    #
    #    for flier in A_bp['fliers']:
    #        flier.set(marker='*',color='purple',alpha=0.5)
    #
    #    for means in A_bp['means']:
    #        means.set(marker='o',markerfacecolor='purple')  
            
    #---------------Box Plot Plant Communities vs. Evapotranspiration--------------- #     
    
    #    E_bp=axE.boxplot([E_f,E_d,E_m,E_w,E_s], patch_artist=True, showmeans=True, showfliers=False)
    #    axE.plot([1,2,3,4,5],[np.mean(E_f),np.mean(E_d),np.mean(E_m),np.mean(E_w),np.mean(E_s)],'-k')
        
        axE.bar([1,2,3],[np.mean(E_d)/1000.0,np.mean(E_m)/1000.0,np.mean(E_w)/1000.0],yerr=[np.std(E_d)/1000.0,np.std(E_m)/1000.0,np.std(E_w)/1000.0],edgecolor='black',linewidth=lw[iii], width=0.2,color='green',alpha=1.0,error_kw={'ecolor':'black', 'lw':2, 'capsize':5, 'capthick':2})
        axE.set_xticks([1, 2, 3])
        axE.set_xticklabels(['Dry Meadow','Moist Meadow','Wet Meadow'],fontname='Times New Roman')
        axE.tick_params(axis='x', labelsize=25)
        axE.tick_params(axis='y', labelsize=16)
        
        
    #    #E boxplot specs
    #    for box in E_bp['boxes']:
    #        #change outline color
    #        box.set(color='black',linewidth=2)
    #        #change fill color
    #        box.set(facecolor='black',alpha=0.5)
    #
    #    for whisker in E_bp['whiskers']:
    #        whisker.set(color='black',linewidth=2)
    #    
    #    for cap in E_bp['caps']:
    #        cap.set(color='black',linewidth=2)
    #
    #    for median in E_bp['medians']:
    #        median.set(color='black', linewidth=2)
    #
    #    for flier in E_bp['fliers']:
    #        flier.set(marker='*',color='black',alpha=0.5)
    #
    #    for means in E_bp['means']:
    #        means.set(marker='o',markerfacecolor='black')  
            
            
    #---------------Regression Plot WUE vs. NUE---------------#   
   
        
        axF.scatter([np.mean(nue_d),np.mean(nue_m),np.mean(nue_w)],[np.mean(wue_d),np.mean(wue_m),np.mean(wue_w)],edgecolors=[color[iii]]*3,facecolors=[color[iii]]*3,marker='o',s=50)
        axF.plot([np.mean(nue_d),np.mean(nue_m),np.mean(nue_w)],[np.mean(wue_d),np.mean(wue_m),np.mean(wue_w)],color=color[iii])
        axF.arrow(np.mean(nue_m),np.mean(wue_m),np.mean(nue_w)-np.mean(nue_m),np.mean(wue_w)-np.mean(wue_m), head_width=0.15, head_length=0.15, fc=color[iii], ec=color[iii])
    #    axF.scatter(nue_tot,wue_tot) 
        
        #axF.plot(np.unique([np.mean(nue_d),np.mean(nue_m),np.mean(nue_w)]), np.poly1d(np.polyfit([np.mean(nue_d),np.mean(nue_m),np.mean(nue_w)], [np.mean(wue_d),np.mean(wue_m),np.mean(wue_w)], 1))(np.unique([np.mean(nue_d),np.mean(nue_m),np.mean(nue_w)])),color=color[iii])
        #axF.plot([np.mean(nue_f),np.mean(nue_d),np.mean(nue_m),np.mean(nue_w),np.mean(nue_s)],[np.mean(wue_f),np.mean(wue_d),np.mean(wue_m),np.mean(wue_w),np.mean(wue_s)])
    
        #axF.tick_params(axis='x', labelsize=25)
        axF.tick_params(axis='y', labelsize=16)
        
    #---------------Make Plot Interactive---------------# 
       
        plt.pause(0.00000001)
        plt.ion()

    
#---------------Finalize Figure---------------#    

    ##---Legend---##
    #axA.legend(bbox_to_anchor=(1, 1), loc='left', prop={'size':15})

    ##---Save Figure--##
#    fig1.savefig('WUE_vs_NUE_var_chl.png') 
#    fig2.savefig('gsw_plantcoms.png') 
#    fig3.savefig('na_plantcoms.png') 
#    fig4.savefig('A_plantcoms.png')
#    fig5.savefig('E_plantcoms.png')
#    fig6.savefig('WUE_vs_NUE_regression_var_chl.png')
