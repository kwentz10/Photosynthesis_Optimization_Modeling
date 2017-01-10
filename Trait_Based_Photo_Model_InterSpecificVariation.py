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

"""

#---------------Import Modules---------------#

import itertools as it
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rcParams

#Import combinations of variable parameters 
from leaf_parameter_inputs import leaf_params

#Import photosynthesis model
from Photosynthesis_Model import photo_bound_meso_eqstom as photo

#Import functions to switch between Pa and umol/mol at sea level
from photo_functions import pa_con_atmfrac


#---------------Photosynthesis + Stomatal Conductance Model---------------#

##---Intercept of Carboxylation Rate vs. Light Absorption Rate---##
j_b=20.9 #intercept coefficient 

##---Maximum Slope of Carboxylation Rate vs. Light Absorption Rate---##
j_m_max=2.9 #slope coefficient 

##---Maximum Chlorophyll---##
chl_max=1000 #Chlorophyll Content of leaves (umol/m2)

##---Convert Moles to Mass of Chlorophyll---##
chl_mass=0.89351 #1 umol of Chlorophyll= 0.89351 mg Chlorophyll

##---Amount of Ribulose Bisphosphate Produced per Chlorophyll---##
rc=120 #nmol RuBP/ mg Chlorophyll

##---Conversion Coefficient to Convert Chlorophyll Content to Ribulose Bisphosphate Content---##
crc=chl_mass*rc

##---Rubisco Maximum Content---##
rub_max=(chl_max*crc)/1000 #(umol RuBP/m2)

##---Air Temperature---##
t=20 #degrees C

##---Constant Parameter Arrays for Model---##

#I have commented out parameters that I am assuming are variable (for the time being)

s_c=np.linspace(0.019-0.0005,0.019+0.0005,2)#specific leaf area (m2 C/g C)
ra=np.zeros(shape=2)+20.7 #specific rubisco activity (umol CO2/g Rub s)
nm_c=((s_c*(100.0**2))*0.077+20.25)/1000.0#leaf nitrogen (g N/ g C)
flnr=np.zeros(shape=2)+0.65 #fraction of leaf nitrogen in rubisco (g N Rub/g N leaf)
frnr=np.zeros(shape=2)+6.25 #weight fraction of nitrogen in rubisco molecule (g Rub/g N Rub) 
ea_str=pa_con_atmfrac(611*np.exp(17.27*t/(t+273.3))) #saturation vapor pressure of air (Pa-->umol h20.mol air)
rh=np.zeros(shape=2)+0.55 #relative humidity (kPa/kPa)
ea=rh*ea_str #vapor pressure deficit (umol h2O/mol air)
ca=np.zeros(shape=2)+410 #ambient carbon dioxide (umol CO2/mol air)
tau25=np.zeros(shape=2)+2904.12 #specifity coefficient of tau at 25 C (unitless) 
ko25=np.zeros(shape=2)+296100 #Michaelis-Menten kinetic coefficient for oxygen at 25 C(umol/mol) 
kc25=np.zeros(shape=2)+ 296 #Michaelis-Menten kinetic coefficient for carbon dioxide at 25 C (umol/mol)
o=np.zeros(shape=2)+210000 #concentration of ambient oxygen (umol/mol)
g0=np.zeros(shape=2)+0.01 #Ball-Berry stomatal conductance intercept parameter
a=np.zeros(shape=2)+1.6 #Conversion Coefficient between stomatal conductance to water and carbon dioxide 
chl_c=np.zeros(shape=2)+400 #Chlorophyll Content of the Leaf (umol chl/m2)
tl_c=np.zeros(shape=2)+(31+273.15) #Temperature of the Leaf (K)
vwc_c=np.zeros(shape=2)+0.15 #Soil Volumetric Water Content (cm3/cm3)
vwc_min=0.08 #minimum soil water content for photosynthesis to occur (permanent wilting point) (cm3/cm3) 
vwc_max=0.3 #maximum soil water content where increases in soil water do not affect photosynthesis (field capacity?) (cm3/cm3)
q=0.2 #parameter for soil water affect on photosynthesis (unitless)
ij_c=np.zeros(shape=2)+0.96 #leaf area & angle index--downregulates jmax
m=np.zeros(shape=2)+9.0 #ball-berry parameter (unitless)
b=1.37 #Conversion Coefficient between boundary layer conductance to water and carbon dioxide 
dia_c=3.5/100. #Mean diameter or size of leaf (m)
u=5.0 #windspeed (m/s)
#gm=?? (mesophyll conductance)

#---------------Determine if I Want to Keep Any of the Variable Parameters Constant---------------#

const_params=[]
for xx in it.combinations(['dia'],1):
    const_params+=[xx]

#do this when I do not put any of the variable parameters as constant. instead I 
#vary each parameter one at a time while keeping the other parameters constant.
if const_params==[()]:
    const_params=['nan']   


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
    #axA.set_xlim([0,30])
    #axA.set_ylim([0,30])
    #axA2.set_ylim([0,9])
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
    axC.set_ylabel('Leaf N (mol H2O/m2s)',fontsize=36, fontname='Times New Roman')
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
    axE.set_ylabel('Evapotranspiration (umol H2O/m2s)',fontsize=36, fontname='Times New Roman')
    #axA.set_xlim([0,20])
    #axA.set_ylim([0,10])
    axE.set_title('Plant Communities vs. Evapotranspiration', fontname='Times New Roman',fontsize=36,fontweight='bold')

   
    #--figure 6--#
    
    #put in correct ax value (e.g. axA, axB)
    fig6,axF = plt.subplots(figsize=(15,15))

    axF.set_xlabel('NUE (umol CO2/g N s)',fontsize=36, fontname='Times New Roman')
    axF.set_ylabel('WUE (umol CO2/mmol H2O)',fontsize=36, fontname='Times New Roman')
    #axA.set_xlim([0,20])
    #axA.set_ylim([0,10])
    axF.set_title('NUE vs. WUE for all Plant Communities', fontname='Times New Roman',fontsize=36,fontweight='bold')


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
        if 'nm' in const_params[ii]:
            nm=nm_c
        if 's' in const_params[ii]:
            s=s_c 
        if 'dia' in const_params[ii]:
            dia=dia_c
        if 'chl' in const_params[ii]:
            chl=chl_c
        if 'tl' in const_params[ii]:
            tl=tl_c
        if 'vwc' in const_params[ii]:
            vwc=vwc_c
        if 'ij' in const_params[ii]:
            ij=ij_c
       
#---------------Photosynthesis Function---------------#

        #alter this line of code for when implementing different photosynthesis functions
        wue, nue, A, E, cs, ci, gsw, gs, gbw, gb =photo(s,nm,tl,ea,chl,crc,rub_max,ij,vwc,kc25,ko25,o,tau25,ca,rh,m,a,frnr,flnr,ra,j_b,j_m_max,q,vwc_min,vwc_max,b)
        
        #test to make sure wue and nue are positive at not 'nan'
        if wue==-999 and nue==-999:
            continue
    
    
#---------------Make Array of Values for Each Meadow---------------#  
        
        #number of simulations per meadow type:
        m_sim=len(leaf_params)/5.0 #meadow simulations 
        
        if i<(m_sim):  
            wue_f+=[wue]
            nue_f+=[nue]
            na_f+=[na]
            gsw_f+=[gsw]
            A_f+=[A]
            E_f+=[E]
        elif i>=(m_sim) and i<(m_sim*2):
            wue_d+=[wue]
            nue_d+=[nue]
            na_d+=[na]
            gsw_d+=[gsw]
            A_d+=[A]
            E_d+=[E]
        elif i>=(m_sim*2) and i<(m_sim*3):
            wue_m+=[wue]
            nue_m+=[nue]
            na_m+=[na]
            gsw_m+=[gsw]
            A_m+=[A]
            E_m+=[E]
        elif i>=(m_sim*3) and i<(m_sim*4):
            wue_w+=[wue]
            nue_w+=[nue]
            na_w+=[na]
            gsw_w+=[gsw]
            A_w+=[A]
            E_w+=[E]
        elif i>=(m_sim*4) and i<(m_sim*5):
            wue_s+=[wue]
            nue_s+=[nue]
            na_s+=[na]
            gsw_s+=[gsw]
            A_s+=[A]
            E_s+=[E]
        
        nue_tot+=[nue]
        wue_tot+=[wue]

#---------------Plot Plant Communities vs. NUE & WUE---------------#      
    
    nue_bp=axA.boxplot([nue_f,nue_d,nue_m,nue_w,nue_s], positions=[0.875,1.875,2.875,3.875,4.875],widths=0.25, patch_artist=True, showmeans=True, showfliers=False)
    wue_bp=axA2.boxplot([wue_f,wue_d,wue_m,wue_w,wue_s], positions=[1.125,2.125,3.125,4.125,5.125],widths=0.25, patch_artist=True, showmeans=True, showfliers=False)
    axA.plot([0.875,1.875,2.875,3.875,4.875],[np.mean(nue_f),np.mean(nue_d),np.mean(nue_m),np.mean(nue_w),np.mean(nue_s)],'-r')
    axA2.plot([1.125,2.125,3.125,4.125,5.125],[np.mean(wue_f),np.mean(wue_d),np.mean(wue_m),np.mean(wue_w),np.mean(wue_s)],'-b')

    axA.set_xticks([1, 2, 3, 4, 5])
    axA.set_xticklabels(['Fellfield','Dry Meadow','Moist Meadow','Wet Meadow','Snowbed'],fontname='Times New Roman')
    axA.tick_params(axis='x', labelsize=28)
    axA.tick_params(axis='y', labelsize=18)
    axA2.tick_params(axis='y', labelsize=18)
    
    #nue boxplot specs
    for box in nue_bp['boxes']:
        #change outline color
        box.set(color='red',linewidth=2)
        #change fill color
        box.set(facecolor='red',alpha=0.2)

    for whisker in nue_bp['whiskers']:
        whisker.set(color='red',linewidth=2)
    
    for cap in nue_bp['caps']:
        cap.set(color='red',linewidth=2)

    for median in nue_bp['medians']:
        median.set(color='red', linewidth=2)

    for flier in nue_bp['fliers']:
        flier.set(marker='*',color='red',alpha=0.5)

    for means in nue_bp['means']:
        means.set(marker='o',markerfacecolor='red')    

    
    #wue boxplot specs    
    for box in wue_bp['boxes']:
        #change outline color
        box.set(color='blue',linewidth=2)
        #change fill color
        box.set(facecolor='blue',alpha=0.2)

    for whisker in wue_bp['whiskers']:
        whisker.set(color='blue',linewidth=2)
    
    for cap in wue_bp['caps']:
        cap.set(color='blue',linewidth=2)

    for median in wue_bp['medians']:
        median.set(color='blue', linewidth=2)

    for flier in wue_bp['fliers']:
        flier.set(marker='*',color='blue',alpha=0.5)

    for means in wue_bp['means']:
        means.set(marker='o',markerfacecolor='blue')
          

#---------------Box Plot Plant Communities vs. Stomatal Condcutance--------------- #     

    gsw_bp=axB.boxplot([gsw_f,gsw_d,gsw_m,gsw_w,gsw_s], patch_artist=True, showmeans=True, showfliers=False)
    axB.plot([1,2,3,4,5],[np.mean(gsw_f),np.mean(gsw_d),np.mean(gsw_m),np.mean(gsw_w),np.mean(gsw_s)],'-c')
  
    axB.set_xticks([1, 2, 3,4,5])
    axB.set_xticklabels(['Fellfield','Dry Meadow','Moist Meadow','Wet Meadow','Snowbed'],fontname='Times New Roman')
    axB.tick_params(axis='x', labelsize=25)
    axB.tick_params(axis='y', labelsize=16)
    
    #gsw boxplot specs
    for box in gsw_bp['boxes']:
        #change outline color
        box.set(color='cyan',linewidth=2)
        #change fill color
        box.set(facecolor='cyan',alpha=0.5)

    for whisker in gsw_bp['whiskers']:
        whisker.set(color='cyan',linewidth=2)
    
    for cap in gsw_bp['caps']:
        cap.set(color='cyan',linewidth=2)

    for median in gsw_bp['medians']:
        median.set(color='cyan', linewidth=2)

    for flier in gsw_bp['fliers']:
        flier.set(marker='*',color='cyan',alpha=0.5)

    for means in gsw_bp['means']:
        means.set(marker='o',markerfacecolor='cyan')    

#---------------Box Plot Plant Communities vs. Leaf N--------------- #     

    na_bp=axC.boxplot([na_f,na_d,na_m,na_w,na_s], patch_artist=True, showmeans=True, showfliers=False)
    axC.plot([1,2,3,4,5],[np.mean(na_f),np.mean(na_d),np.mean(na_m),np.mean(na_w),np.mean(na_s)],'-',color='orange')
  
    axC.set_xticks([1, 2, 3,4,5])
    axC.set_xticklabels(['Fellfield','Dry Meadow','Moist Meadow','Wet Meadow','Snowbed'],fontname='Times New Roman')
    axC.tick_params(axis='x', labelsize=25)
    axC.tick_params(axis='y', labelsize=16)
    
    #na boxplot specs
    for box in na_bp['boxes']:
        #change outline color
        box.set(color='orange',linewidth=2)
        #change fill color
        box.set(facecolor='orange',alpha=0.5)

    for whisker in na_bp['whiskers']:
        whisker.set(color='orange',linewidth=2)
    
    for cap in na_bp['caps']:
        cap.set(color='orange',linewidth=2)

    for median in na_bp['medians']:
        median.set(color='orange', linewidth=2)

    for flier in na_bp['fliers']:
        flier.set(marker='*',color='orange',alpha=0.5)

    for means in na_bp['means']:
        means.set(marker='o',markerfacecolor='orange')    

        
#---------------Box Plot Plant Communities vs. Assimilation--------------- #     

    A_bp=axD.boxplot([A_f,A_d,A_m,A_w,A_s], patch_artist=True, showmeans=True, showfliers=False)
    axD.plot([1,2,3,4,5],[np.mean(A_f),np.mean(A_d),np.mean(A_m),np.mean(A_w),np.mean(A_s)],'-',color='purple')
    
    axD.set_xticks([1, 2, 3,4,5])
    axD.set_xticklabels(['Fellfield','Dry Meadow','Moist Meadow','Wet Meadow','Snowbed'],fontname='Times New Roman')
    axD.tick_params(axis='x', labelsize=25)
    axD.tick_params(axis='y', labelsize=16)
    
    #A boxplot specs
    for box in A_bp['boxes']:
        #change outline color
        box.set(color='purple',linewidth=2)
        #change fill color
        box.set(facecolor='purple',alpha=0.5)

    for whisker in A_bp['whiskers']:
        whisker.set(color='purple',linewidth=2)
    
    for cap in A_bp['caps']:
        cap.set(color='purple',linewidth=2)

    for median in A_bp['medians']:
        median.set(color='purple', linewidth=2)

    for flier in A_bp['fliers']:
        flier.set(marker='*',color='purple',alpha=0.5)

    for means in A_bp['means']:
        means.set(marker='o',markerfacecolor='purple')  
        
#---------------Box Plot Plant Communities vs. Evapotranspiration--------------- #     

    E_bp=axE.boxplot([E_f,E_d,E_m,E_w,E_s], patch_artist=True, showmeans=True, showfliers=False)
    axE.plot([1,2,3,4,5],[np.mean(E_f),np.mean(E_d),np.mean(E_m),np.mean(E_w),np.mean(E_s)],'-k')
  
    axE.set_xticks([1, 2, 3,4,5])
    axE.set_xticklabels(['Fellfield','Dry Meadow','Moist Meadow','Wet Meadow','Snowbed'],fontname='Times New Roman')
    axE.tick_params(axis='x', labelsize=25)
    axE.tick_params(axis='y', labelsize=16)
    
    
    #E boxplot specs
    for box in E_bp['boxes']:
        #change outline color
        box.set(color='black',linewidth=2)
        #change fill color
        box.set(facecolor='black',alpha=0.5)

    for whisker in E_bp['whiskers']:
        whisker.set(color='black',linewidth=2)
    
    for cap in E_bp['caps']:
        cap.set(color='black',linewidth=2)

    for median in E_bp['medians']:
        median.set(color='black', linewidth=2)

    for flier in E_bp['fliers']:
        flier.set(marker='*',color='black',alpha=0.5)

    for means in E_bp['means']:
        means.set(marker='o',markerfacecolor='black')  
        
        
#---------------Regression Plot WUE vs. NUE---------------#   

    axF.plot([np.mean(nue_f),np.mean(nue_d),np.mean(nue_m),np.mean(nue_w),np.mean(nue_s)],[np.mean(wue_f),np.mean(wue_d),np.mean(wue_m),np.mean(wue_w),np.mean(wue_s)])
    axF.scatter(nue_tot,wue_tot) 
    #axF.tick_params(axis='x', labelsize=25)
    axF.tick_params(axis='y', labelsize=16)
    
#---------------Make Plot Interactive---------------# 
   
#        plt.pause(0.0001)
#        plt.ion()

    
#---------------Finalize Figure---------------#    

    ##---Legend---##
    #axA.legend(bbox_to_anchor=(1, 1), loc='left', prop={'size':15})

    ##---Save Figure--##
    fig1.savefig('NUE_vs_WUE_var_allcoms.png') 
    fig2.savefig('gsw_plantcoms.png') 
    fig3.savefig('na_plantcoms.png') 
    fig4.savefig('A_plantcoms.png')
    fig5.savefig('E_plantcoms.png')
    fig6.savefig('WUE_vs_NUE_regression.png')
