# -*- coding: utf-8 -*-
"""
Photosynthesis and Stomatal Conductance Model 
Created 9/27/2016
Katherine Wentz

This is a program that runs photosynthesis and
stomatal conductance given changes in leaf-
level traits. I derive photosynthesis from a
stomatal conductance model. That way I am 
able to void the ci term. I am breaking up the code
into 2 different models. The first model
pretends that there is no intercept term in
the Ball-Berry stomatal conductance model. The
second model contains the intercept term.
The end product is graphs of NUE vs. WUE.


Update: I am going to run the model for plants with 
traits that are distinctive of the meadow moisture 
gradient in the alpine tundra.

"""

#Chlorophyll, temp of leaf, m, s,nm change the relationship between wue and nue. 
#but nm and sla are dependent variables, so they either have to both be constant or both be varying

#---------------Import Modules---------------#

import itertools as it
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rcParams

#The line of code below is for if I want to input all combinations of changed parameters into my model:
from leaf_parameter_inputs_intraspecific_var import leaf_params_plas

#The line of code below imports the photosynthesis model 
from Photosynthesis_Model import photo

from photo_functions import arr_temp, bol_temp, pa_con_atmfrac

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
#lamb=np.zeros(shape=3)+0.0074 #marginal WUE (umol CO2/umol H2O)
b=np.zeros(shape=2)+0.0 #Ball-Berry stomatal conductance intercept parameter
a=np.zeros(shape=2)+1.6 #Conversion Coefficient between stomatal conductance to water and carbon dioxide 
chl_c=np.zeros(shape=2)+400 #Chlorophyll Content of the Leaf (umol chl/m2)
tl_c=np.zeros(shape=2)+(31+273.15) #Temperature of the Leaf (K)
vwc_c=np.zeros(shape=2)+0.15 #Soil Volumetric Water Content (cm3/cm3)
vwc_min=0.08 #minimum soil water content for photosynthesis to occur (permanent wilting point) (cm3/cm3) 
vwc_max=0.3 #maximum soil water content where increases in soil water do not affect photosynthesis (field capacity?) (cm3/cm3)
q=0.2 #parameter for soil water affect on photosynthesis (unitless)
ij_c=np.zeros(shape=2)+0.96 #leaf area & angle index--downregulates jmax
m_c=np.zeros(shape=2)+15.0 #ball-berry parameter (unitless)

##--What Variable Parameters are Constant?
dict_params=[]
for xx in it.combinations(['nm','chl','tl','m','ij','vwc'],0): 
    dict_params+=[xx]

if dict_params==[()]:
    dict_params=['nan']   

for ii in range(len(dict_params)):
    
    #---------------Initialize Plot---------------#

    ##---Figure With Subplots Blueprint---##

    fb1=plt.figure(figsize=(30,40)) 
#    axA = fb1.add_subplot(321)
#    axB = fb1.add_subplot(322)
#    axC = fb1.add_subplot(323)
#    axD = fb1.add_subplot(324)
#    axE = fb1.add_subplot(325)
    axA = plt.subplot2grid((2,6),(0,0),colspan=2)
    axA2=axA.twinx()
    axB = plt.subplot2grid((2,6),(0,2),colspan=2)
    axB2=axB.twinx()
    axC = plt.subplot2grid((2,6),(0,4),colspan=2)
    axC2=axC.twinx()
    axD = plt.subplot2grid((2,6),(1,1),colspan=2)
    axD2=axD.twinx()
    axE = plt.subplot2grid((2,6),(1,3),colspan=2)
    axE2=axE.twinx()
    
    ##---Figure Without Subplots Blueprint---##

    #put in correct ax value (e.g. axA, axB)
#    fig,axA = plt.subplots(figsize=(10,5))


    ##---Define Plot Parameters Based on Graph Interests---##

    axA.set_xlabel('Percent Change in Trait Plasticity (%)',fontsize=12, fontname='Times New Roman')
    axA.set_ylabel('NUE (umol CO2/g N s)',fontsize=12, fontname='Times New Roman')
    axA2.set_ylabel('WUE (umol CO2/mmol H2O)',fontsize=12, fontname='Times New Roman')
    #    axA.set_xlim([0,25])
#    axA.set_ylim([0,10])
    axA.set_title(' Fellfield: Variation in Leaf Index of Light Range', fontname='Times New Roman',fontsize=14,fontweight='bold')
#    axA.set_title('Growth Response: constant %s, %s, %s, %s, %s' % (dict_params[ii][0],dict_params[ii][1],dict_params[ii][2],dict_params[ii][3],dict_params[ii][4]), fontname='Times New Roman',fontsize=23,fontweight='bold')

    axB.set_xlabel('Percent Change in Trait Plasticity (%)',fontsize=12, fontname='Times New Roman')
    axB.set_ylabel('NUE (umol CO2/g N s)',fontsize=12, fontname='Times New Roman')
    axB2.set_ylabel('WUE (umol CO2/mmol H2O)',fontsize=12, fontname='Times New Roman')
#    axA.set_xlim([0,25])
#    axA.set_ylim([0,10])
    axB.set_title('Dry Meadow: Variation in Leaf Index of Light Range', fontname='Times New Roman',fontsize=14,fontweight='bold')
#    axA.set_title('Growth Response: constant %s, %s, %s, %s, %s' % (dict_params[ii][0],dict_params[ii][1],dict_params[ii][2],dict_params[ii][3],dict_params[ii][4]), fontname='Times New Roman',fontsize=23,fontweight='bold')

    axC.set_xlabel('Percent Change in Trait Plasticity (%)',fontsize=12, fontname='Times New Roman')
    axC.set_ylabel('NUE (umol CO2/g N s)',fontsize=12, fontname='Times New Roman')
    axC2.set_ylabel('WUE (umol CO2/mmol H2O)',fontsize=12, fontname='Times New Roman')
#    axA.set_xlim([0,25])
#    axA.set_ylim([0,10])
    axC.set_title('Moist Meadow: Variation in Leaf Index of Light Range', fontname='Times New Roman',fontsize=14,fontweight='bold')
#    axA.set_title('Growth Response: constant %s, %s, %s, %s, %s' % (dict_params[ii][0],dict_params[ii][1],dict_params[ii][2],dict_params[ii][3],dict_params[ii][4]), fontname='Times New Roman',fontsize=23,fontweight='bold')

    axD.set_xlabel('Percent Change in Trait Plasticity (%)',fontsize=12, fontname='Times New Roman')
    axD.set_ylabel('NUE (umol CO2/g N s)',fontsize=12, fontname='Times New Roman')
    axD2.set_ylabel('WUE (umol CO2/mmol H2O)',fontsize=12, fontname='Times New Roman')
#    axA.set_xlim([0,25])
#    axA.set_ylim([0,10])
    axD.set_title('Wet Meadow: Variation in Leaf Index of Light Range', fontname='Times New Roman',fontsize=14,fontweight='bold')
#    axA.set_title('Growth Response: constant %s, %s, %s, %s, %s' % (dict_params[ii][0],dict_params[ii][1],dict_params[ii][2],dict_params[ii][3],dict_params[ii][4]), fontname='Times New Roman',fontsize=23,fontweight='bold')

    axE.set_xlabel('Percent Change in Trait Plasticity (%)',fontsize=12, fontname='Times New Roman')
    axE.set_ylabel('NUE (umol CO2/g N s)',fontsize=12, fontname='Times New Roman')
    axE2.set_ylabel('WUE (umol CO2/mmol H2O)',fontsize=12, fontname='Times New Roman')
#    axA.set_xlim([0,25])
#    axA.set_ylim([0,10])
    axE.set_title('Snowbed: Variation in Leaf Index of Light Range', fontname='Times New Roman',fontsize=14,fontweight='bold')
#    axA.set_title('Growth Response: constant %s, %s, %s, %s, %s' % (dict_params[ii][0],dict_params[ii][1],dict_params[ii][2],dict_params[ii][3],dict_params[ii][4]), fontname='Times New Roman',fontsize=23,fontweight='bold')


    ##---Line Type for Each Plant---##
#    n=32 #number of variable parameter combinations for each meadow type
#
#    color=['darkblue','slateblue','blueviolet','mediumorchid','purple']
#    legend=['90% Trait Variation','70% Trait Variation','50% Trait Variation','30% Trait Variation','10% Trait Variation']
#    
#    linewidth=[8.0,6.0,4.0,2.0,0.5]

    wue_f_per1=[]
    nue_f_per1=[]
    wue_f_per2=[]
    nue_f_per2=[]
    wue_f_per3=[]
    nue_f_per3=[]

    wue_d_per1=[]
    nue_d_per1=[]
    wue_d_per2=[]
    nue_d_per2=[]
    wue_d_per3=[]
    nue_d_per3=[]

    wue_m_per1=[]
    nue_m_per1=[]
    wue_m_per2=[]
    nue_m_per2=[]
    wue_m_per3=[]
    nue_m_per3=[]

    wue_w_per1=[]
    nue_w_per1=[]
    wue_w_per2=[]
    nue_w_per2=[]
    wue_w_per3=[]
    nue_w_per3=[]

    wue_s_per1=[]
    nue_s_per1=[]
    wue_s_per2=[]
    nue_s_per2=[]
    wue_s_per3=[]
    nue_s_per3=[]
  
    ##---Variable Parameter Arrays for Model---##
    
    for iii in range(len(leaf_params_plas)):
        ##---Initialize Arrays for Each Meadow---##
#        wue_f=[]
#        nue_f=[]
#        wue_d=[]
#        nue_d=[]
#        wue_m=[]
#        nue_m=[]
#        wue_w=[]
#        nue_w=[]    
#        wue_s=[]
#        nue_s=[]

        for i in range(len(leaf_params_plas[iii])):
            for key,val in leaf_params_plas[iii][i].items():
                exec(key + '=val')
            #set variables constant
            if 'nm' in dict_params[ii]:
                nm=nm_c
                s=s_c #nm and s are dependent variables
            if 'm' in dict_params[ii]:
                m=m_c
            if 'chl' in dict_params[ii]:
                chl=chl_c
            if 'tl' in dict_params[ii]:
                tl=tl_c
            if 'vwc' in dict_params[ii]:
                vwc=vwc_c
            if 'ij' in dict_params[ii]:
                ij=ij_c
            
            
            wue, nue, A, gsw, E, na=photo(s,nm,tl,ea,chl,crc,rub_max,ij,vwc,kc25,ko25,o,tau25,ca,rh,m,a,frnr,flnr,ra,j_b,j_m_max,q,vwc_min,vwc_max,b)
            
            if wue==-999 and nue==-999:
                continue

    
#---------------Test for Low NUE Values---------------#  

#    if any(nue<15):
#        break
    

    
#---------------Make Array of Values for Each Meadow---------------#  

            if i+1<33:
                if iii==0:
                    wue_f_per1+=[wue]
                    nue_f_per1+=[nue]
                if iii==1:
                    wue_f_per2+=[wue]
                    nue_f_per2+=[nue]
                if iii==2:
                    wue_f_per3+=[wue]
                    nue_f_per3+=[nue]
            elif i+1>=33 and i+1<65:
                if iii==0:
                    wue_d_per1+=[wue]
                    nue_d_per1+=[nue]
                if iii==1:
                    wue_d_per2+=[wue]
                    nue_d_per2+=[nue]
                if iii==2:
                    wue_d_per3+=[wue]
                    nue_d_per3+=[nue]
            elif i+1>=65 and i+1<97:
                if iii==0:
                    wue_m_per1+=[wue]
                    nue_m_per1+=[nue]
                if iii==1:
                    wue_m_per2+=[wue]
                    nue_m_per2+=[nue]
                if iii==2:
                    wue_m_per3+=[wue]
                    nue_m_per3+=[nue]             
            elif i+1>=97 and i+1<129:
                if iii==0:
                    wue_w_per1+=[wue]
                    nue_w_per1+=[nue]
                if iii==1:
                    wue_w_per2+=[wue]
                    nue_w_per2+=[nue]
                if iii==2:
                    wue_w_per3+=[wue]
                    nue_w_per3+=[nue] 
            elif i+1>=129 and i+1<161:
                if iii==0:
                    wue_s_per1+=[wue]
                    nue_s_per1+=[nue]
                if iii==1:
                    wue_s_per2+=[wue]
                    nue_s_per2+=[nue]
                if iii==2:
                    wue_s_per3+=[wue]
                    nue_s_per3+=[nue]
            
#            wue_f_tot+=wue_f
#            nue_f_tot+=nue_f
#            wue_d_tot+=wue_d
#            nue_d_tot+=nue_d
#            wue_m_tot+=wue_m
#            nue_m_tot+=nue_m
#            wue_w_tot+=wue_w
#            nue_w_tot+=nue_w
#            wue_s_tot+=wue_s
#            nue_s_tot+=nue_s
   
#---------------Plot NUE vs. WUE---------------#      

    #I am plotting NUE vs. WUE for each plant trait
    
    ##---Separate Plots Into Different Figures: Set Up---##
    #fb=plt.figure(i+1,figsize=(6,6)) #figure blueprint
    #fig,ax1 = plt.subplots()
    #ax1.set_xlabel('NUE (umol CO2/g N s)',fontsize=12)
    #ax1.set_ylabel('WUE (umol CO2/umol H20)',fontsize=12)
    #ax1.set_title('NUE vs. WUE',fontsize=14)
    
    ##---Plot---##
    #will need to change ax value if plotting separate graphs for each iteration, e.g. ax1
#    axA.plot(nue,wue,label='%s' %trait[i], color='%s' %color[i],marker='%s' %marker[i],linestyle='%s' %style[i]) 
   

    nue_bp_f=axA.boxplot([nue_f_per1,nue_f_per2,nue_f_per3], positions=[0.875,1.875,2.875],widths=0.25, patch_artist=True, showmeans=True,showfliers=False)
    wue_bp_f=axA2.boxplot([wue_f_per1,wue_f_per2,wue_f_per3], positions=[1.125,2.125,3.125],widths=0.25, patch_artist=True, showmeans=True,showfliers=False)
    axA.plot([0.875,1.875,2.875],[np.mean(nue_f_per1),np.mean(nue_f_per2),np.mean(nue_f_per3)],'-r')
    axA2.plot([1.125,2.125,3.125],[np.mean(wue_f_per1),np.mean(wue_f_per2),np.mean(wue_f_per3)],'-b')
    
    axA.plot([0.875,1.875,2.875],[np.max(nue_f_per1),np.max(nue_f_per2),np.max(nue_f_per3)],'--r')
    axA.plot([0.875,1.875,2.875],[np.min(nue_f_per1),np.min(nue_f_per2),np.min(nue_f_per3)],'--r')
    axA2.plot([1.125,2.125,3.125],[np.max(wue_f_per1),np.max(wue_f_per2),np.max(wue_f_per3)],'--b')
    axA2.plot([1.125,2.125,3.125],[np.min(wue_f_per1),np.min(wue_f_per2),np.min(wue_f_per3)],'--b')

    axA.set_xticks([1, 2, 3])
    axA.set_xticklabels(['-40%','0','+40%'])
    rcParams['xtick.labelsize']=10   
    
    nue_bp_d=axB.boxplot([nue_d_per1,nue_d_per2,nue_d_per3], positions=[0.875,1.875,2.875],widths=0.25, patch_artist=True, showmeans=True,showfliers=False)
    wue_bp_d=axB2.boxplot([wue_d_per1,wue_d_per2,wue_d_per3], positions=[1.125,2.125,3.125],widths=0.25, patch_artist=True, showmeans=True,showfliers=False)
    axB.plot([0.875,1.875,2.875],[np.mean(nue_d_per1),np.mean(nue_d_per2),np.mean(nue_d_per3)],'-r')
    axB2.plot([1.125,2.125,3.125],[np.mean(wue_d_per1),np.mean(wue_d_per2),np.mean(wue_d_per3)],'-b')

    axB.plot([0.875,1.875,2.875],[np.max(nue_d_per1),np.max(nue_d_per2),np.max(nue_d_per3)],'--r')
    axB.plot([0.875,1.875,2.875],[np.min(nue_d_per1),np.min(nue_d_per2),np.min(nue_d_per3)],'--r')
    axB2.plot([1.125,2.125,3.125],[np.max(wue_d_per1),np.max(wue_d_per2),np.max(wue_d_per3)],'--b')
    axB2.plot([1.125,2.125,3.125],[np.min(wue_d_per1),np.min(wue_d_per2),np.min(wue_d_per3)],'--b')
    
    axB.set_xticks([1, 2, 3])
    axB.set_xticklabels(['-40%','0','+40%'])    
    rcParams['xtick.labelsize']=10   
    
    nue_bp_m=axC.boxplot([nue_m_per1,nue_m_per2,nue_m_per3], positions=[0.875,1.875,2.875],widths=0.25, patch_artist=True, showmeans=True,showfliers=False)
    wue_bp_m=axC2.boxplot([wue_m_per1,wue_m_per2,wue_m_per3], positions=[1.125,2.125,3.125],widths=0.25, patch_artist=True, showmeans=True,showfliers=False)
    axC.plot([0.875,1.875,2.875],[np.mean(nue_m_per1),np.mean(nue_m_per2),np.mean(nue_m_per3)],'-r')
    axC2.plot([1.125,2.125,3.125],[np.mean(wue_m_per1),np.mean(wue_m_per2),np.mean(wue_m_per3)],'-b')

    axC.plot([0.875,1.875,2.875],[np.max(nue_m_per1),np.max(nue_m_per2),np.max(nue_m_per3)],'--r')
    axC.plot([0.875,1.875,2.875],[np.min(nue_m_per1),np.min(nue_m_per2),np.min(nue_m_per3)],'--r')
    axC2.plot([1.125,2.125,3.125],[np.max(wue_m_per1),np.max(wue_m_per2),np.max(wue_m_per3)],'--b')
    axC2.plot([1.125,2.125,3.125],[np.min(wue_m_per1),np.min(wue_m_per2),np.min(wue_m_per3)],'--b')

    axC.set_xticks([1, 2, 3])
    axC.set_xticklabels(['-40%','0','+40%'])     
    rcParams['xtick.labelsize']=10   
    
    nue_bp_w=axD.boxplot([nue_w_per1,nue_w_per2,nue_w_per3], positions=[0.875,1.875,2.875],widths=0.25, patch_artist=True, showmeans=True,showfliers=False)
    wue_bp_w=axD2.boxplot([wue_w_per1,wue_w_per2,wue_w_per3], positions=[1.125,2.125,3.125],widths=0.25, patch_artist=True, showmeans=True,showfliers=False)
    axD.plot([0.875,1.875,2.875],[np.mean(nue_w_per1),np.mean(nue_w_per2),np.mean(nue_w_per3)],'-r')
    axD2.plot([1.125,2.125,3.125],[np.mean(wue_w_per1),np.mean(wue_w_per2),np.mean(wue_w_per3)],'-b')

    axD.plot([0.875,1.875,2.875],[np.max(nue_w_per1),np.max(nue_w_per2),np.max(nue_w_per3)],'--r')
    axD.plot([0.875,1.875,2.875],[np.min(nue_w_per1),np.min(nue_w_per2),np.min(nue_w_per3)],'--r')
    axD2.plot([1.125,2.125,3.125],[np.max(wue_w_per1),np.max(wue_w_per2),np.max(wue_w_per3)],'--b')
    axD2.plot([1.125,2.125,3.125],[np.min(wue_w_per1),np.min(wue_w_per2),np.min(wue_w_per3)],'--b')
 
    axD.set_xticks([1, 2, 3])
    axD.set_xticklabels(['-40%','0','+40%'])        
    rcParams['xtick.labelsize']=10   
    
    nue_bp_s=axE.boxplot([nue_s_per1,nue_s_per2,nue_s_per3], positions=[0.875,1.875,2.875],widths=0.25, patch_artist=True, showmeans=True,showfliers=False)
    wue_bp_s=axE2.boxplot([wue_s_per1,wue_s_per2,wue_s_per3], positions=[1.125,2.125,3.125],widths=0.25, patch_artist=True, showmeans=True,showfliers=False)
    axE.plot([0.875,1.875,2.875],[np.mean(nue_s_per1),np.mean(nue_s_per2),np.mean(nue_s_per3)],'-r')
    axE2.plot([1.125,2.125,3.125],[np.mean(wue_s_per1),np.mean(wue_s_per2),np.mean(wue_s_per3)],'-b')

    axE.plot([0.875,1.875,2.875],[np.max(nue_s_per1),np.max(nue_s_per2),np.max(nue_s_per3)],'--r')
    axE.plot([0.875,1.875,2.875],[np.min(nue_s_per1),np.min(nue_s_per2),np.min(nue_s_per3)],'--r')
    axE2.plot([1.125,2.125,3.125],[np.max(wue_s_per1),np.max(wue_s_per2),np.max(wue_s_per3)],'--b')
    axE2.plot([1.125,2.125,3.125],[np.min(wue_s_per1),np.min(wue_s_per2),np.min(wue_s_per3)],'--b')

    axE.set_xticks([1, 2, 3])
    axE.set_xticklabels(['-40%','0','+40%'])        
    rcParams['xtick.labelsize']=10   
    
    nue_bp=[nue_bp_f,nue_bp_d,nue_bp_m,nue_bp_w,nue_bp_s]     
    wue_bp=[wue_bp_f,wue_bp_d,wue_bp_m,wue_bp_w,wue_bp_s]           

    #nue boxplot specs
    for jj in range(len(nue_bp)):
        for box in nue_bp[jj]['boxes']:
            #change outline color
            box.set(color='red',linewidth=2)
            #change fill color
            box.set(facecolor='red',alpha=0.2)
    
        for whisker in nue_bp[jj]['whiskers']:
            whisker.set(color='red',linewidth=2,linestyle='-')
        
        for cap in nue_bp[jj]['caps']:
            cap.set(color='red',linewidth=2)
    
        for median in nue_bp[jj]['medians']:
            median.set(color='red', linewidth=2)
    
        for flier in nue_bp[jj]['fliers']:
            flier.set(marker='*',color='red',alpha=0.5)
    
        for means in nue_bp[jj]['means']:
            means.set(marker='o',markerfacecolor='red')    


    #wue boxplot specs  
    for jj in range(len(wue_bp)):
        for box in wue_bp[jj]['boxes']:
            #change outline color
            box.set(color='blue',linewidth=2)
            #change fill color
            box.set(facecolor='blue',alpha=0.2)
    
        for whisker in wue_bp[jj]['whiskers']:
            whisker.set(color='blue',linewidth=2,linestyle='-')
        
        for cap in wue_bp[jj]['caps']:
            cap.set(color='blue',linewidth=2)
    
        for median in wue_bp[jj]['medians']:
            median.set(color='blue', linewidth=2)
    
        for flier in wue_bp[jj]['fliers']:
            flier.set(marker='*',color='blue',alpha=0.5)
    
        for means in wue_bp[jj]['means']:
            means.set(marker='o',markerfacecolor='blue')
     
        
        
#        axA.plot([max(nue_f),max(nue_f),min(nue_f),min(nue_f),max(nue_f)],[max(wue_f),min(wue_f),min(wue_f),max(wue_f),max(wue_f)],color=color[iii], linestyle='-', linewidth=linewidth[iii], label=legend[iii]) 
##    axA.scatter([max(nue_f),max(nue_f),min(nue_f),min(nue_f),max(nue_f)],[max(wue_f),min(wue_f),min(wue_f),max(wue_f),max(wue_f)],color='k', marker='^', label='fellfield') 
#        axA.fill_between([np.min(nue_f),np.max(nue_f)],np.min(wue_f),np.max(wue_f),color=color[iii],alpha=0.2) 
#        axA.set_ylim([np.min(wue_f_tot)-np.min(wue_f_tot)*0.01,np.max(wue_f_tot)+np.max(wue_f_tot)*0.01])
#        axA.set_xlim([np.min(nue_f_tot)-0.5,np.max(nue_f_tot)+0.5])
# 
#    
#        axB.plot([max(nue_d),max(nue_d),min(nue_d),min(nue_d),max(nue_d)],[max(wue_d),min(wue_d),min(wue_d),max(wue_d),max(wue_d)],color=color[iii], linestyle='-',linewidth=linewidth[iii], label=legend[iii])
###    axA.scatter([max(nue_d),max(nue_d),min(nue_d),min(nue_d),max(nue_d)],[max(wue_d),min(wue_d),min(wue_d),max(wue_d),max(wue_d)],color='r', marker='d',label='dry meadow')     
#        axB.fill_between([np.min(nue_d),np.max(nue_d)],np.min(wue_d),np.max(wue_d),color=color[iii],alpha=0.2) 
#        axB.set_ylim([np.min(wue_d_tot)-np.min(wue_d_tot)*0.01,np.max(wue_d_tot)+np.max(wue_d_tot)*0.01])
#        axB.set_xlim([np.min(nue_d_tot)-0.5,np.max(nue_d_tot)+0.5])
# 
#        
#        axC.plot([max(nue_m),max(nue_m),min(nue_m),min(nue_m),max(nue_m)],[max(wue_m),min(wue_m),min(wue_m),max(wue_m),max(wue_m)],color=color[iii], linestyle='-',linewidth=linewidth[iii], label=legend[iii])
###    axA.scatter([max(nue_m),max(nue_m),min(nue_m),min(nue_m),max(nue_m)],[max(wue_m),min(wue_m),min(wue_m),max(wue_m),max(wue_m)],color='y', marker='o',s=40,label='moist meadow') 
#        axC.fill_between([np.min(nue_m),np.max(nue_m)],np.min(wue_m),np.max(wue_m),color=color[iii],alpha=0.2) 
#        axC.set_ylim([np.min(wue_m_tot)-np.min(wue_m_tot)*0.01,np.max(wue_m_tot)+np.max(wue_m_tot)*0.01])
#        axC.set_xlim([np.min(nue_m_tot)-0.5,np.max(nue_m_tot)+0.5])
#        
#        axD.plot([max(nue_w),max(nue_w),min(nue_w),min(nue_w),max(nue_w)],[max(wue_w),min(wue_w),min(wue_w),max(wue_w),max(wue_w)],color=color[iii], linestyle='-',linewidth=linewidth[iii],label=legend[iii])
###    axA.scatter([max(nue_w),max(nue_w),min(nue_w),min(nue_w),max(nue_w)],[max(wue_w),min(wue_w),min(wue_w),max(wue_w),max(wue_w)],facecolor='none',edgecolor='g', marker='o',s=120,label='wet meadow')    
#        axD.fill_between([np.min(nue_w),np.max(nue_w)],np.min(wue_w),np.max(wue_w),color=color[iii],alpha=0.2) 
#        axD.set_ylim([np.min(wue_w_tot)-np.min(wue_w_tot)*0.01,np.max(wue_w_tot)+np.max(wue_w_tot)*0.01])
#        axD.set_xlim([np.min(nue_w_tot)-0.5,np.max(nue_w_tot)+0.5])
#        
#        axE.plot([max(nue_s),max(nue_s),min(nue_s),min(nue_s),max(nue_s)],[max(wue_s),min(wue_s),min(wue_s),max(wue_s),max(wue_s)],color=color[iii], linestyle='-',linewidth=linewidth[iii], label=legend[iii])
###    axA.scatter([max(nue_s),max(nue_s),min(nue_s),min(nue_s),max(nue_s)],[max(wue_s),min(wue_s),min(wue_s),max(wue_s),max(wue_s)],color='b', marker='*', label='snowbed')         
#        axE.fill_between([np.min(nue_s),np.max(nue_s)],np.min(wue_s),np.max(wue_s),color=color[iii],alpha=0.2) 
#        axE.set_ylim([np.min(wue_s_tot)-np.min(wue_s_tot)*0.01,np.max(wue_s_tot)+np.max(wue_s_tot)*0.01])
#        axE.set_xlim([np.min(nue_s_tot)-0.5,np.max(nue_s_tot)+0.5])
#   
    
    ##---Separate Plots Into Different Figures: Legend---##
    #ax1.legend(loc=4)

##---------------Nutrient Use Efficiency and Water Use Efficiency--------------- #     

##I am plotting the change in plant trait vs. two y axes: WUE and NUE

        #fb3=plt.figure(3,figsize=(6,6)) #figure blueprint
        #fig, ax1 = plt.subplots()
#        axB=axA.twinx() #second y axis
#        axA.set_xlabel('m (umol CO2/m2s)',fontsize=12)
#        axA.set_ylabel('NUE (umol CO2/g N s)',fontsize=12)
#        axB.set_ylabel('WUE (umol CO2/mol H20)',fontsize=12)
#        axA.scatter(m[0:2],nue,color=color[i],label='NUE') 
#        axB.scatter(m[0:2],wue,color=color[i],label='WUE') 
#        ax1.set_title('NUE vs. WUE with Increasing Photosynthesis',fontsize=14)
##        ax1.legend(loc=2)
#        ax2.legend(loc=1)
        

#---------------Make Plot Interactive---------------# 
#    
#        plt.pause(0.0001)
#        plt.ion()
    #end of sensitivity analysis iterations
    
#---------------Finalize Figure---------------#    

    #axA refers to first figure in subplot; axB refers to second figure in subplot
    #if only one axis is run then the figure is just one plot

    ##---Legend---##
#    axA.legend(bbox_to_anchor=(2, -2), loc='upper', prop={'size':11})
   
#    axA.set_ylim([0,10])
    fb1.tight_layout() #makes it so subplots text does not overlap

    ##---Save Figure--##
    plt.savefig('NUE_vs_WUE_variable_ij_range.png') 
