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
from leaf_parameter_inputs import leaf_params

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

nue_tot=[]
wue_tot=[]
for ii in range(len(dict_params)):
    
    #---------------Initialize Plot---------------#

    ##---Figure With Subplots Blueprint---##

    #fb1=plt.figure(1,figsize=(12,2)) 
    #axA = fb1.add_subplot(121)
    #axB = fb1.add_subplot(122)

    ##---Figure Without Subplots Blueprint---##
    
    #--figure 1--
    
    #put in correct ax value (e.g. axA, axB)
    fig1,axA = plt.subplots(figsize=(30,15))
    
    #twin axis
    axA2=axA.twinx()
    
    ##---Define Plot Parameters Based on Graph Interests---##

#    axA.set_xlabel('Plant Communities',fontsize=20, fontname='Times New Roman')
    axA.set_ylabel('NUE (umol CO2/g N s)',fontsize=36, fontname='Times New Roman')
    axA2.set_ylabel('WUE (umol CO2/mmol H2O)', fontsize=36, fontname='Times New Roman')  
#    axA.set_xlim([0,20])
#    axA.set_ylim([0,10])
    axA.set_title('WUE and NUE for Plant Communities', fontname='Times New Roman',fontsize=36,fontweight='bold')

    
    #--figure 2--
    
    #put in correct ax value (e.g. axA, axB)
    fig2,axB = plt.subplots(figsize=(15,15))
    
    ##---Define Plot Parameters Based on Graph Interests---##

#    axA.set_xlabel('Plant Communities',fontsize=20, fontname='Times New Roman')
    axB.set_ylabel('gsw (mol H2O/m2s)',fontsize=36, fontname='Times New Roman')
#    axA.set_xlim([0,20])
#    axA.set_ylim([0,10])
    axB.set_title('Plant Communities vs. Stomatal Conductance', fontname='Times New Roman',fontsize=36,fontweight='bold')

    #--figure 3--
    
    #put in correct ax value (e.g. axA, axB)
    fig3,axC = plt.subplots(figsize=(15,15))
    
    ##---Define Plot Parameters Based on Graph Interests---##

#    axA.set_xlabel('Plant Communities',fontsize=20, fontname='Times New Roman')
    axC.set_ylabel('Leaf N (mol H2O/m2s)',fontsize=36, fontname='Times New Roman')
#    axA.set_xlim([0,20])
#    axA.set_ylim([0,10])
    axC.set_title('Plant Communities vs. Leaf Nitrogen', fontname='Times New Roman',fontsize=36,fontweight='bold')

    #--figure 4--
    
    #put in correct ax value (e.g. axA, axB)
    fig4,axD = plt.subplots(figsize=(15,15))
    
    ##---Define Plot Parameters Based on Graph Interests---##

#    axA.set_xlabel('Plant Communities',fontsize=20, fontname='Times New Roman')
    axD.set_ylabel('Assimilation (umol CO2/m2s)',fontsize=36, fontname='Times New Roman')
#    axA.set_xlim([0,20])
#    axA.set_ylim([0,10])
    axD.set_title('Plant Communities vs. Assimilation', fontname='Times New Roman',fontsize=36,fontweight='bold')

  
      #--figure 5--
    
    #put in correct ax value (e.g. axA, axB)
    fig5,axE = plt.subplots(figsize=(15,15))
    
    ##---Define Plot Parameters Based on Graph Interests---##

#    axA.set_xlabel('Plant Communities',fontsize=20, fontname='Times New Roman')
    axE.set_ylabel('Evapotranspiration (umol H2O/m2s)',fontsize=36, fontname='Times New Roman')
#    axA.set_xlim([0,20])
#    axA.set_ylim([0,10])
    axE.set_title('Plant Communities vs. Evapotranspiration', fontname='Times New Roman',fontsize=36,fontweight='bold')

   
    #--figure 6--
    
    #put in correct ax value (e.g. axA, axB)
    fig6,axF = plt.subplots(figsize=(15,15))
    
    ##---Define Plot Parameters Based on Graph Interests---##

    axF.set_xlabel('NUE (umol CO2/g N s)',fontsize=36, fontname='Times New Roman')
    axF.set_ylabel('WUE (umol CO2/mmol H2O)',fontsize=36, fontname='Times New Roman')
#    axA.set_xlim([0,20])
#    axA.set_ylim([0,10])
    axF.set_title('NUE vs. WUE for all Plant Communities', fontname='Times New Roman',fontsize=36,fontweight='bold')

    
    
  
    ##---Line Type for Each Plant---##
#    n=16 #number of variable parameter combinations for each meadow type
#
#    color=['k']*n+['r']*n+['y']*n+['g']*n+['b']*n
#    
#    marker=['d']*n+['o']*n+['v']*n+['*']*n+['^']*n

    ##---Initialize Arrays for Each Meadow---##
    
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

    ##---Variable Parameter Arrays for Model---##
    for i in range(len(leaf_params)):
        for key,val in leaf_params[i].items():
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
            wue_f+=[wue]
            nue_f+=[nue]
            na_f+=[na]
            gsw_f+=[gsw]
            A_f+=[A]
            E_f+=[E]
        elif i+1>=33 and i+1<65:
            wue_d+=[wue]
            nue_d+=[nue]
            na_d+=[na]
            gsw_d+=[gsw]
            A_d+=[A]
            E_d+=[E]
        elif i+1>=65 and i+1<97:
            wue_m+=[wue]
            nue_m+=[nue]
            na_m+=[na]
            gsw_m+=[gsw]
            A_m+=[A]
            E_m+=[E]
        elif i+1>=97 and i+1<129:
            wue_w+=[wue]
            nue_w+=[nue]
            na_w+=[na]
            gsw_w+=[gsw]
            A_w+=[A]
            E_w+=[E]
        elif i+1>=129 and i+1<161:
            wue_s+=[wue]
            nue_s+=[nue]
            na_s+=[na]
            gsw_s+=[gsw]
            A_s+=[A]
            E_s+=[E]
        
        nue_tot+=[nue]
        wue_tot+=[wue]

#---------------Plot Plant Communities NUE vs. WUE---------------#      

    #I am plotting NUE vs. WUE for each plant community
    
    ##---Plot---##
    
    nue_bp=axA.boxplot([nue_f,nue_d,nue_m,nue_w,nue_s], positions=[0.875,1.875,2.875,3.875,4.875],widths=0.25, patch_artist=True, showmeans=True, showfliers=False)
    wue_bp=axA2.boxplot([wue_f,wue_d,wue_m,wue_w,wue_s], positions=[1.125,2.125,3.125,4.125,5.125],widths=0.25, patch_artist=True, showmeans=True, showfliers=False)
    axA.plot([0.875,1.875,2.875,3.875,4.875],[np.mean(nue_f),np.mean(nue_d),np.mean(nue_m),np.mean(nue_w),np.mean(nue_s)],'-r')
    axA2.plot([1.125,2.125,3.125,4.125,5.125],[np.mean(wue_f),np.mean(wue_d),np.mean(wue_m),np.mean(wue_w),np.mean(wue_s)],'-b')

    axA.set_xticks([1, 2, 3,4,5])
    axA.set_xticklabels(['Fellfield','Dry Meadow','Moist Meadow','Wet Meadow','Snowbed'],fontname='Times New Roman')
    axA.tick_params(axis='x', labelsize=28)
    
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
        
    #will need to change ax value if plotting separate graphs for each iteration, e.g. ax1
#    axA.plot(nue,wue,label='%s' %trait[i], color='%s' %color[i],marker='%s' %marker[i],linestyle='%s' %style[i]) 
   
#    axA.plot([max(nue_f),max(nue_f),min(nue_f),min(nue_f),max(nue_f)],[max(wue_f),min(wue_f),min(wue_f),max(wue_f),max(wue_f)],color='red', linestyle='-', linewidth=8.5, label='fellfield') 
##    axA.scatter([max(nue_f),max(nue_f),min(nue_f),min(nue_f),max(nue_f)],[max(wue_f),min(wue_f),min(wue_f),max(wue_f),max(wue_f)],color='k', marker='^', label='fellfield') 
#    axA.fill_between([np.min(nue_f),np.max(nue_f)],np.min(wue_f),np.max(wue_f),color='red',alpha=0.5) 
##    axA.set_ylim([np.min(wue_f)-np.min(wue_f)*0.01,np.max(wue_f)+np.max(wue_f)*0.01])
##    axA.set_xlim([np.min(nue_f)-0.5,np.max(nue_f)+0.5])  
#    
#    axA.plot([max(nue_d),max(nue_d),min(nue_d),min(nue_d),max(nue_d)],[max(wue_d),min(wue_d),min(wue_d),max(wue_d),max(wue_d)],color='brown', linestyle='-',linewidth=6.5, label='dry meadow')
##    axA.scatter([max(nue_d),max(nue_d),min(nue_d),min(nue_d),max(nue_d)],[max(wue_d),min(wue_d),min(wue_d),max(wue_d),max(wue_d)],color='r', marker='d',label='dry meadow')     
#    axA.fill_between([np.min(nue_d),np.max(nue_d)],np.min(wue_d),np.max(wue_d),color='brown',alpha=0.5) 
##    axA.set_ylim([np.min(wue_d)-np.min(wue_d)*0.01,np.max(wue_d)+np.max(wue_d)*0.01])
##    axA.set_xlim([np.min(nue_d)-0.5,np.max(nue_d)+0.5])
#    
#    axA.plot([max(nue_m),max(nue_m),min(nue_m),min(nue_m),max(nue_m)],[max(wue_m),min(wue_m),min(wue_m),max(wue_m),max(wue_m)],color='orange', linestyle='-',linewidth=4.5, label='moist meadow')
##    axA.scatter([max(nue_m),max(nue_m),min(nue_m),min(nue_m),max(nue_m)],[max(wue_m),min(wue_m),min(wue_m),max(wue_m),max(wue_m)],color='y', marker='o',s=40,label='moist meadow') 
#    axA.fill_between([np.min(nue_m),np.max(nue_m)],np.min(wue_m),np.max(wue_m),color='orange',alpha=0.5) 
##    axA.set_ylim([np.min(wue_m)-np.min(wue_m)*0.01,np.max(wue_m)+np.max(wue_m)*0.01])
##    axA.set_xlim([np.min(nue_m)-0.5,np.max(nue_m)+0.5])   
#    
#    axA.plot([max(nue_w),max(nue_w),min(nue_w),min(nue_w),max(nue_w)],[max(wue_w),min(wue_w),min(wue_w),max(wue_w),max(wue_w)],color='g', linestyle='-',linewidth=2.0,label='wet meadow')
##    axA.scatter([max(nue_w),max(nue_w),min(nue_w),min(nue_w),max(nue_w)],[max(wue_w),min(wue_w),min(wue_w),max(wue_w),max(wue_w)],facecolor='none',edgecolor='g', marker='o',s=120,label='wet meadow')    
#    axA.fill_between([np.min(nue_w),np.max(nue_w)],np.min(wue_w),np.max(wue_w),color='green',alpha=0.5) 
##    axA.set_ylim([np.min(wue_w)-np.min(wue_w)*0.01,np.max(wue_w)+np.max(wue_w)*0.01])
##    axA.set_xlim([np.min(nue_w)-0.5,np.max(nue_w)+0.5])    
#    
#    axA.plot([max(nue_s),max(nue_s),min(nue_s),min(nue_s),max(nue_s)],[max(wue_s),min(wue_s),min(wue_s),max(wue_s),max(wue_s)],color='c', linestyle='-',linewidth=1.0, label='snowbed')
##    axA.scatter([max(nue_s),max(nue_s),min(nue_s),min(nue_s),max(nue_s)],[max(wue_s),min(wue_s),min(wue_s),max(wue_s),max(wue_s)],color='b', marker='*', label='snowbed')         
#    axA.fill_between([np.min(nue_s),np.max(nue_s)],np.min(wue_s),np.max(wue_s),color='cyan',alpha=0.5) 
##    axA.set_ylim([np.min(wue_s)-np.min(wue_s)*0.01,np.max(wue_s)+np.max(wue_s)*0.01])
##    axA.set_xlim([np.min(nue_s)-0.5,np.max(nue_s)+0.5])    
   
    
    ##---Separate Plots Into Different Figures: Legend---##
    #ax1.legend(loc=4)
  

##---------------Box Plot Plant communities vs. Stomatal Condcutance--------------- #     

    gsw_bp=axB.boxplot([gsw_f,gsw_d,gsw_m,gsw_w,gsw_s], patch_artist=True, showmeans=True, showfliers=False)
  
    axB.set_xticks([1, 2, 3,4,5])
    axB.set_xticklabels(['Fellfield','Dry Meadow','Moist Meadow','Wet Meadow','Snowbed'],fontname='Times New Roman')
    axB.plot([1,2,3,4,5],[np.mean(gsw_f),np.mean(gsw_d),np.mean(gsw_m),np.mean(gsw_w),np.mean(gsw_s)],'-c')
    
    axB.tick_params(axis='x', labelsize=25)
    
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

 ##---------------Box Plot Plant communities vs. Leaf N--------------- #     

    na_bp=axC.boxplot([na_f,na_d,na_m,na_w,na_s], patch_artist=True, showmeans=True, showfliers=False)
  
    axC.set_xticks([1, 2, 3,4,5])
    axC.set_xticklabels(['Fellfield','Dry Meadow','Moist Meadow','Wet Meadow','Snowbed'],fontname='Times New Roman')
    axC.plot([1,2,3,4,5],[np.mean(na_f),np.mean(na_d),np.mean(na_m),np.mean(na_w),np.mean(na_s)],'-',color='orange')
  
    axC.tick_params(axis='x', labelsize=25)
    
    #gsw boxplot specs
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

        
 ##---------------Box Plot Plant communities vs. Assimilation--------------- #     

    A_bp=axD.boxplot([A_f,A_d,A_m,A_w,A_s], patch_artist=True, showmeans=True, showfliers=False)
  
    axD.set_xticks([1, 2, 3,4,5])
    axD.set_xticklabels(['Fellfield','Dry Meadow','Moist Meadow','Wet Meadow','Snowbed'],fontname='Times New Roman')
    axD.plot([1,2,3,4,5],[np.mean(A_f),np.mean(A_d),np.mean(A_m),np.mean(A_w),np.mean(A_s)],'-',color='purple')
  
    axD.tick_params(axis='x', labelsize=25)
    
    #gsw boxplot specs
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
        
 ##---------------Box Plot Plant communities vs. Evapotranspiration--------------- #     

    E_bp=axE.boxplot([E_f,E_d,E_m,E_w,E_s], patch_artist=True, showmeans=True, showfliers=False)
  
    axE.set_xticks([1, 2, 3,4,5])
    axE.set_xticklabels(['Fellfield','Dry Meadow','Moist Meadow','Wet Meadow','Snowbed'],fontname='Times New Roman')
    axE.plot([1,2,3,4,5],[np.mean(E_f),np.mean(E_d),np.mean(E_m),np.mean(E_w),np.mean(E_s)],'-k')
  
    axE.tick_params(axis='x', labelsize=25)
    
    
    #gsw boxplot specs
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
        
        
##---------------Regression Plot WUE vs. NUE--------------- #   

    axF.plot([np.mean(nue_f),np.mean(nue_d),np.mean(nue_m),np.mean(nue_w),np.mean(nue_s)],[np.mean(wue_f),np.mean(wue_d),np.mean(wue_m),np.mean(wue_w),np.mean(wue_s)])
    axF.scatter(nue_tot,wue_tot) 

    
    #---------------Make Plot Interactive---------------# 
#    
#        plt.pause(0.0001)
#        plt.ion()
    #end of sensitivity analysis iterations
    
#---------------Finalize Figure---------------#    

    #axA refers to first figure in subplot; axB refers to second figure in subplot
    #if only one axis is run then the figure is just one plot

    ##---Legend---##
#    axA.legend(bbox_to_anchor=(1, 1), loc='left', prop={'size':15})

    ##---Save Figure--##
    fig1.savefig('NUE_vs_WUE_var_allcoms.png') 
    fig2.savefig('gsw_plantcoms.png') 
    fig3.savefig('na_plantcoms.png') 
    fig4.savefig('A_plantcoms.png')
    fig5.savefig('E_plantcoms.png')
    fig6.savefig('WUE_vs_NUE_regression.png')
