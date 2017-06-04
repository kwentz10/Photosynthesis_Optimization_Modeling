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
import operator
from matplotlib import pyplot as plt
from matplotlib import rcParams
from scipy.optimize import curve_fit

#Import combinations of variable parameters 
from uncertain_params import params

#Import photosynthesis model
from Photosynthesis_Model import photo_bound_meso_eqstom as photo

#Import functions to switch between Pa and umol/mol at sea level
from photo_functions import pa_con_atmfrac


#---------------Determine if I Want to Keep Any of the Variable Parameters Constant---------------#

const_params=[]
for xx in it.combinations(['ht','t'],0): #keep ht and t constant for constant vpd
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
    
    axA.set_ylabel('NUE (umol CO2/g N s)',fontsize=28, fontname='Times New Roman')
    axA.set_ylim([0,5])
    axA.set_title('NUE in Alpine Tundra Plant Communities', fontname='Times New Roman',fontsize=36,fontweight='bold')

    
    #--figure 2--#
    
    #put in correct ax value (e.g. axA, axB)
    fig2,axB = plt.subplots(figsize=(15,15))

    #axA.set_xlabel('Plant Communities',fontsize=20, fontname='Times New Roman')
    axB.set_ylabel('gsw (mol H2O/m2s)',fontsize=28, fontname='Times New Roman')
    #axA.set_xlim([0,20])
    #axA.set_ylim([0,10])
    axB.set_title('Plant Communities vs. Stomatal Conductance', fontname='Times New Roman',fontsize=36,fontweight='bold')


    #--figure 4--#
    
    #put in correct ax value (e.g. axA, axB)
    fig4,axD = plt.subplots(figsize=(15,15))

    #axA.set_xlabel('Plant Communities',fontsize=20, fontname='Times New Roman')
    axD.set_ylabel('Assimilation (umol CO2/m2s)',fontsize=28, fontname='Times New Roman')
    #axA.set_xlim([0,20])
    #axA.set_ylim([0,10])
    axD.set_title('Assimilation in Alpine Tundra Plant Communities', fontname='Times New Roman',fontsize=36,fontweight='bold')

  
    #--figure 5--#
    
    #put in correct ax value (e.g. axA, axB)
    fig5,axE = plt.subplots(figsize=(15,15))

    #axA.set_xlabel('Plant Communities',fontsize=20, fontname='Times New Roman')
    axE.set_ylabel('Evapotranspiration (mmol H2O/m2s)',fontsize=28, fontname='Times New Roman')
    #axA.set_xlim([0,20])
    #axA.set_ylim([0,10])
    axE.set_title('Plant Communities vs. Evapotranspiration', fontname='Times New Roman',fontsize=36,fontweight='bold')

   
    #--figure 6--#
    
    #put in correct ax value (e.g. axA, axB)
    fig6,axF = plt.subplots(figsize=(15,15))

    axF.set_xlabel('NUE (umol CO2/g N s)',fontsize=28, fontname='Times New Roman')
    axF.set_ylabel('WUE (umol CO2/mmol H2O)',fontsize=36, fontname='Times New Roman')
    axF.set_xlim([0,5])
    axF.set_ylim([0,5])
    axF.set_title('NUE vs. WUE for all Plant Communities', fontname='Times New Roman',fontsize=36,fontweight='bold')
  
    #--figure 8--#
    fig8,axG = plt.subplots(figsize=(15,15))

    axG.set_xlabel('VPD (cmol H2O/mol air)',fontsize=28, fontname='Times New Roman')
   
    axG.set_ylabel('WUE (umol CO2/mmol H2O)',fontsize=28, fontname='Times New Roman')
#    axG.set_xlim([0,9])
#    axG.set_ylim([0,3])
    axG.set_title('VPD vs. WUE in Alpine Tundra Plant Communities', fontname='Times New Roman',fontsize=36,fontweight='bold')

    
    #-----figure 7----#
    
    fig7,axH = plt.subplots(figsize=(30,15))
    
    axH.set_ylabel('WUE (umol CO2/mmol H2O)',fontsize=28, fontname='Times New Roman')
    axH.set_ylim([0,5])
    axH.set_title('WUE in Alpine Tundra Plant Communities', fontname='Times New Roman',fontsize=36,fontweight='bold')



# 
#     #-----figure 9----#
#         
#    #put in correct ax value (e.g. axA, axB)
#    fig9,axI = plt.subplots(figsize=(30,15))
#    
#    #twin axis
##    axI2=axI.twinx()
#
#    #axA.set_xlabel('Plant Communities',fontsize=20, fontname='Times New Roman')
#    axI.set_ylabel('NUE (umol CO2/g N s)',fontsize=36, fontname='Times New Roman')
##    axI2.set_ylabel('Extended Summer NUE (umol CO2/g N s)', fontsize=36, fontname='Times New Roman')  
#    axI.set_ylim([0,5])
##    axI2.set_ylim([0,4])
#    axI.set_title('NUE Variability in an Extended Summer', fontname='Times New Roman',fontsize=36,fontweight='bold')
#
#     #-----figure 10----#
#         
#    #put in correct ax value (e.g. axA, axB)
#    fig10,axJ = plt.subplots(figsize=(30,15))
#    
##    #twin axis
##    axJ2=axJ.twinx()
#
#    #axA.set_xlabel('Plant Communities',fontsize=20, fontname='Times New Roman')
#    axJ.set_ylabel('WUE (umol CO2/mmol H2O)',fontsize=36, fontname='Times New Roman')
##    axJ2.set_ylabel('Extended Summer WUE (umol CO2/mmol H2O)', fontsize=36, fontname='Times New Roman')  
#    axJ.set_ylim([0,5])
##    axJ2.set_ylim([0,4])
#    axJ.set_title('WUE Variability in an Extended Summer', fontname='Times New Roman',fontsize=36,fontweight='bold')


        #-----figure 11----#

    #put in correct ax value (e.g. axA, axB)
    fig11,axS = plt.subplots(figsize=(15,15))

    axS.set_xlabel('Assimilation (umol CO2/m2s)',fontsize=28, fontname='Times New Roman')
    axS.set_ylabel('Stomatal Conductance (mol CO2/m2s)',fontsize=28, fontname='Times New Roman')
    axS.set_title('Assimilation vs. Stomatal Conductance in Alpine Tundra Plant Communities', fontname='Times New Roman',fontsize=30,fontweight='bold')
 
    #---------------Initialize Arrays for Each Meadow---------------#
        
    #total nue and wue
    nue_tot=[]
    wue_tot=[]
   
    #wue and nue arrays
    wue_d=[]
    nue_d=[]
    wue_d_tms=[]
    nue_d_tms=[]
    wue_m=[]
    nue_m=[]
    wue_m_tms=[]
    nue_m_tms=[]        
    wue_w=[]
    nue_w=[]    
    wue_w_tms=[]
    nue_w_tms=[]        

    
    
    #gsw arrays

    gsw_d=[]
    gsw_d_tms=[]        
    gsw_m=[]
    gsw_m_tms=[]
    gsw_w=[]
    gsw_w_tms=[]
    
    gs_d=[]
    gs_m=[]
    gs_w=[]    
    #assimilation arrays

    A_d=[]
    A_d_tms=[]
    A_m=[]
    A_m_tms=[]
    A_w=[]
    A_w_tms=[]

    #evapo arrays

    E_d=[]
    E_d_tms=[]        
    E_m=[]
    E_m_tms=[]
    E_w=[]
    E_w_tms=[]  
    
    #vapor pressure deficit arrays
    vpd_d=[]
    vpd_d_tms=[]
    vpd_m=[]
    vpd_m_tms=[]
    vpd_w=[]
    vpd_w_tms=[]    


    #leaf temp
    tl_d=[]
    tl_m=[]
    tl_w=[]




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
            if 't' in const_params[ii]:
                temp=t_c
    
                
    
            
            #------calculate vapor pressure-----#
            pa_v=611*np.exp((17.27*temp)/(temp+237.3)) #saturation vapor pressure of air (Pa)
            ea_str=pa_con_atmfrac(pa_v,3528) #saturation vapor pressure of air (Pa-->umol h20/mol air)
            ea=rh*ea_str #vapor pressure (umol h2O/mol air)                
    
    
            #correct for leaf temperatures using leaf height
       
            t_diff=18-0.4*ht
        
            tl=temp+t_diff
    
            
    
            
            if xx==0: 
    
                
                #---------------Photosynthesis Function---------------#
            
                #alter this line of code for when implementing different photosynthesis functions
                wue, nue, A, E, cs, ci, gsw, gs, gbw, gb, gm, cc,dd =photo(tk_25,ekc,eko,etau,ev,ej,toptv,toptj,na, qeff, PAR,tl,ea,chl,ij,kc25,ko25,o,ca,rh,m,a,frnr,flnr,ra,jm,g0,b,dia,u,q,vwc_min,vwc_max,vwc)
     
           
                #test to make sure wue and nue are positive at not 'nan'
                if wue[0]==-999 and nue[0]==-999:
               
                    continue
            
                
                    
                wue_d+=[wue[0]]
                nue_d+=[nue[0]]
                gsw_d+=[gsw[0]]
                A_d+=[A[0]]
                E_d+=[E[0]]
                vpd_d+=[dd[0]]
                tl_d+=[tl]
                gs_d+=[gs[0]]
                
    #            elif xx==1:
    #                wue_d_tms+=[wue.tolist()[0]]
    #                nue_d_tms+=[nue.tolist()[0]]
    #                gsw_d_tms+=[gsw.tolist()[0]]
    #                A_d_tms+=[A[0]]
    #                E_d_tms+=[E.tolist()[0]]
    #                vpd_d_tms+=[dd.tolist()[0]]
                    
                
            elif xx==1:
                
    
                #---------------Photosynthesis Function---------------#
            
                #alter this line of code for when implementing different photosynthesis functions
                wue, nue, A, E, cs, ci, gsw, gs, gbw, gb, gm, cc,dd =photo(tk_25,ekc,eko,etau,ev,ej,toptv,toptj,na, qeff, PAR,tl,ea,chl,ij,kc25,ko25,o,ca,rh,m,a,frnr,flnr,ra,jm,g0,b,dia,u,q,vwc_min,vwc_max,vwc)
     
           
                #test to make sure wue and nue are positive at not 'nan'
                if wue[0]==-999 and nue[0]==-999:
               
                    continue                
                
                
    
                    
                wue_m+=[wue[0]]
                nue_m+=[nue[0]]
                gsw_m+=[gsw[0]]
                A_m+=[A[0]]
                E_m+=[E[0]]
                vpd_m+=[dd[0]]
                tl_m+=[tl]
                gs_m+=[gs[0]]
                
    #            elif xx==1:
    #                wue_m_tms+=[wue.tolist()[0]]
    #                nue_m_tms+=[nue.tolist()[0]]
    #                gsw_m_tms+=[gsw.tolist()[0]]
    #                A_m_tms+=[A[0]]
    #                E_m_tms+=[E.tolist()[0]]
    #                vpd_m_tms+=[dd.tolist()[0]]
    
       
        
            elif xx==2:
    
                #---------------Photosynthesis Function---------------#
            
                #alter this line of code for when implementing different photosynthesis functions
                wue, nue, A, E, cs, ci, gsw, gs, gbw, gb, gm, cc,dd =photo(tk_25,ekc,eko,etau,ev,ej,toptv,toptj,na, qeff, PAR,tl,ea,chl,ij,kc25,ko25,o,ca,rh,m,a,frnr,flnr,ra,jm,g0,b,dia,u,q,vwc_min,vwc_max,vwc)
     
           
                #test to make sure wue and nue are positive at not 'nan'
                if wue[0]==-999 and nue[0]==-999:
               
                    continue                                
                
    
               
                wue_w+=[wue[0]]
                nue_w+=[nue[0]]
                gsw_w+=[gsw[0]]
                A_w+=[A[0]]
                E_w+=[E[0]]
                vpd_w+=[dd[0]]
                tl_w+=[tl]
                gs_w+=[gs[0]]
                
    #            elif xx==1:
    #                wue_w_tms+=[wue.tolist()[0]]
    #                nue_w_tms+=[nue.tolist()[0]]
    #                gsw_w_tms+=[gsw.tolist()[0]]
    #                A_w_tms+=[A[0]]
    #                E_w_tms+=[E.tolist()[0]]
    #                vpd_w_tms+=[dd.tolist()[0]]

    
#---------------Plot Plant Communities vs. NUE & WUE---------------#      

#        axA.bar([1.5,2.5,3.5],[np.mean(nue_d),np.mean(nue_m),np.mean(nue_w)], yerr=[np.std(nue_d),np.std(nue_m),np.std(nue_w)], edgecolor='black',linewidth=lw[xx], width=0.2, align="center",color='red',alpha=1.0,error_kw={'ecolor':'black', 'lw':2, 'capsize':5, 'capthick':2})
#        axA.set_xticks([1.5, 2.5, 3.5])
#        axA.set_xticklabels(['Dry Meadow','Moist Meadow','Wet Meadow'],fontname='Times New Roman')
#        axA.tick_params(axis='x', labelsize=28)
#        axA.tick_params(axis='y', labelsize=18)
# 
#        axH.bar([1.5,2.5,3.5],[np.mean(wue_d),np.mean(wue_m),np.mean(wue_w)], yerr=[np.std(wue_d),np.std(wue_m),np.std(wue_w)], edgecolor='black',linewidth=lw[xx], width=0.2, align="center",color='blue',alpha=1.0,error_kw={'ecolor':'black', 'lw':2, 'capsize':5, 'capthick':2})      
#        axH.set_xticks([1.5, 2.5, 3.5])
#        axH.set_xticklabels(['Dry Meadow','Moist Meadow','Wet Meadow'],fontname='Times New Roman')
#        axH.tick_params(axis='x', labelsize=28)
#        axH.tick_params(axis='y', labelsize=18)
#  
#
#        axI.bar([0.8,1.8,2.8],[np.mean(nue_d),np.mean(nue_m),np.mean(nue_w)], yerr=[np.std(nue_d),np.std(nue_m),np.std(nue_w)], edgecolor='black',linewidth=lw[xx], width=0.2, color='red',alpha=1.0,error_kw={'ecolor':'black', 'lw':2, 'capsize':5, 'capthick':2})
#        axI2.bar([1,2,3],[np.mean(wue_d),np.mean(wue_m),np.mean(wue_w)], yerr=[np.std(wue_d),np.std(wue_m),np.std(wue_w)], edgecolor='black',linewidth=lw[xx], width=0.2, color='blue',alpha=1.0,error_kw={'ecolor':'black', 'lw':2, 'capsize':5, 'capthick':2})
#        
#        axI.set_xticks([1, 2, 3])
#        axI.set_xticklabels(['Dry Meadow','Moist Meadow','Wet Meadow'],fontname='Times New Roman')
#        axI.tick_params(axis='x', labelsize=28)
#        axI.tick_params(axis='y', labelsize=18)
#        axI2.tick_params(axis='y', labelsize=18)



    
#    #nue tms
#    nue_0_bp=axI.boxplot([nue_d,nue_m,nue_w], positions=[0.875,1.875,2.875],widths=0.25, patch_artist=True, showmeans=True, showfliers=False)
#    nue_tms_bp=axI.boxplot([nue_d_tms,nue_m_tms,nue_w_tms], positions=[1.125,2.125,3.125],widths=0.25, patch_artist=True, showmeans=True, showfliers=False)
#
#    axI.set_xticks([1,2,3])
#    axI.set_xticklabels(['Dry Meadow','Moist Meadow','Wet Meadow'],fontname='Times New Roman')
#    rcParams['xtick.labelsize']=24   
#
# 
#    #nue initial boxplot specs
#    for box in nue_0_bp['boxes']:
#        #change outline color
#        box.set(color='black',linewidth=2)
#        #change fill color
#        box.set(facecolor='black',alpha=0.2)
#
#    for whisker in nue_0_bp['whiskers']:
#        whisker.set(color='black',linewidth=2)
#    
#    for cap in nue_0_bp['caps']:
#        cap.set(color='black',linewidth=2)
#
#    for median in nue_0_bp['medians']:
#        median.set(color='black', linewidth=2)
#
#    for flier in nue_0_bp['fliers']:
#        flier.set(marker='*',color='black',alpha=0.5)
#
#    for means in nue_0_bp['means']:
#        means.set(marker='o',markerfacecolor='black')    
#   
#    #nue tms boxplot specs    
#    for box in nue_tms_bp['boxes']:
#        #change outline color
#        box.set(color='black',linewidth=2)
#        #change fill color
#        box.set(facecolor='black',alpha=0.5)
#
#    for whisker in nue_tms_bp['whiskers']:
#        whisker.set(color='black',linewidth=5)
#    
#    for cap in nue_tms_bp['caps']:
#        cap.set(color='black',linewidth=5)
#
#    for median in nue_tms_bp['medians']:
#        median.set(color='black', linewidth=5)
#
#    for flier in nue_tms_bp['fliers']:
#        flier.set(marker='*',color='black',alpha=0.5)
#
#    for means in nue_tms_bp['means']:
#        means.set(marker='o',markerfacecolor='black')
#
#    #wue tms
#    wue_0_bp=axJ.boxplot([wue_d,wue_m,wue_w], positions=[0.875,1.875,2.875],widths=0.25, patch_artist=True, showmeans=True, showfliers=False)
#    wue_tms_bp=axJ.boxplot([wue_d_tms,wue_m_tms,wue_w_tms], positions=[1.125,2.125,3.125],widths=0.25, patch_artist=True, showmeans=True, showfliers=False)
#
#
#    axJ.set_xticks([1,2,3])
#    axJ.set_xticklabels(['Dry Meadow','Moist Meadow','Wet Meadow'],fontname='Times New Roman')
#    rcParams['xtick.labelsize']=24   
#    
#    
#    #wue initial boxplot specs
#    for box in wue_0_bp['boxes']:
#        #change outline color
#        box.set(color='black',linewidth=2)
#        #change fill color
#        box.set(facecolor='black',alpha=0.2)
#
#    for whisker in wue_0_bp['whiskers']:
#        whisker.set(color='black',linewidth=2)
#    
#    for cap in wue_0_bp['caps']:
#        cap.set(color='black',linewidth=2)
#
#    for median in wue_0_bp['medians']:
#        median.set(color='black', linewidth=2)
#
#    for flier in wue_0_bp['fliers']:
#        flier.set(marker='*',color='black',alpha=0.5)
#
#    for means in wue_0_bp['means']:
#        means.set(marker='o',markerfacecolor='black')    
#   
#    #nue tms boxplot specs    
#    for box in wue_tms_bp['boxes']:
#        #change outline color
#        box.set(color='black',linewidth=2)
#        #change fill color
#        box.set(facecolor='black',alpha=0.5)
#
#    for whisker in wue_tms_bp['whiskers']:
#        whisker.set(color='black',linewidth=5)
#    
#    for cap in wue_tms_bp['caps']:
#        cap.set(color='black',linewidth=5)
#
#    for median in wue_tms_bp['medians']:
#        median.set(color='black', linewidth=5)
#
#    for flier in wue_tms_bp['fliers']:
#        flier.set(marker='*',color='black',alpha=0.5)
#
#    for means in wue_tms_bp['means']:
#        means.set(marker='o',markerfacecolor='black')
#   


    #just nue
    nue_bp=axA.boxplot([nue_d,nue_m,nue_w], patch_artist=True, showmeans=True, showfliers=False)
    
    axA.set_xticks([1,2,3])
    axA.set_xticklabels(['Dry Meadow','Moist Meadow','Wet Meadow'],fontname='Times New Roman')
    rcParams['xtick.labelsize']=28   

    
    #nue boxplot specs
    for box in nue_bp['boxes']:
        #change outline color
        box.set(color='black',linewidth=2)
        #change fill color
        box.set(facecolor='black',alpha=0.2)

    for whisker in nue_bp['whiskers']:
        whisker.set(color='black',linewidth=2)
    
    for cap in nue_bp['caps']:
        cap.set(color='black',linewidth=2)

    for median in nue_bp['medians']:
        median.set(color='black', linewidth=2)

    for flier in nue_bp['fliers']:
        flier.set(marker='*',color='black',alpha=0.5)

    for means in nue_bp['means']:
        means.set(marker='o',markerfacecolor='black')    

    #just wue
    wue_bp=axH.boxplot([wue_d,wue_m,wue_w], patch_artist=True, showmeans=True, showfliers=False)
    
    axH.set_xticks([1,2,3])
    axH.set_xticklabels(['Dry Meadow','Moist Meadow','Wet Meadow'],fontname='Times New Roman')
    rcParams['xtick.labelsize']=28   

 
    #wue boxplot specs    
    for box in wue_bp['boxes']:
        #change outline color
        box.set(color='black',linewidth=2)
        #change fill color
        box.set(facecolor='black',alpha=0.2)

    for whisker in wue_bp['whiskers']:
        whisker.set(color='black',linewidth=2)
    
    for cap in wue_bp['caps']:
        cap.set(color='black',linewidth=2)

    for median in wue_bp['medians']:
        median.set(color='black', linewidth=2)

    for flier in wue_bp['fliers']:
        flier.set(marker='*',color='black',alpha=0.5)

    for means in wue_bp['means']:
        means.set(marker='o',markerfacecolor='black')
      

#---------------Box Plot Plant Communities vs. Stomatal Condcutance--------------- #     
#    axB.plot([1,2,3,4,5],[np.mean(gsw_f),np.mean(gsw_d),np.mean(gsw_m),np.mean(gsw_w),np.mean(gsw_s)],'-c')
    
    
#        axB.bar([1.5,2.5,3.5],[np.mean(gsw_d),np.mean(gsw_m),np.mean(gsw_w)], yerr=[np.std(gsw_d),np.std(gsw_m),np.std(gsw_w)],edgecolor='black',linewidth=lw[xx], width=0.2, align="center",color='cyan',alpha=1.0,error_kw={'ecolor':'black', 'lw':2, 'capsize':5, 'capthick':2})
#        axB.set_xticks([1.5, 2.5, 3.5])
#        axB.set_xticklabels(['Dry Meadow','Moist Meadow','Wet Meadow'],fontname='Times New Roman')
#        axB.tick_params(axis='x', labelsize=25)
#        axB.tick_params(axis='y', labelsize=16)

    gsw_bp=axB.boxplot([gsw_d,gsw_m,gsw_w], patch_artist=True, showmeans=True, showfliers=False)

    axB.set_xticks([1,2,3])
    axB.set_xticklabels(['Dry Meadow','Moist Meadow','Wet Meadow'],fontname='Times New Roman')
    rcParams['xtick.labelsize']=28   
    
    #gsw boxplot specs
    for box in gsw_bp['boxes']:
        #change outline color
        box.set(color='black',linewidth=2)
        #change fill color
        box.set(facecolor='black',alpha=0.2)

    for whisker in gsw_bp['whiskers']:
        whisker.set(color='black',linewidth=2)
    
    for cap in gsw_bp['caps']:
        cap.set(color='black',linewidth=2)

    for median in gsw_bp['medians']:
        median.set(color='black', linewidth=2)

    for flier in gsw_bp['fliers']:
        flier.set(marker='*',color='black',alpha=0.5)

    for means in gsw_bp['means']:
        means.set(marker='o',markerfacecolor='black')    
  
       
#---------------Box Plot Plant Communities vs. Assimilation--------------- #     

    A_bp=axD.boxplot([A_d,A_m,A_w], patch_artist=True, showmeans=True, showfliers=False)

    axD.set_xticks([1,2,3])
    axD.set_xticklabels(['Dry Meadow','Moist Meadow','Wet Meadow'],fontname='Times New Roman')
    rcParams['xtick.labelsize']=28   

#        axD.plot([1,2,3,4,5],[np.mean(A_f),np.mean(A_d),np.mean(A_m),np.mean(A_w),np.mean(A_s)],'-',color='purple')
    
#        axD.bar([1.5,2.5,3.5],[np.mean(A_d),np.mean(A_m),np.mean(A_w)],yerr=[np.std(A_d),np.std(A_m),np.std(A_w)],edgecolor='black',linewidth=lw[xx], width=0.2, align="center",color='purple',alpha=1.0,error_kw={'ecolor':'black', 'lw':2, 'capsize':5, 'capthick':2})
#        axD.set_xticks([1.5, 2.5, 3.5])
#        axD.set_xticklabels(['Dry Meadow','Moist Meadow','Wet Meadow'],fontname='Times New Roman')
#        axD.tick_params(axis='x', labelsize=25)
#        axD.tick_params(axis='y', labelsize=16)
    
    #A boxplot specs
    for box in A_bp['boxes']:
        #change outline color
        box.set(color='black',linewidth=2)
        #change fill color
        box.set(facecolor='black',alpha=0.2)

    for whisker in A_bp['whiskers']:
        whisker.set(color='black',linewidth=2)
    
    for cap in A_bp['caps']:
        cap.set(color='black',linewidth=2)

    for median in A_bp['medians']:
        median.set(color='black', linewidth=2)

    for flier in A_bp['fliers']:
        flier.set(marker='*',color='black',alpha=0.5)

    for means in A_bp['means']:
        means.set(marker='o',markerfacecolor='black')  
        
#---------------Box Plot Plant Communities vs. Evapotranspiration--------------- #     

    E_bp=axE.boxplot([E_d,E_m,E_w], patch_artist=True, showmeans=True, showfliers=False)

    axE.set_xticks([1,2,3])
    axE.set_xticklabels(['Dry Meadow','Moist Meadow','Wet Meadow'],fontname='Times New Roman')
    rcParams['xtick.labelsize']=28   

#    axE.plot([1,2,3,4,5],[np.mean(E_f),np.mean(E_d),np.mean(E_m),np.mean(E_w),np.mean(E_s)],'-k')
    
#        axE.bar([1.5,2.5,3.5],[np.mean(E_d)/1000.0,np.mean(E_m)/1000.0,np.mean(E_w)/1000.0],yerr=[np.std(E_d)/1000.0,np.std(E_m)/1000.0,np.std(E_w)/1000.0],edgecolor='black',linewidth=lw[xx], width=0.2,align="center",color='green',alpha=1.0,error_kw={'ecolor':'black', 'lw':2, 'capsize':5, 'capthick':2})
#        axE.set_xticks([1.5,2.5, 3.5])
#        axE.set_xticklabels(['Dry Meadow','Moist Meadow','Wet Meadow'],fontname='Times New Roman')
#        axE.tick_params(axis='x', labelsize=25)
#        axE.tick_params(axis='y', labelsize=16)
#        
    
    #E boxplot specs
    for box in E_bp['boxes']:
        #change outline color
        box.set(color='black',linewidth=2)
        #change fill color
        box.set(facecolor='black',alpha=0.2)

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
   
    
    axF.scatter([np.mean(nue_d),np.mean(nue_m),np.mean(nue_w)],[np.mean(wue_d),np.mean(wue_m),np.mean(wue_w)],edgecolors='black',facecolors='black',marker='o',s=50)
    axF.plot([np.mean(nue_d),np.mean(nue_m),np.mean(nue_w)],[np.mean(wue_d),np.mean(wue_m),np.mean(wue_w)],color='black')
    axF.arrow(np.mean(nue_m),np.mean(wue_m),np.mean(nue_w)-np.mean(nue_m),np.mean(wue_w)-np.mean(wue_m), head_width=0.15, head_length=0.15, fc='black', ec='black')
#    axF.scatter(nue_tot,wue_tot) 
    
    #axF.plot(np.unique([np.mean(nue_d),np.mean(nue_m),np.mean(nue_w)]), np.poly1d(np.polyfit([np.mean(nue_d),np.mean(nue_m),np.mean(nue_w)], [np.mean(wue_d),np.mean(wue_m),np.mean(wue_w)], 1))(np.unique([np.mean(nue_d),np.mean(nue_m),np.mean(nue_w)])),color=color[iii])
    #axF.plot([np.mean(nue_f),np.mean(nue_d),np.mean(nue_m),np.mean(nue_w),np.mean(nue_s)],[np.mean(wue_f),np.mean(wue_d),np.mean(wue_m),np.mean(wue_w),np.mean(wue_s)])

    #axF.tick_params(axis='x', labelsize=25)
#    axF.tick_params(axis='y', labelsize=16)
    

#---------------Regression Plot Assimilation vs. Stomatal Conductance---------------#   
   
    
    axS.scatter(A_d+A_m+A_w,gs_d+gs_m+gs_w,edgecolors='black',facecolors='black',marker='o',s=50)
#    axF.plot([np.mean(nue_d),np.mean(nue_m),np.mean(nue_w)],[np.mean(wue_d),np.mean(wue_m),np.mean(wue_w)],color='black')
    
    axS.plot(np.unique(A_d+A_m+A_w), np.poly1d(np.polyfit(A_d+A_m+A_w, gs_d+gs_m+gs_w, 1))(np.unique(A_d+A_m+A_w)),color='black',linewidth=6)
    #axF.plot([np.mean(nue_f),np.mean(nue_d),np.mean(nue_m),np.mean(nue_w),np.mean(nue_s)],[np.mean(wue_f),np.mean(wue_d),np.mean(wue_m),np.mean(wue_w),np.mean(wue_s)])

    axS.set_ylim([0,0.55])
    axS.set_xlim([0,25])
#    axS.tick_params(axis='x', labelsize=16)
#    axS.tick_params(axis='y', labelsize=16)



    
#---------------Plot VPD vs. WUE for validation (use all points rather than mean of points)-------------#
    def func(x, a, b, c):
        return a * np.exp(-b * x) + c
    xdata=np.array(vpd_d+vpd_m+vpd_w)/10000.
    ydata=np.array(wue_d+wue_m+wue_w)
    L=sorted(zip(xdata,ydata),key=operator.itemgetter(0))
    new_x,new_y=zip(*L)
    popt, pcov = curve_fit(func, np.array(new_x), np.array(new_y))
    
    axG.plot(np.array(new_x), func(np.array(new_x), *popt), 'k-', label='fit',linewidth=6)
    axG.scatter(xdata,ydata,edgecolors='black',facecolors='black',marker='o',s=50)

#    #axG.plot([np.mean(vpd_d),np.mean(vpd_m),np.mean(vpd_w)],[np.mean(wue_d),np.mean(wue_m),np.mean(wue_w)],color=color[xx])
#    axG.plot(np.unique(vpd_d+vpd_m+vpd_w), np.poly1d(np.polyfit(vpd_d+vpd_m+vpd_w, wue_d+wue_m+wue_w, 1))(np.unique(vpd_d+vpd_m+vpd_w)),color='black')
#    axG.legend()
    
    #validation data:
#    xdata=np.array([13734.61196533664,13734.61196533664,13734.61196533664,13734.61196533664,21364.95194607922,21364.95194607922,21364.95194607922,21364.95194607922,28995.2919268218,28995.2919268218,28995.2919268218,28995.2919268218,13734.61196533664,13734.61196533664,13734.61196533664,13734.61196533664,21364.95194607922,21364.95194607922,21364.95194607922,21364.95194607922,28995.2919268218,28995.2919268218,28995.2919268218,28995.2919268218])/10000.
#    ydata=np.array([11,12,11,3,5,4.5,1,5,4,2,1,.5,9,9,10,11,5,5,5,5,5,5,5,2])*.4091
#    axG.scatter(xdata,ydata,edgecolors="black",facecolors="black",marker='o',s=50)   
##    #---------------Make Plot Interactive---------------# 
#       
#        plt.pause(0.00000001)
#        plt.ion()

    axA.tick_params(axis='y', labelsize=20)
#    axA.yaxis.label.set_color('red')



    axB.tick_params(axis='y', labelsize=20)
#    axB.yaxis.label.set_color('red')  


    axE.tick_params(axis='y', labelsize=20)
#    axC.yaxis.label.set_color('red')  
 

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



    axE.tick_params(axis='y', labelsize=20)
#    axB.yaxis.label.set_color('red')  


    axF.tick_params(axis='y', labelsize=20)
#    axC.yaxis.label.set_color('red')  


    axG.tick_params(axis='y', labelsize=20)
#    axC.yaxis.label.set_color('red')
    axG.tick_params(axis='x', pad=15,labelsize=20)


    axH.tick_params(axis='y', labelsize=20)
#    axB.yaxis.label.set_color('red')  


    axS.tick_params(axis='y', labelsize=20)
#    axC.yaxis.label.set_color('red')  
    axS.tick_params(axis='x', pad=15,labelsize= 20)



#---------------Finalize Figure---------------#    

    ##---Legend---##
    #axA.legend(bbox_to_anchor=(1, 1), loc='left', prop={'size':15})

    ##---Save Figure--##
    fig1.savefig('NUE.png') 
    fig7.savefig('WUE.png')
    fig4.savefig('Assimilation.png')
#    fig6.savefig('NUE_vs_WUE.png')
    fig8.savefig('VPD_vs_WUE.png')
#    fig9.savefig('NUE_tms.png')
#    fig10.savefig('WUE_tms.png')
    fig11.savefig('A_vs_gsw.png')
