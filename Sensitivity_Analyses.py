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



#Import photosynthesis model
from Photosynthesis_Model import photo_bound_meso_eqstom as photo

#Import functions to switch between Pa and umol/mol at sea level
from photo_functions import pa_con_atmfrac



#---------------Determine if I Want to Keep Any of the Variable Parameters Constant---------------#

var_params=[]
for xy in it.combinations(['dia','chl','na','ht'],4):
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
#    axA.set_xlabel('Surface Temperature ($^\circ$C)',fontsize=36, fontname='Times New Roman')
    axA.set_xlabel('Total Abiotic Variability (%)',fontsize=36, fontname='Times New Roman')
#    axA.set_xlabel('Leaf Chlorophyll (umol chl/m2)',fontsize=36, fontname='Times New Roman')
#    axA.set_xlabel('Leaf Diameter (m2)',fontsize=36, fontname='Times New Roman')
#    axA.set_xlabel('Leaf Height (cm)',fontsize=36, fontname='Times New Roman')
#    axA.set_xlabel('Leaf Nitrogen (gN/m2)',fontsize=36, fontname='Times New Roman')
#    axA.set_xlabel('Total Physiological Variability (%)',fontsize=36, fontname='Times New Roman')


    axA.set_ylabel('NUE (umol CO2/g N s)',fontsize=36, fontname='Times New Roman')
    axA.set_ylim([0,4])
    axA.set_title('NUE & WUE in the Dry Meadow', fontname='Times New Roman',fontsize=36,fontweight='bold')
    #twin axis
    axA2=axA.twinx()
    axA2.set_ylabel('WUE (umol CO2/mmol H2O)', fontsize=36, fontname='Times New Roman')  
    axA2.set_ylim([0,4])


    #--figure 2--#
    fig2,axB = plt.subplots(figsize=(15,15))

#    axB.set_xlabel('Volumetric Water Content (m3/m3)',fontsize=36, fontname='Times New Roman')
#    axB.set_xlabel('Surface Temperature ($^\circ$C)',fontsize=36, fontname='Times New Roman')
    axB.set_xlabel('Total Abiotic Variability (%)',fontsize=36, fontname='Times New Roman')
#    axB.set_xlabel('Leaf Chlorophyll (umol chl/m2)',fontsize=36, fontname='Times New Roman')
#    axB.set_xlabel('Leaf Diameter (m2)',fontsize=36, fontname='Times New Roman')
#    axB.set_xlabel('Leaf Height (cm)',fontsize=36, fontname='Times New Roman')
#    axB.set_xlabel('Leaf Nitrogen (gN/m2)',fontsize=36, fontname='Times New Roman')
#    axB.set_xlabel('Total Physiological Variability (%)',fontsize=36, fontname='Times New Roman')


    axB.set_ylabel('NUE (umol CO2/g N s)',fontsize=36, fontname='Times New Roman')
    axB.set_ylim([0,4])
    axB.set_title('NUE & WUE in the Moist Meadow', fontname='Times New Roman',fontsize=36,fontweight='bold')
    #twin axis
    axB2=axB.twinx()
    axB2.set_ylabel('WUE (umol CO2/mmol H2O)', fontsize=36, fontname='Times New Roman')  
    axB2.set_ylim([0,4])

    #--figure 3--#

    fig3,axC = plt.subplots(figsize=(15,15))

#    axC.set_xlabel('Volumetric Water Content (m3/m3)',fontsize=36, fontname='Times New Roman')
#    axC.set_xlabel('Surface Temperature ($^\circ$C)',fontsize=36, fontname='Times New Roman')
    axC.set_xlabel('Total Abiotic Variability (%)',fontsize=36, fontname='Times New Roman')
#    axC.set_xlabel('Leaf Chlorophyll (umol chl/m2)',fontsize=36, fontname='Times New Roman')
#    axC.set_xlabel('Leaf Diameter (m)',fontsize=36, fontname='Times New Roman')
#    axC.set_xlabel('Leaf Height (cm)',fontsize=36, fontname='Times New Roman')
#    axC.set_xlabel('Leaf Nitrogen (gN/m2)',fontsize=36, fontname='Times New Roman')
#    axC.set_xlabel('Total Physiological Variability (%)',fontsize=36, fontname='Times New Roman')
    
    axC.set_ylabel('NUE (umol CO2/g N s)',fontsize=36, fontname='Times New Roman')
    axC.set_ylim([0,4])
    axC.set_title('NUE & WUE in the Wet Meadow', fontname='Times New Roman',fontsize=36,fontweight='bold')
    #twin axis
    axC2=axC.twinx()
    axC2.set_ylabel('WUE (umol CO2/mmol H2O)', fontsize=36, fontname='Times New Roman')  
    axC2.set_ylim([0,4])    

#3D Plots
   
    fig4 = plt.figure()
    axD = fig4.add_subplot(111, projection='3d') 
    axD.set_title('Dry Meadow NUE As a Function of Abiotic and Physiological Variability', fontname='Times New Roman',fontsize=16,fontweight='bold')
    axD.set_xlabel('Total Abiotic Variability (%)',fontsize=12, fontname='Times New Roman')
    axD.set_ylabel('Total Physiological Variability (%)',fontsize=12, fontname='Times New Roman')
    axD.set_zlabel('NUE (umol CO2/g N s)',fontsize=12, fontname='Times New Roman')


    fig5 = plt.figure()
    axE = fig5.add_subplot(111, projection='3d') 
    axE.set_title('Dry Meadow WUE As a Function of Abiotic and Physiological Variability', fontname='Times New Roman',fontsize=16,fontweight='bold')
    axE.set_xlabel('Total Abiotic Variability (%)',fontsize=12, fontname='Times New Roman')
    axE.set_ylabel('Total Physiological Variability (%)',fontsize=12, fontname='Times New Roman')
    axE.set_zlabel('WUE (umol CO2/mmol H2O)',fontsize=12, fontname='Times New Roman')

    fig6 = plt.figure()
    axF = fig6.add_subplot(111, projection='3d') 
    axF.set_title('Moist Meadow NUE As a Function of Abiotic and Physiological Variability', fontname='Times New Roman',fontsize=16,fontweight='bold')
    axF.set_xlabel('Total Abiotic Variability (%)',fontsize=12, fontname='Times New Roman')
    axF.set_ylabel('Total Physiological Variability (%)',fontsize=12, fontname='Times New Roman')
    axF.set_zlabel('NUE (umol CO2/g N s)',fontsize=12, fontname='Times New Roman')

    fig7 = plt.figure()
    axG = fig7.add_subplot(111, projection='3d') 
    axG.set_title('Moist Meadow WUE As a Function of Abiotic and Physiological Variability', fontname='Times New Roman',fontsize=16,fontweight='bold')
    axG.set_xlabel('Total Abiotic Variability (%)',fontsize=12, fontname='Times New Roman')
    axG.set_ylabel('Total Physiological Variability (%)',fontsize=12, fontname='Times New Roman')
    axG.set_zlabel('WUE (umol CO2/mmol H2O)',fontsize=12, fontname='Times New Roman')

    fig8 = plt.figure()
    axH = fig8.add_subplot(111, projection='3d') 
    axH.set_title('Wet Meadow NUE As a Function of Abiotic and Physiological Variability', fontname='Times New Roman',fontsize=16,fontweight='bold')
    axH.set_xlabel('Total Abiotic Variability (%)',fontsize=12, fontname='Times New Roman')
    axH.set_ylabel('Total Physiological Variability (%)',fontsize=12, fontname='Times New Roman')
    axH.set_zlabel('NUE (umol CO2/g N s)',fontsize=12, fontname='Times New Roman')

    fig9 = plt.figure()
    axI = fig9.add_subplot(111, projection='3d') 
    axI.set_title('Wet Meadow WUE As a Function of Abiotic and Physiological Variability', fontname='Times New Roman',fontsize=16,fontweight='bold')
    axI.set_xlabel('Total Abiotic Variability (%)',fontsize=12, fontname='Times New Roman')
    axI.set_ylabel('Total Physiological Variability (%)',fontsize=12, fontname='Times New Roman')
    axI.set_zlabel('WUE (umol CO2/mmol H2O)',fontsize=12, fontname='Times New Roman')

    
    ####-----------run through various temperatures---------######

   
    #make 40% changes
    
    #abiotic  variability: temp and moisture
#    temp_per=np.linspace(1.0,1.4,10)
#    vwc_per=np.linspace(1.0,0.6,10)

##    #abiotic  variability: temp 
#    temp_per=np.linspace(1.0,1.4,10)
#    vwc_per=np.linspace(1.0,1.0,10)

##    #abiotic  variability: moisture
#    temp_per=np.linspace(1.0,1.0,10)
#    vwc_per=np.linspace(1.0,0.6,10)


#    abiotic and physiological variability (set all traits to varying)
#    temp_per=np.linspace(1.0,1.4,10)
#    vwc_per=np.linspace(1.0,0.6,10)
    
#    #physiological variability --remember to set plant traits as varying (for total--set all plant traits to varying)
    temp_per=np.linspace(1.0,1.0,1)
    vwc_per=np.linspace(1.0,1.0,1)

    abiotic_index=np.linspace(0,40,len(temp_per))
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
    abi_d_x=[]
 
    t_m_x=[]
    vwc_m_x=[]
    abi_m_x=[]  
 
    t_w_x=[]
    vwc_w_x=[]
    abi_w_x=[]
    
    #x arrays
    chl_d_x=[]
    dia_d_x=[]
    ht_d_x=[]
    sla_d_x=[]
    na_d_x=[]
    bi_d_x=[]

    chl_m_x=[]
    dia_m_x=[]
    ht_m_x=[]
    sla_m_x=[]
    na_m_x=[]
    bi_m_x=[]

    chl_w_x=[]
    dia_w_x=[]
    ht_w_x=[]
    sla_w_x=[]
    na_w_x=[]
    bi_w_x=[]

 

    for xx in range(len(temp_per)):

        #---------------Photosynthesis + Stomatal Conductance Model---------------#

        
        ##---Constant Parameter Arrays for Model---##
     
        #----Params Used in Model Currently----#
          
        tk_25=298.16; #absolute temperature at 25 C
        ekc=80500.0 #Activation energy for K of CO2 (J mol-1) VARIABLE
        eko=14500.0 #Activation energy for K of O2 (J mol-1) VARIABLE
        etau=-29000.0  #Activation energy for tau (???) (J mol-1) VARIABLE
        ev=55000.0 #Activation energy for carboxylation (J mol-1) VARIABLE
        ej=55000.0 #Activation energy for electron transport (J mol-1) VARIABLE
        toptv=303.0 #Optimum temperature for maximum carboxylation (K) VARIABLE
        toptj=303.0 #Optimum temperature for maximum electron transport (K) VARIABLE
        ra=np.zeros(shape=1)+20.7 #specific rubisco activity (umol CO2/g Rub s) VARIABLE
        flnr=np.zeros(shape=1)+0.1 #fraction of leaf nitrogen in rubisco (g N Rub/g N leaf) VARIABLE
        frnr=np.zeros(shape=1)+6.25 #weight fraction of nitrogen in rubisco molecule (g Rub/g N Rub) VARIABLE
        rh=np.zeros(shape=1)+0.5 #relative humidity (kPa/kPa) VARIABLE
        ca=np.zeros(shape=1)+405 #ambient carbon dioxide (umol CO2/mol air)
        ko25=np.zeros(shape=1)+30000 #Michaelis-Menten kinetic coefficient for oxygen at 25 C(Pa) VARIABLE
        kc25=np.zeros(shape=1)+30 #Michaelis-Menten kinetic coefficient for carbon dioxide at 25 C (Pa) VARIABLE
        o=np.zeros(shape=1)+210000 #concentration of ambient oxygen (umol/mol)
        g0=np.zeros(shape=1)+0.002 #Ball-Berry stomatal conductance intercept parameter (mol H2O/m2s) VARIABLE
        a=np.zeros(shape=1)+1.6 #Conversion Coefficient between stomatal conductance to water and carbon dioxide (unitless)
        ij=np.zeros(shape=1)+1.0 #leaf angle index--downregulates jmax
        m=np.zeros(shape=1)+9.0 #ball-berry parameter (unitless) VARIABLE
        b=1.37 #Conversion Coefficient between boundary layer conductance to water and carbon dioxide 
        u=5.0 #windspeed (m/s) VARIABLE
        qeff=0.32 #leaf quantum yield, electrons VARIABLE
        PAR=2000 #photosynthetic active radiation (umol/m2s) VARIABLE
        jm=2.68 #slope coefficient  VARIABLE
        
        vwc_min=0.08 #minimum soil water content for photosynthesis to occur (permanent wilting point) (cm3/cm3)  VARIABLE
        vwc_max=0.68 #maximum soil water content where increases in soil water do not affect photosynthesis (field capacity?) (cm3/cm3) VARIABLE
        q=0.2 #parameter for soil water affect on photosynthesis (unitless) VARIABLE

        #----abiotic parameters that vary for each meadow type----#
        t=[(17)*temp_per[xx],(15)*temp_per[xx],(13)*temp_per[xx]] #surface temperature for dry, moist, and wet meadows VARIABLE
        vwc=[0.16*vwc_per[xx],0.28*vwc_per[xx],0.54*vwc_per[xx]] #volumetric water content for dry, moist, and wet meadows (June through August) VARIABLE
        
        
        #----biotic parameters that vary for each meadow type----#
        na_d=2.45
        dia_d=1.4
        chl_d=396
        ht_d=9.2

        na_m=5.0
        dia_m=2.3
        chl_m=465
        ht_m=19.2
        
        na_w=6.25
        dia_w=2.6
        chl_w=476
        ht_w=20.0
                
        
        #set variable parameters constant if I specify this above
        if 'na' in var_params[ii]:
            na=[np.linspace(na_d*1.0,na_d*0.6,10),np.linspace(na_m*1.0,na_m*0.6,10),np.linspace(na_w*1.0,na_w*0.6,10)]
        else:
            na=[np.linspace(na_d,na_d,10),np.linspace(na_m,na_m,10),np.linspace(na_w,na_w,10)]
        
        if 'dia' in var_params[ii]:
            dia=[np.linspace(dia_d*1.0,dia_d*0.6,10),np.linspace(dia_m*1.0,dia_m*0.6,10),np.linspace(dia_w*1.0,dia_w*0.6,10)]
        else:
            dia=[np.linspace(dia_d,dia_d,10),np.linspace(dia_m,dia_m,10),np.linspace(dia_w,dia_w,10)]
        
        if 'chl' in var_params[ii]:
            chl=[np.linspace(chl_d*1.0,chl_d*0.6,10),np.linspace(chl_m*1.0,chl_m*0.6,10),np.linspace(chl_w*1.0,chl_w*0.6,10)]
        else:
            chl=[np.linspace(chl_d,chl_d,10),np.linspace(chl_m,chl_m,10),np.linspace(chl_w,chl_w,10)]
       
        if 'ht' in var_params[ii]:
            ht=[np.linspace(ht_d*1.0,ht_d*0.6,10),np.linspace(ht_m*1.0,ht_m*0.6,10),np.linspace(ht_w*1.0,ht_w*0.6,10)]
        else:
            ht=[np.linspace(ht_d,ht_d,10),np.linspace(ht_m,ht_m,10),np.linspace(ht_w,ht_w,10)]
           
            


            #---------------Make Array of Values for Each Meadow---------------#  
            
        #biotic index
        biotic_index=np.linspace(0,-40,len(na[0]))
        

        for iii in range(len(na[0])):
        
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

 
   
            t_diff=18-0.4*ht[0][iii]
        
            tl=t_dm+t_diff
    

        

            #---------------Photosynthesis Function---------------#
        
            #alter this line of code for when implementing different photosynthesis functions
            wue, nue, A, E, cs, ci, gsw, gs, gbw, gb, gm, cc,dd =photo(tk_25,ekc,eko,etau,ev,ej,toptv,toptj,na[0][iii], qeff, PAR,tl,ea,chl[0][iii],ij,kc25,ko25,o,ca,rh,m,a,frnr,flnr,ra,jm,g0,b,dia[0][iii],u,q,vwc_min,vwc_max,vwc_dm)     

            
            #test to make sure wue and nue are positive at not 'nan'
            if wue[0]==-999 and nue[0]==-999:
           
                continue
    
 
            wue_d+=[wue[0]]
            nue_d+=[nue[0]]
            gsw_d+=[gsw[0]]
            A_d+=[A[0]]
            E_d+=[E[0]]
            vpd_d+=[dd[0]]
            
            t_d_x+=[t_dm]
            vwc_d_x+=[vwc_dm]
            chl_d_x+=[chl[0][iii]]
            dia_d_x+=[dia[0][iii]]
            ht_d_x+=[ht[0][iii]]
            na_d_x+=[na[0][iii]]
            abi_d_x+=[abiotic_index[xx]]
            bi_d_x+=[b_i]

        
        for iii in range(len(na[0])):

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
   
            t_diff=18-0.4*ht[1][iii]


       
            tl=t_mm+t_diff
            



            #---------------Photosynthesis Function---------------#
        
            #alter this line of code for when implementing different photosynthesis functions
            wue, nue, A, E, cs, ci, gsw, gs, gbw, gb, gm, cc,dd =photo(tk_25,ekc,eko,etau,ev,ej,toptv,toptj,na[1][iii], qeff, PAR,tl,ea,chl[1][iii],ij,kc25,ko25,o,ca,rh,m,a,frnr,flnr,ra,jm,g0,b,dia[1][iii],u,q,vwc_min,vwc_max,vwc_mm)


              
            #test to make sure wue and nue are positive at not 'nan'
            if wue[0]==-999 and nue[0]==-999:
           
                continue
        
            wue_m+=[wue[0]]
            nue_m+=[nue[0]]
            gsw_m+=[gsw[0]]
            A_m+=[A[0]]
            E_m+=[E[0]]
            vpd_m+=[dd[0]]
   
            t_m_x+=[t_mm]
            vwc_m_x+=[vwc_mm]
            chl_m_x+=[chl[1][iii]]
            dia_m_x+=[dia[1][iii]]
            ht_m_x+=[ht[1][iii]]
            na_m_x+=[na[1][iii]]
            abi_m_x+=[abiotic_index[xx]]
            bi_m_x+=[b_i]



        for iii in range(len(na[0])):

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
   
            t_diff=18-0.4*ht[2][iii]
            
            tl=t_wm+t_diff
 

        

            #---------------Photosynthesis Function---------------#
        
            #alter this line of code for when implementing different photosynthesis functions
            wue, nue, A, E, cs, ci, gsw, gs, gbw, gb, gm, cc,dd =photo(tk_25,ekc,eko,etau,ev,ej,toptv,toptj,na[2][iii], qeff, PAR,tl,ea,chl[2][iii],ij,kc25,ko25,o,ca,rh,m,a,frnr,flnr,ra,jm,g0,b,dia[2][iii],u,q,vwc_min,vwc_max,vwc_wm)

            
            #test to make sure wue and nue are positive at not 'nan'
            if wue[0]==-999 and nue[0]==-999:
           
                continue                
        
            wue_w+=[wue[0]]
            nue_w+=[nue[0]]
            gsw_w+=[gsw[0]]
            A_w+=[A[0]]
            E_w+=[E[0]]
            vpd_w+=[dd[0]]
            
            t_w_x+=[t_wm]
            vwc_w_x+=[vwc_wm]
            chl_w_x+=[chl[2][iii]]
            dia_w_x+=[dia[2][iii]]
            ht_w_x+=[ht[2][iii]]
            na_w_x+=[na[2][iii]]
            abi_w_x+=[abiotic_index[xx]]
            bi_w_x+=[b_i]
    
#---------------Regression Plot WUE vs. NUE---------------#   
   

    axA.plot(abi_d_x,nue_d,color='red',linewidth=10,alpha=0.5,label="$\Delta$ NUE: %0.1f" %(max(nue_d)-min(nue_d)))
    axA2.plot(abi_d_x,wue_d,color='blue',linewidth=10,alpha=0.5,label="$\Delta$ WUE: %0.1f" %(max(wue_d)-min(wue_d)))
#    axA.invert_xaxis()
    axA.legend(loc="upper left",fontsize=25)
    axA2.legend(loc="upper right",fontsize=25)

    axB.plot(abi_m_x,nue_m,color='red',linewidth=10,alpha=0.5,label="$\Delta$ NUE: %0.1f" %(max(nue_m)-min(nue_m)))
    axB2.plot(abi_m_x,wue_m,color='blue',linewidth=10,alpha=0.5,label="$\Delta$ WUE: %0.1f" %(max(wue_m)-min(wue_m)))
 #   axB.invert_xaxis()
    axB.legend(loc="upper left",fontsize=25)
    axB2.legend(loc="upper right",fontsize=25)
    
    axC.plot(abi_w_x,nue_w,color='red',linewidth=10,alpha=0.5,label="$\Delta$ NUE: %0.1f" %(max(nue_w)-min(nue_w)))
    axC2.plot(abi_w_x,wue_w,color='blue',linewidth=10,alpha=0.5,label="$\Delta$ WUE: %0.1f" %(max(wue_w)-min(wue_w)))
#    axC.invert_xaxis()  
    axC.legend(loc="upper left",fontsize=25)
    axC2.legend(loc="upper right",fontsize=25)     
    
#3D Plots

    axD.plot_trisurf(abi_d_x,bi_d_x,nue_d,linewidth=0.2, antialiased=False,color="red",alpha=0.5)
#    axD.set_zticks(np.arange(min(nue_d),max(nue_d)+0.2,0.2))
    axE.plot_trisurf(abi_d_x,bi_d_x,wue_d,linewidth=0.2, antialiased=False,color="blue",alpha=0.5)

    axF.plot_trisurf(abi_m_x,bi_m_x,nue_m,linewidth=0.2, antialiased=False,color="red",alpha=0.5)
#    axF.set_zticks(np.arange(min(nue_m),max(nue_m)+0.2,0.2))
    axG.plot_trisurf(abi_m_x,bi_m_x,wue_m,linewidth=0.2, antialiased=False,color="blue",alpha=0.5)

    axH.plot_trisurf(abi_w_x,bi_w_x,nue_w,linewidth=0.2, antialiased=False,color="red",alpha=0.5)
#    axH.set_zticks(np.arange(min(nue_w),max(nue_w)+0.2,0.2))
    axI.plot_trisurf(abi_w_x,bi_w_x,wue_w,linewidth=0.2, antialiased=False,color="blue",alpha=0.5)


    axA.tick_params(axis='y', labelsize=24)
    axA.yaxis.label.set_color('red')
    axA.tick_params(axis='y', colors='red')
    axA.tick_params(axis='x', pad=15)


    axB.tick_params(axis='y', labelsize=24)
    axB.yaxis.label.set_color('red')
    axB.tick_params(axis='y', colors='red')    
    axB.tick_params(axis='x', pad=15)

    axC.tick_params(axis='y', labelsize=24)
    axC.yaxis.label.set_color('red')
    axC.tick_params(axis='y', colors='red')    
    axC.tick_params(axis='x', pad=15)

    axA2.tick_params(axis='y', labelsize=24)
    axA2.yaxis.label.set_color('blue')
    axA2.tick_params(axis='y', colors='blue')
    axA2.tick_params(axis='x', pad=15)

    axB2.tick_params(axis='y', labelsize=24)
    axB2.yaxis.label.set_color('blue')    
    axB2.tick_params(axis='y', colors='blue')
    axB2.tick_params(axis='x', pad=15)
    
    axC2.tick_params(axis='y', labelsize=24)
    axC2.yaxis.label.set_color('blue')   
    axC2.tick_params(axis='y', colors='blue')
    axC2.tick_params(axis='x', pad=15)
    
    axA.tick_params(axis='x', labelsize=24)
    axB.tick_params(axis='x', labelsize=24)
    axC.tick_params(axis='x', labelsize=24)
             
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
#    fig1.tight_layout()
#    fig2.tight_layout()
#    fig3.tight_layout()    
#    fig1.savefig('NUE_vs_WUE_dm_abiotic.png')
#    fig2.savefig('NUE_vs_WUE_mm_abiotic.png')
#    fig3.savefig('NUE_vs_WUE_wm_abiotic.png')
  
#3D Figures
    fig4.savefig('NUE_dm_a_biotic')
    fig5.savefig('WUE_dm_a_biotic')
    fig6.savefig('NUE_mm_a_biotic')
    fig7.savefig('WUE_mm_a_biotic')
    fig8.savefig('NUE_wm_a_biotic')
    fig9.savefig('WUE_wm_a_biotic')    