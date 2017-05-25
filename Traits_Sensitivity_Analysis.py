#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 25 11:15:02 2017

@author: Katherine
"""

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


#---------------Initialize Plots---------------#

   
#--figure 1--#
fig1,axA = plt.subplots(figsize=(15,15))
axA.set_xlabel('Trait Variability (%)',fontsize=30, fontname='Times New Roman')
axA.set_ylabel('Assimilation (umol CO2/m2 s)',fontsize=30, fontname='Times New Roman')
axA.set_title('Dry Meadow: Assimilation Sensitivity to Leaf Trait Variability', fontname='Times New Roman',fontsize=30,fontweight='bold')



#--figure 2--#

fig2,axB = plt.subplots(figsize=(15,15))
axB.set_xlabel('Trait Variability (%)',fontsize=30, fontname='Times New Roman')
axB.set_ylabel('Assimilation (umol CO2/m2 s)',fontsize=30, fontname='Times New Roman')
axB.set_title('Moist Meadow: Assimilation Sensitivity to Leaf Trait Variability', fontname='Times New Roman',fontsize=30,fontweight='bold')


#--figure 3--#

fig3,axC = plt.subplots(figsize=(15,15))
axC.set_xlabel('Trait Variability (%)',fontsize=30, fontname='Times New Roman')
axC.set_ylabel('Assimilation (umol CO2/m2 s)',fontsize=30, fontname='Times New Roman')
axC.set_title('Wet Meadow: Assimilation Sensitivity to Leaf Trait Variability', fontname='Times New Roman',fontsize=30,fontweight='bold')


#set x axis for graphs
percent=np.linspace(-20,20,10)


#initialize assimilation arrays
A_tot_dm=[]
A_tot_mm=[]
A_tot_wm=[]
  

#---------------Loop Through Variables---------------#

var_params=['na','dia','chl','ht']

for ii in range(len(var_params)):
    


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
    t=[17,15,13] #surface temperature for dry, moist, and wet meadows VARIABLE
    vwc=[0.12,0.15,0.27] #volumetric water content for dry, moist, and wet meadows (June through August) VARIABLE
    
    
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
        na=[np.linspace(na_d*0.8,na_d*1.2,10),np.linspace(na_m*0.8,na_m*1.2,10),np.linspace(na_w*0.8,na_w*1.2,10)]
    else:
        na=[np.linspace(na_d,na_d,10),np.linspace(na_m,na_m,10),np.linspace(na_w,na_w,10)]
    
    if 'dia' in var_params[ii]:
        dia=[np.linspace(dia_d*0.8,dia_d*1.2,10),np.linspace(dia_m*0.8,dia_m*1.2,10),np.linspace(dia_w*0.8,dia_w*1.2,10)]
    else:
        dia=[np.linspace(dia_d,dia_d,10),np.linspace(dia_m,dia_m,10),np.linspace(dia_w,dia_w,10)]
    
    if 'chl' in var_params[ii]:
        chl=[np.linspace(chl_d*0.8,chl_d*1.2,10),np.linspace(chl_m*0.8,chl_m*1.2,10),np.linspace(chl_w*0.8,chl_w*1.2,10)]
    else:
        chl=[np.linspace(chl_d,chl_d,10),np.linspace(chl_m,chl_m,10),np.linspace(chl_w,chl_w,10)]
   
    if 'ht' in var_params[ii]:
        ht=[np.linspace(ht_d*0.8,ht_d*1.2,10),np.linspace(ht_m*0.8,ht_m*1.2,10),np.linspace(ht_w*0.8,ht_w*1.2,10)]
    else:
        ht=[np.linspace(ht_d,ht_d,10),np.linspace(ht_m,ht_m,10),np.linspace(ht_w,ht_w,10)]
       
        
    #------initialize sensitivity array----#
    
    A_dm=[]
    A_mm=[]
    A_wm=[]

    #---------------Make Array of Values for Each Meadow---------------#  
        
    

    for iii in range(10):

        
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


        A_dm+=[A[0]]


    
    for iii in range(10):


        
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
    
        A_mm+=[A[0]]



    for iii in range(10):



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
    
        A_wm+=[A[0]]
 
    A_tot_dm+=[A_dm]   
    A_tot_mm+=[A_mm]
    A_tot_wm+=[A_wm]         
#---------------Plot---------------#   
   
axA.plot(percent,A_tot_dm[0],color='red',linewidth=10,alpha=0.5,label="Leaf Nitrogen Content")
axA.plot(percent,A_tot_dm[1],color='blue',linewidth=10,alpha=0.5,label="Leaf Diameter")
axA.plot(percent,A_tot_dm[2],color='green',linewidth=10,alpha=0.5,label="Leaf Chlorophyll Content")
axA.plot(percent,A_tot_dm[3],color='black',linewidth=10,alpha=0.5,label="Leaf Height")    
axA.legend(loc="upper left",fontsize=20)

axB.plot(percent,A_tot_mm[0],color='red',linewidth=10,alpha=0.5,label="Leaf Nitrogen Content")
axB.plot(percent,A_tot_mm[1],color='blue',linewidth=10,alpha=0.5,label="Leaf Diameter")
axB.plot(percent,A_tot_mm[2],color='green',linewidth=10,alpha=0.5,label="Leaf Chlorophyll Content")
axB.plot(percent,A_tot_mm[3],color='black',linewidth=10,alpha=0.5,label="Leaf Height")    
axB.legend(loc="upper left",fontsize=20)

axC.plot(percent,A_tot_wm[0],color='red',linewidth=10,alpha=0.5,label="Leaf Nitrogen Content")
axC.plot(percent,A_tot_wm[1],color='blue',linewidth=10,alpha=0.5,label="Leaf Diameter")
axC.plot(percent,A_tot_wm[2],color='green',linewidth=10,alpha=0.5,label="Leaf Chlorophyll Content")
axC.plot(percent,A_tot_wm[3],color='black',linewidth=10,alpha=0.5,label="Leaf Height")    
axC.legend(loc="upper left",fontsize=20)


axA.tick_params(axis='y', labelsize=24)
axA.tick_params(axis='x', pad=15,labelsize=24)


axB.tick_params(axis='y', labelsize=24)   
axB.tick_params(axis='x', pad=15,labelsize=24)

axC.tick_params(axis='y', labelsize=24)  
axC.tick_params(axis='x', pad=15,labelsize=24)


         



#---------------Finalize Figure---------------#    


#---Save Figure--##
fig1.tight_layout()
fig2.tight_layout()
fig3.tight_layout()    
fig1.savefig('DM_Assimilation_Sensitivity.png')
fig2.savefig('MM_Assimilation_Sensitivity.png')
fig3.savefig('WM_Assimilation_Sensitivity.png')
  
  