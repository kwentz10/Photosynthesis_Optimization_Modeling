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

#The line of code below is for if I want to input all combinations of changed parameters into my model:
from leaf_parameter_inputs import leaf_params

#The line of code below imports the temperature functions used in this model
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

    #fb1=plt.figure(1,figsize=(12,2)) 
    #axA = fb1.add_subplot(121)
    #axB = fb1.add_subplot(122)

    ##---Figure Without Subplots Blueprint---##

    #put in correct ax value (e.g. axA, axB)
    fig,axA = plt.subplots(figsize=(10,5))


    ##---Define Plot Parameters Based on Graph Interests---##

    axA.set_xlabel('NUE (umol CO2/g N s)',fontsize=23, fontname='Times New Roman')
    axA.set_ylabel('WUE (umol CO2/mmol H2O)',fontsize=23, fontname='Times New Roman')
    axA.set_xlim([0,20])
    axA.set_ylim([0,10])
    axA.set_title('Growth Response Across Four Plant Trait Assemblages: all params vary', fontname='Times New Roman',fontsize=23,fontweight='bold')
#    axA.set_title('Growth Response: constant %s, %s, %s, %s, %s' % (dict_params[ii][0],dict_params[ii][1],dict_params[ii][2],dict_params[ii][3],dict_params[ii][4]), fontname='Times New Roman',fontsize=23,fontweight='bold')

    ##---Line Type for Each Plant---##
    n=16 #number of variable parameter combinations for each meadow type

    color=['k']*n+['r']*n+['y']*n+['g']*n+['b']*n
    
    marker=['d']*n+['o']*n+['v']*n+['*']*n+['^']*n

    ##---Initialize Arrays for Each Meadow---##
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
#        
#        print dict_params[ii]
#        print 's:',s
#        print 'nm:',nm
#        print 'm:',m
#        print 'chl:',chl
#        print 'tl:',tl
#        print 'vwc:',vwc
#        print 'ij:',ij
##---Calculated Parameter Arrays for Model(Constant+Variable Plant Trait(s))---##
        es_str=pa_con_atmfrac(611*np.exp(17.27*(tl-273.15)/((tl-273.15)+273.3))) #calculate saturation vapor pressure of surface (Pa)
        d=es_str-ea #calculate vapor pressure deficit (umol H2O/mol air)
    
        l=1/s #leaf mass per unit area (g C/m2 C)
        na=nm*l #leaf nitrogen (g N/ m2 C)
    
    #below is commented out because I am no longer using a variable lambda parameter
#    m=ca/(rh*d*lamb) ##Ball-Berry stomatal conductance slope parameter (unitless)

        rub=(chl*crc)/1000 # find ribulose bisphosphate content (umol RuBP/m2)
        
        if all(rub<rub_max):
            j_m=j_m_max*(rub/rub_max)*ij #find j_m slope based on ribulose bisphosphate content & leaf area/angle index
        else:
            j_m=j_m_max*ij
            
        vopt=frnr*flnr*ra*na #optimal carboxylation rate, limited by CO2 (umol CO2/m2s)
        jopt=vopt*j_m+j_b #optimal carboxylation rate, limited by RuBP (umol CO2/m2s)

##---Temperature Effects on Parameters---##
    
        #parameters
        tk_25=298.16; #absolute temperature at 25 C
        ekc=80500.0 #Activation energy for K of CO2 (J mol-1)
        eko=14500.0 #Activation energy for K of O2 (J mol-1)
        etau=-29000.0  #Activation energy for tau (???) (J mol-1)
        ev=55000.0 #Activation energy for carboxylation (J mol-1)
        ej=55000.0 #Activation energy for electron transport (J mol-1)
        toptv=298.0 #Optimum temperature for maximum carboxylation (K)
        toptj=298.0 #Optimum temperature for maximum electron transport (K)
    
        #calculated parameters due to temperature
        kc=arr_temp(kc25,ekc,tk_25,tl) #Michaelis-Menten kinetic coefficient for carbon dioxide at leaf temperature (umol/mol)
        ko=arr_temp(ko25,eko,tk_25,tl) #Michaelis-Menten kinetic coefficient for oxygen at leaf temperature (umol/mol) 
        tau=arr_temp(tau25,etau,tk_25,tl) #specifity coefficient of tau at leaf temperature (unitless) 
        gamma=o/(2*tau) #carbon dioxide compensation point (umol/mol)
        vmax1=bol_temp(vopt,ev,toptv,tl) #carboxylation rate at leaf temperature, limited by CO2 (umol CO2/m2s)
        jmax1=bol_temp(jopt,ej,toptj,tl) #carboxylation rate at leaf temperature, limited by RuBP (umol CO2/m2s)

##---Soil Moisture Effect on Parameters---##
    
        if all(vwc>=vwc_max):
            Wfac=1
        elif all(vwc<vwc_max):
            Wfac=((vwc-vwc_min)/(vwc_max-vwc_min))**q
    
        vmax=Wfac*vmax1
        jmax=Wfac*jmax1
    
##---Define a1 and a2 depending on whether plant is rubisco limited or light limited---##

        #rubisco limited
        a1_r=vmax
        a2_r=kc*(1+(o/ko))
        #light limited
        a1_l=jmax/4
        a2_l=2*gamma

        
##---(1)Photosynthesis and Stomatal Conductance Models (b is not taken into account)---##

        if any(b==0.0):
    
            #In order to generate this model I combined the following equations:
            #A=gsc*(ca-ci)
            #gsc=gsw/a
            #gsw=mArh/ca
            #solve for ci
            #plug into A=a1(ci-gamma)/ci+a2
            #Rubisco Limiting: a1=vcmax; a2=kc(1+o/ko)
            #Light Limiting: a1=2.2*vcmax/4; a2=2*gamma

            #Solve for Assimilation
        
            ci=ca-((a*ca)/(m*rh)) #internal carbon dioxide (umol/mol)
        
            ##---Rubisco-Limited Assimilation---##
            A_r=(a1_r*(ci-gamma))/(ci+a2_r) #rubisco limited assimilation rate (umol CO2/m2s)
        
            ##---Light-Limited Assimilation---##
            A_l=(a1_l*(ci-gamma))/(ci+a2_l) #light limited assimilation rate (umol CO2/m2s)
        
            ##---Determine Rate-Limiting Assimilation---##
            A=[]
            for xx in range(len(A_r)):
                if A_r[xx]<A_l[xx]:
                    A+=[A_r[xx]] #rubisco limited
                elif A_l[xx]<A_r[xx]:
                    A+=[A_l[xx]] #light limited
                else: 
                    A+=[A_l[xx]] #both light and rubisco limited
                
            ##---Solve for Stomatal Conductance to Water---##
            gsw=m*A*rh/ca #stomatal conductance to water (mol air/m2s)
        
            ##---Solve for Evapotranspiration---##
            E=gsw*d #(umol H2O/m2s)
    

##---(2)Photosynthesis and Stomatal Conductance Models (with b)---##

        elif any(b>0.0):
    
            #In order to generate this model I combined the following equations:
            #A=gsc*(ca-ci)
            #gsc=gsw/a
            #gsw=mArh/ca+b
            #solve for ci
            #plug into A=a1(ci-gamma)/ci+a2
            #Rubisco Limiting: a1=vcmax; a2=kc(1+o/ko)
            #Light Limiting: a1=2.2*vcmax/4; a2=2*gamma

            #Solve for Assimilation Using Quadratic Equation
        
            ##---Rubisco-Limited Assimilation---##
            aa_r=m*rh*ca-a*ca+m*rh*a2_r
            bb_r=b*(ca**2)+b*ca*a2_r-a1_r*m*rh*ca+a*ca*a1_r+a1_r*m*rh*gamma
            cc_r=a1_r*b*(ca**2)+gamma*b*ca*a1_r

            A1_r=(-bb_r+np.sqrt(bb_r**2-4*aa_r*cc_r))/(2*aa_r)
            A2_r=(-bb_r-np.sqrt(bb_r**2-4*aa_r*cc_r))/(2*aa_r)
            
            #Choose Highest Values for Assimilation and Conductance
            A_r=[]
            for j in range(len(A1_r)):
                if A1_r[j]>A2_r[j]:
                    A_r+=[A1_r[j]]
                elif A2_r[j]>A1_r[j]:
                    A_r+=[A2_r[j]]
                else:
                    A_r+=[A1_r[j]]
        
            ##---Light-Limited Assimilation---##
            aa_l=m*rh*ca-a*ca+m*rh*a2_l
            bb_l=b*(ca**2)+b*ca*a2_l-a1_l*m*rh*ca+a*ca*a1_l+a1_l*m*rh*gamma
            cc_l=a1_l*b*(ca**2)+gamma*b*ca*a1_l

            A1_l=(-bb_l+np.sqrt(bb_l**2-4*aa_l*cc_l))/(2*aa_l)
            A2_l=(-bb_l-np.sqrt(bb_l**2-4*aa_l*cc_l))/(2*aa_l)
            
            #Choose Highest Values for Assimilation and Conductance
            A_l=[]
            for j in range(len(A1_l)):
                if A1_l[j]>A2_l[j]:
                    A_l+=[A1_l[j]]
                elif A2_l[j]>A1_l[j]:
                    A_l+=[A2_l[j]]
                else:
                    A_l+=[A1_l[j]]

            ##---Determine Rate-Limiting Assimilation---##
            A=[]
            for xx in range(len(A_r)):
                if A_r[xx]<A_l[xx]:
                    A+=[A_r[xx]] #rubisco limited
                elif A_l[xx]<A_r[xx]:
                    A+=[A_l[xx]] #light limited
                else: 
                    A+=[A_l[xx]] #both light and rubisco limited         
        
            ##---Solve for Stomatal Conductance to Water---##
            gsw=m*A*rh/ca #stomatal conductance to water (mol H2O/m2s) #make array from list
        
            ##---Solve for Evapotranspiration---##
            E=gsw*d #(umol H2O/m2s)

        
#---------------Test for Nan or Negative Values---------------#       
    
        for xxx in range(len(A)):
            if np.isnan(A[xxx]):
                print "A array contains nan values"
                break
            if A[xxx]<0.0:
                print "A array contains negative values"
                break
            if np.isnan(gsw[xxx]):
                print "gsw array contains nan values"
                break
            if gsw[xxx]<0.0:
                print "gsw array contains negative values"
                break
            if np.isnan(E[xxx]):
                print "E array contains nan values"
                break
            if E[xxx]<0.0:
                print "E array contains negative values"
                break

        
#---------------WUE vs. NUE---------------#    
    
    
        wue=np.diff(A)/np.diff(E)*1000.0 #multiply by 1000 to get from umol CO2/umol H20 to umol CO2/mmol H20
        nue=np.diff(A)/np.diff(na)

    
#---------------Test for Low NUE Values---------------#  

#    if any(nue<15):
#        break
    

    
#---------------Make Array of Values for Each Meadow---------------#  

        if i+1<17:
            wue_f+=[wue]
            nue_f+=[nue]
        elif i+1>=17 and i+1<33:
            wue_d+=[wue]
            nue_d+=[nue]
        elif i+1>=33 and i+1<49:
            wue_m+=[wue]
            nue_m+=[nue]
        elif i+1>=49 and i+1<65:
            wue_w+=[wue]
            nue_w+=[nue]
        elif i+1>=65 and i+1<81:
            wue_s+=[wue]
            nue_s+=[nue]

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
   
    axA.plot([max(nue_f),max(nue_f),min(nue_f),min(nue_f),max(nue_f)],[max(wue_f),min(wue_f),min(wue_f),max(wue_f),max(wue_f)],color='red', linestyle='-', linewidth=8.5, label='fellfield') 
#    axA.scatter([max(nue_f),max(nue_f),min(nue_f),min(nue_f),max(nue_f)],[max(wue_f),min(wue_f),min(wue_f),max(wue_f),max(wue_f)],color='k', marker='^', label='fellfield') 
    axA.fill_between([np.min(nue_f),np.max(nue_f)],np.min(wue_f),np.max(wue_f),color='red',alpha=0.5) 
    
    
    axA.plot([max(nue_d),max(nue_d),min(nue_d),min(nue_d),max(nue_d)],[max(wue_d),min(wue_d),min(wue_d),max(wue_d),max(wue_d)],color='brown', linestyle='-',linewidth=6.5, label='dry meadow')
#    axA.scatter([max(nue_d),max(nue_d),min(nue_d),min(nue_d),max(nue_d)],[max(wue_d),min(wue_d),min(wue_d),max(wue_d),max(wue_d)],color='r', marker='d',label='dry meadow')     
    axA.fill_between([np.min(nue_d),np.max(nue_d)],np.min(wue_d),np.max(wue_d),color='brown',alpha=0.5) 
    
    
    axA.plot([max(nue_m),max(nue_m),min(nue_m),min(nue_m),max(nue_m)],[max(wue_m),min(wue_m),min(wue_m),max(wue_m),max(wue_m)],color='orange', linestyle='-',linewidth=4.5, label='moist meadow')
#    axA.scatter([max(nue_m),max(nue_m),min(nue_m),min(nue_m),max(nue_m)],[max(wue_m),min(wue_m),min(wue_m),max(wue_m),max(wue_m)],color='y', marker='o',s=40,label='moist meadow') 
    axA.fill_between([np.min(nue_m),np.max(nue_m)],np.min(wue_m),np.max(wue_m),color='orange',alpha=0.5) 
    
    
    axA.plot([max(nue_w),max(nue_w),min(nue_w),min(nue_w),max(nue_w)],[max(wue_w),min(wue_w),min(wue_w),max(wue_w),max(wue_w)],color='g', linestyle='-',linewidth=2.0,label='wet meadow')
#    axA.scatter([max(nue_w),max(nue_w),min(nue_w),min(nue_w),max(nue_w)],[max(wue_w),min(wue_w),min(wue_w),max(wue_w),max(wue_w)],facecolor='none',edgecolor='g', marker='o',s=120,label='wet meadow')    
    axA.fill_between([np.min(nue_w),np.max(nue_w)],np.min(wue_w),np.max(wue_w),color='green',alpha=0.5) 
    
    
    axA.plot([max(nue_s),max(nue_s),min(nue_s),min(nue_s),max(nue_s)],[max(wue_s),min(wue_s),min(wue_s),max(wue_s),max(wue_s)],color='c', linestyle='-',linewidth=1.0, label='snowbed')
#    axA.scatter([max(nue_s),max(nue_s),min(nue_s),min(nue_s),max(nue_s)],[max(wue_s),min(wue_s),min(wue_s),max(wue_s),max(wue_s)],color='b', marker='*', label='snowbed')         
    axA.fill_between([np.min(nue_s),np.max(nue_s)],np.min(wue_s),np.max(wue_s),color='cyan',alpha=0.5) 
    
   
    
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
    axA.legend(bbox_to_anchor=(1, 1), loc='left', prop={'size':11})

    ##---Save Figure--##
    plt.savefig('NUE_vs_WUE_variable_test.png') 
