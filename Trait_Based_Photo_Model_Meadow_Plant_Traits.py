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
The end product is graphs of 
the plant trait vs. NUE and WUE.


Update: I am running a sensitivity analysis of 
different plant traits in model

"""

#---------------Import Modules---------------#

import numpy as np
from matplotlib import pyplot as plt

#The line of code below is for if I want to input all combinations of changed parameters into my model:
from leaf_parameter_inputs import leaf_params

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
#axA.set_xlim([19,44])
axA.set_title('Growth Response Across Four Plant Trait Assemblages', fontname='Times New Roman',fontsize=23,fontweight='bold')

##---Line Type for Each Plant---##

color=['g','g',
'y','y',
'b','b',
'r','r']

#marker=[]
#linestyle=[]

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

##---Temperature---##
t=25 #degrees C

##---Parameter Arrays for Model (Constant)---##

#I have commented out parameters that I am assuming are variable (for the time being)

#s=np.zeros(shape=100)+0.02 #specific leaf area (m2 C/g C)
#ra=np.zeros(shape=100)+20.7 #specific rubisco activity (umol CO2/g Rub s)
#nm=np.zeros(shape=100)+0.03 #leaf nitrogen (g N/ g C)
#flnr=np.zeros(shape=100)+0.6 #fraction of leaf nitrogen in rubisco (g N Rub/g N leaf)
frnr=np.zeros(shape=3)+6.25 #weight fraction of nitrogen in rubisco molecule (g Rub/g N Rub) 
e_str=611*np.exp(17.27*t/(t+273.3)) #saturation vapor pressure (Pa)
rh=np.zeros(shape=3)+0.6 #relative humidity (%)
d=(e_str*(1-rh))*9.9 #vapor pressure deficit (umol/mol)-->multiply by 9.9 to get from Pa to umol/mol
ca=np.zeros(shape=3)+410 #ambient carbon dioxide (umol/mol)
gamma=np.zeros(shape=3)+29.6 #carbon dioxide compensation point (umol/mol)
ko=np.zeros(shape=3)+296078 #325685 farquhar #kinetic coefficient for oxygen (umol/mol)
kc=np.zeros(shape=3)+ 296 #454 farquhar #kinetic coefficient for carbon dioxide (umol/mol)
o=np.zeros(shape=3)+210000 #concentration of ambient oxygen (umol/mol)
#m=np.zeros(shape=9)+9 #Ball-Berry stomatal conductance slope parameter
b=np.zeros(shape=3)+0.01 #Ball-Berry stomatal conductance intercept parameter
a=np.zeros(shape=3)+1.6 #Conversion Coefficient between stomatal conductance to water and carbon dioxide 
chl=np.zeros(shape=3)+450 #Chlorophyll Content of leaves (umol/m2)


##---Parameter Arrays for Model (Variable)---##

for i in range(len(leaf_params)):
    for key,val in leaf_params[i].items():
        exec(key + '=val')


##---Calculated Parameter Arrays for Model(Constant+Variable Plant Trait(s))---##

    l=1/s #leaf mass per unit area (g C/m2 C)
    na=nm*l #leaf nitrogen (g N/ m2 C)
    rub=(chl*crc)/1000 # find ribulose bisphosphate content (umol RuBP/m2)
    j_m=j_m_max*(rub/rub_max) #find j_m slope based on ribulose bisphosphate content


##---Define a1 and a2 depending on whether plant is rubisco limited or light limited---##

    #rubisco limited
    a1_r=frnr*flnr*ra*na
    a2_r=kc*(1+(o/ko))
    #light limited
    a1_l=((frnr*flnr*ra*na)*j_m+j_b)/4
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
        gsw=m*A*rh/ca #stomatal conductance to water (mol H2O/m2s)
        
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

        
#---------------Test for Nan Values---------------#       

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
   
    axA.scatter(nue,wue,color=color[i]) 
    
    
    ##---Separate Plots Into Different Figures: Legend---##
    #ax1.legend(loc=4)


#---------------Make Plot Interactive---------------# 
    
    plt.pause(0.0001)
    plt.ion()
    #end of sensitivity analysis iterations
    
#---------------Finalize Figure---------------#    

#axA refers to first figure in subplot; axB refers to second figure in subplot
#if only one axis is run then the figure is just one plot

##---Legend---##
axA.legend(loc=2,prop={'size':11})

##---Save Figure--##
plt.savefig('NUE_vs_WUE.png') 
