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
"""

import numpy as np
from matplotlib import pyplot as plt

#---------------Photosynthesis + Stomatal Conductance Model---------------#

##---Rubisco or Light Limiting?---##

r_l="true"
l_l="false"

##---Parameter Arrays for Model (Constant)---##

s=np.zeros(shape=100)+0.02 #specific leaf area (m2 C/g C)
r_a=np.zeros(shape=100)+20.7 #specific rubisco activity (umol CO2/g Rub s)
n_m=np.zeros(shape=100)+0.03 #leaf nitrogen (g N/ g C)
f_lnr=np.zeros(shape=100)+0.5 #fraction of leaf nitrogen in rubisco (g N Rub/g N leaf)
f_rnr=np.zeros(shape=100)+6.25 #weight fraction of nitrogen in rubisco molecule (g Rub/g N Rub) 
rh=np.zeros(shape=100)+100 #relative humidity (%)
ca=np.zeros(shape=100)+410 #ambient carbon dioxide (umol/mol)
gamma=np.zeros(shape=100)+29.6 #carbon dioxide compensation point (umol/mol)
ko=np.zeros(shape=100)+296078 #kinetic coefficient for oxygen (umol/mol)
kc=np.zeros(shape=100)+296 #kinetic coefficient for carbon dioxide (umol/mol)
o=np.zeros(shape=100)+210000 #concentration of ambient oxygen (umol/mol)
m=np.zeros(shape=100)+5 #Ball-Berry stomatal conductance slope parameter
b=np.zeros(shape=100)+0 #Ball-Berry stomatal conductance intercept parameter
a=np.zeros(shape=100)+1.6 #Conversion Coefficient between stomatal conductance to water and carbon dioxide 

##---Parameter Arrays for Model(Variable Plant Trait(s))---##


#CHANGE-----##################################

#s=np.linspace(0.01,0.06,100) #specific leaf area (m2 C/g C)
#n_m=np.linspace(0.056,0.0028,100) #leaf nitrogen gN/gC
#f_lnr=np.linspace(0.1,0.9,100)
r_a=np.linspace(15,25,100)
m=np.linspace(1,10,100)
#b=np.linspace(0,5,100)

#CHANGE-----##################################


##---Calculated Parameter Arrays for Model(Constant+Variable Plant Trait(s))---##

l=1/s #leaf mass per unit area (g C/m2 C)
n_a=n_m*l #leaf nitrogen (g N/ m2 C)


##---Define a1 and a2 depending on whether plant is rubisco limited or light limited---##

if r_l=="true" and l_l=="false":
    a1=f_rnr*f_lnr*r_a*n_a
    a2=kc*(1+(o/ko))
elif l_l=="true" and r_l=="false":
    a1=(f_rnr*f_lnr*r_a*n_a)*2.2/4
    a2=2*gamma
else:
    print "undefined limitation"

##---(1)Photosynthesis and Stomatal Conductance Models (b is not taken into account)---##

if b.any()==0.0:
    
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
    A=(a1*(ci-gamma))/(ci+a2) #assimilation rate (umol CO2/m2s)

    #Solve for Stomatal Conductance to Water
    gsw=m*A*rh/ca #stomatal conductance to water (mol H2O/m2s)
    

##---(2)Photosynthesis and Stomatal Conductance Models (with b)---##

elif b.any()==1.0:
    
    #In order to generate this model I combined the following equations:
    #A=gsc*(ca-ci)
    #gsc=gsw/a
    #gsw=mArh/ca+b
    #solve for ci
    #plug into A=a1(ci-gamma)/ci+a2
    #Rubisco Limiting: a1=vcmax; a2=kc(1+o/ko)
    #Light Limiting: a1=2.2*vcmax/4; a2=2*gamma

    #Solve for Assimilation Using Quadratic Equation
    aa=m*rh*ca-a*ca+m*rh*a2
    bb=b*(ca**2)+b*ca*a2-a1*m*rh*ca+a*ca*a1+a1*m*rh*gamma
    cc=a1*b*(ca**2)+gamma*b*ca*a1

    A1=(-bb+np.sqrt(bb**2-4*aa*cc))/(2*aa)
    A2=(-bb-np.sqrt(bb**2-4*aa*cc))/(2*aa)

    #Solve for Stomatal Conductance to Water
    gsw1=[(m[x]*A1[x]*rh[x]/ca[x])+b[x] for x in range(len(A1))]
    gsw2=[(m[x]*A2[x]*rh[x]/ca[x])+b[x] for x in range(len(A2))]
    
    A=[]
    gsw=[]      
    #Choose Highest Values for Assimilation and Conductance
    for i in range(len(A1)):
        if A1[i]>A2[i]:
            A+=[A1[i]]
            gsw+=[gsw1[i]]
        elif A2[i]>A1[i]:
            A+=[A2[i]]
            gsw+=[gsw2[i]]
        else:
            A+=[A1[i]]
            gsw+=[gsw1[i]]

    A=np.array(A) #make array from list
    gsw=np.array(gsw) #make array from list
#---------------WUE vs. NUE---------------#    
wue=A/gsw
nue=A/n_a


#---------------Plot NUE vs. WUE---------------#      

#I am plotting NUE vs. WUE for each plant trait
#because I am interested in this relationship
      
fb1=plt.figure(1,figsize=(6,6)) #figure blueprint
fig,ax1 = plt.subplots()
ax1.set_xlabel('NUE (umol CO2/g N s)',fontsize=12)
ax1.set_ylabel('WUE (umol CO2/mol H20)',fontsize=12)
ax1.plot(nue,wue) 
ax1.set_title('NUE vs. WUE',fontsize=14)
#plt.savefig('NUE_vs._WUE.png')
     

##---------------Plot Nitrogen vs. Stomatal Conductance---------------#      
#
##I am plotting leaf nitrogen vs. stomatal conductance for each plant trait
##because I am interested in this relationship
#      
#fb2=plt.figure(2,figsize=(6,6)) #figure blueprint
#fig,ax1 = plt.subplots()
#ax1.set_xlabel('Leaf Nitrogen (g N/m2 C)',fontsize=12)
#ax1.set_ylabel('Stomatal Conductance (mol H2O/m2 s)',fontsize=12)
#ax1.plot(n_a,gsw) 
#ax1.set_title('Leaf Nitrogen vs. Stomatal Conductance',fontsize=14)
##plt.savefig('Leaf_N_vs._Cond_Increasing_SLA.png')
#     
#      
##---------------Nutrient Use Efficiency and Water Use Efficiency--------------- #     
#
##I am plotting the change in plant trait vs. two y axes: WUE and NUE
#
#fb3=plt.figure(3,figsize=(6,6)) #figure blueprint
#fig, ax1 = plt.subplots()
#ax2=ax1.twinx() #second y axis
#ax1.set_xlabel('Photosynthesis (umol CO2/m2s)',fontsize=12)
#ax1.set_ylabel('NUE (umol CO2/g N s)',fontsize=12)
#ax2.set_ylabel('WUE (umol CO2/mol H20)',fontsize=12)
#ax1.plot(A,nue,'-r',label='NUE') 
#ax2.plot(A,wue,'-b',label='WUE') 
#ax1.set_title('NUE vs. WUE with Increasing Photosynthesis',fontsize=14)
#ax1.legend(loc=2)
#ax2.legend(loc=1)
##plt.savefig('NUE_vs._WUE_Increasing_Photo_Increasing_SLA.png')
