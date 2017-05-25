#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat May  6 09:17:57 2017

@author: Katherine
"""

from matplotlib import pyplot as plt
import numpy as np

#NUE
#put in correct ax value (e.g. axA, axB)
fig1,axA = plt.subplots(figsize=(30,15))
    
axA.set_ylabel('NUE (g biomass/g N)',fontsize=36, fontname='Times New Roman')
axA.set_ylim([0,110])
axA.set_title('NUE for Alpine Tundra Plant Communities', fontname='Times New Roman',fontsize=36,fontweight='bold')

axA.bar([1.5,2.5,3.5],[72,76,88], yerr=[2.08,4.33,3.81], edgecolor='black', align="center",width=0.2, color='black',alpha=0.5,error_kw={'ecolor':'black', 'lw':2, 'capsize':5, 'capthick':2})
axA.set_xticks([1.5,2.5,3.5])
axA.set_xticklabels(['Dry Meadow','Moist Meadow','Wet Meadow'],fontname='Times New Roman')
axA.tick_params(axis='x', labelsize=28)
axA.tick_params(axis='y', labelsize=18)
axA.text(1.47,78,"a",fontsize=30)
axA.text(2.47,83,"a",fontsize=30)
axA.text(3.47,94,"b",fontsize=30)

fig1.savefig('NUE_Validation.png') 


#WUE Bowman et al. 1995
#put in correct ax value (e.g. axA, axB)
fig2,axA = plt.subplots(figsize=(30,15))
    
axA.set_ylabel('WUE (umol CO2/mmol H2O)',fontsize=36, fontname='Times New Roman')
axA.set_ylim([0,2])
axA.set_title('WUE for Alpine Tundra Plant Communities', fontname='Times New Roman',fontsize=36,fontweight='bold')

axA.bar([1.5,2.5],[1.49,1.57], yerr=[0.06,0.04], edgecolor='black', align="center",width=0.2, color='black',alpha=0.5,error_kw={'ecolor':'black', 'lw':2, 'capsize':5, 'capthick':2})
axA.set_xticks([1.5,2.5])
axA.set_xticklabels(['Dry Meadow','Moist Meadow'],fontname='Times New Roman')
axA.tick_params(axis='x', labelsize=28)
axA.tick_params(axis='y', labelsize=18)
axA.text(1.485,1.65,"a",fontsize=30)
axA.text(2.485,1.7,"a",fontsize=30)
fig2.savefig('WUE_Validation_1.png') 

##WUE my data--decide whether to put this in (might open a can of worms)
##put in correct ax value (e.g. axA, axB)
#fig1,axA = plt.subplots(figsize=(30,15))
#    
#axA.set_ylabel('WUE (umol CO2/mmol H2O)',fontsize=36, fontname='Times New Roman')
#axA.set_ylim([0,100])
#axA.set_title('WUE for Alpine Tundra Plant Communities', fontname='Times New Roman',fontsize=36,fontweight='bold')
#
#axA.bar([1.5,2.5,3.5],[72,76,88], yerr=[2.08,4.33,3.81], edgecolor='black', align="center",width=0.2, color='red',alpha=1.0,error_kw={'ecolor':'black', 'lw':2, 'capsize':5, 'capthick':2})
#axA.set_xticks([1.5,2.5,3.5])
#axA.set_xticklabels(['Dry Meadow','Moist Meadow','Wet Meadow'],fontname='Times New Roman')
#axA.tick_params(axis='x', labelsize=28)
#axA.tick_params(axis='y', labelsize=18)
#
#fig1.savefig('WUE_Validation_2.png') 

#Assimilation
#put in correct ax value (e.g. axA, axB)
fig3,axA = plt.subplots(figsize=(30,15))
    
axA.set_ylabel('Growth Rate (g C/m2 day)',fontsize=36, fontname='Times New Roman')
axA.set_ylim([0,5])
axA.set_title('Growth Rates of Alpine Tundra Plant Communities', fontname='Times New Roman',fontsize=36,fontweight='bold')

axA.bar([1.5,2.5,3.5],[1.022067, 2.351528, 1.995846], yerr=[0.6347155, 1.8942902, 0.8296679], edgecolor='black', align="center",width=0.2, color='black',alpha=0.5,error_kw={'ecolor':'black', 'lw':2, 'capsize':5, 'capthick':2})
axA.set_xticks([1.5,2.5,3.5])
axA.set_xticklabels(['Dry Meadow','Moist Meadow','Wet Meadow'],fontname='Times New Roman')
axA.tick_params(axis='x', labelsize=28)
axA.tick_params(axis='y', labelsize=18)
axA.text(1.47,1.8,"a",fontsize=30)
axA.text(2.47,4.3,"b",fontsize=30)
axA.text(3.47,3.0,"b",fontsize=30)

fig3.savefig('Assimilation_Validation.png') 

#VPD vs. WUE


fig4,axG = plt.subplots(figsize=(15,15))
axG.set_xlabel('VPD (umol H2O/mol air)',fontsize=36, fontname='Times New Roman')
axG.set_ylabel('WUE (umol CO2/mmol H2O)',fontsize=36, fontname='Times New Roman')
axG.set_title('VPD vs. WUE for Alpine Plants', fontname='Times New Roman',fontsize=36,fontweight='bold')
axG.scatter([13734.61196533664,13734.61196533664,13734.61196533664,13734.61196533664,21364.95194607922,21364.95194607922,21364.95194607922,21364.95194607922,28995.2919268218,28995.2919268218,28995.2919268218,28995.2919268218,13734.61196533664,13734.61196533664,13734.61196533664,13734.61196533664,21364.95194607922,21364.95194607922,21364.95194607922,21364.95194607922,28995.2919268218,28995.2919268218,28995.2919268218,28995.2919268218],np.array([11,12,11,3,5,4.5,1,5,4,2,1,.5,9,9,10,11,5,5,5,5,5,5,5,2])*.4091,edgecolors="black",facecolors="black",marker='o',s=50)
axG.plot(np.unique([13734.61196533664,13734.61196533664,13734.61196533664,13734.61196533664,21364.95194607922,21364.95194607922,21364.95194607922,21364.95194607922,28995.2919268218,28995.2919268218,28995.2919268218,28995.2919268218,13734.61196533664,13734.61196533664,13734.61196533664,13734.61196533664,21364.95194607922,21364.95194607922,21364.95194607922,21364.95194607922,28995.2919268218,28995.2919268218,28995.2919268218,28995.2919268218]), np.poly1d(np.polyfit([13734.61196533664,13734.61196533664,13734.61196533664,13734.61196533664,21364.95194607922,21364.95194607922,21364.95194607922,21364.95194607922,28995.2919268218,28995.2919268218,28995.2919268218,28995.2919268218,13734.61196533664,13734.61196533664,13734.61196533664,13734.61196533664,21364.95194607922,21364.95194607922,21364.95194607922,21364.95194607922,28995.2919268218,28995.2919268218,28995.2919268218,28995.2919268218], np.array([11,12,11,3,5,4.5,1,5,4,2,1,.5,9,9,10,11,5,5,5,5,5,5,5,2])*.4091, 1))(np.unique([13734.61196533664,13734.61196533664,13734.61196533664,13734.61196533664,21364.95194607922,21364.95194607922,21364.95194607922,21364.95194607922,28995.2919268218,28995.2919268218,28995.2919268218,28995.2919268218,13734.61196533664,13734.61196533664,13734.61196533664,13734.61196533664,21364.95194607922,21364.95194607922,21364.95194607922,21364.95194607922,28995.2919268218,28995.2919268218,28995.2919268218,28995.2919268218])),color="black")

fig4.savefig("VPD_vs_WUE_Validation")


####-NOT A VALIDATION FIGURE BUT NEEDED FOR DERIVING EQUATION-#####

#Leaf Height vs. Temp Differential (Tleaf-Tair) (r^2=0.2209 deltaT=-0.4366*leaf height +18.3381)

fig5,axG = plt.subplots(figsize=(15,15))
axG.set_xlabel('Leaf Height (cm)',fontsize=36, fontname='Times New Roman')
axG.set_ylabel('Leaf-Air Temperature (C)',fontsize=36, fontname='Times New Roman')
axG.set_title('Leaf Height vs. Leaf-Air Temperature Difference', fontname='Times New Roman',fontsize=36,fontweight='bold')
axG.scatter([25,25,25,9,9,5,5,5,5,3,3,3,3,10,10,10,2],[3.6000,5.8000,12.5000,17.2000,20.8000,13.1000,7.2000,5.8000,11.2000,14.0000,17.7000,15.5000,24.0000,14.5000,6.3000,24.0000,30.0000],edgecolors="black",facecolors="black",marker='o',s=50)
axG.plot(np.unique([25,25,25,9,9,5,5,5,5,3,3,3,3,10,10,10,2]), np.poly1d(np.polyfit([25,25,25,9,9,5,5,5,5,3,3,3,3,10,10,10,2], [3.6000,5.8000,12.5000,17.2000,20.8000,13.1000,7.2000,5.8000,11.2000,14.0000,17.7000,15.5000,24.0000,14.5000,6.3000,24.0000,30.0000], 1))(np.unique([25,25,25,9,9,5,5,5,5,3,3,3,3,10,10,10,2])),color="black")

fig5.savefig("LeafHeight_vs_TempDiff")

####-NOT A VALIDATION FIGURE BUT NEEDED FOR DERIVING EQUATION-#####

#SLA vs. Nitrogen Allocated to Rubisco (r^2=0.5164 flnr=0.004*SLA+0.0703)
fig6,axG = plt.subplots(figsize=(15,15))
axG.set_xlabel('SLA (m2/kg)',fontsize=36, fontname='Times New Roman')
axG.set_ylabel('Fraction of Leaf N in Rubisco (g N in Rubisco/g N in Leaf)',fontsize=36, fontname='Times New Roman')
axG.set_title('SLA vs. Fraction of Leaf N in Rubisco', fontname='Times New Roman',fontsize=36,fontweight='bold')
axG.scatter([13.0319,13.2181,15.2660,18.0585,22.4335,25.8777,28.7633,29.3218,26.6223,29.0426],[0.1157,0.1385,0.1035,0.1834,0.1353,0.1520,0.1680,0.1753,0.1958,0.2247],edgecolors="black",facecolors="black",marker='o',s=50)
axG.plot(np.unique([13.0319,13.2181,15.2660,18.0585,22.4335,25.8777,28.7633,29.3218,26.6223,29.0426]), np.poly1d(np.polyfit([13.0319,13.2181,15.2660,18.0585,22.4335,25.8777,28.7633,29.3218,26.6223,29.0426], [0.1157,0.1385,0.1035,0.1834,0.1353,0.1520,0.1680,0.1753,0.1958,0.2247], 1))(np.unique([13.0319,13.2181,15.2660,18.0585,22.4335,25.8777,28.7633,29.3218,26.6223,29.0426])),color="black")

fig6.savefig("SLA_vs_FLNR")