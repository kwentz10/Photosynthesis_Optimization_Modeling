#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue May 23 16:04:21 2017

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
import matplotlib as mpl

mpl.rc('font',family='Times New Roman')

#Import combinations of variable parameters 
from uncertain_params import monte_carlo


#Import photosynthesis model
from Photosynthesis_Model import photo_bound_meso_eqstom as photo

#Import functions to switch between Pa and umol/mol at sea level
from photo_functions import pa_con_atmfrac, gs_smooth

#for 3d figure and mulivariate linear regression
#from mpl_toolkits.mplot3d import Axes3D
#from sklearn import linear_model


#import timeseries of vwc and temp
from time_dep_params import na_mm_min_inter,na_mm_max_inter,surtemp_dm, surtemp_mm, surtemp_wm, vwc_dm, vwc_mm, vwc_wm, vwc_dm_dy,vwc_wm_dy,surtemp_dm_dy,surtemp_wm_dy,na_dm_min_inter,na_wm_min_inter,na_dm_max_inter, na_wm_max_inter,na_dm_min_inter_dy,na_dm_max_inter_dy,na_wm_min_inter_dy,na_wm_max_inter_dy,gs0_d,gs0_w,gs0_ddy,gs0_wdy,gsf,gs0_m,vwc_mm_dy,surtemp_mm_dy




surtemp_dm_new=list(surtemp_dm)

surtemp_dm_new[181:213]=np.array(surtemp_dm_new[181:213])+4.0


surtemp_wm_new=list(surtemp_wm)
surtemp_wm_new[181:213]=np.array(surtemp_wm_new[181:213])+4.0


#for averaging two lists
def avg(x,y):
    list=[]
    for i in range(len(x)):
        list+=[np.mean(x[i],y[i])]
    return list

#or [np.mean([x,y]) for x, y in zip(list1,list2)]
        

#for doing letters on graphs for multiple plots in same figure
def get_axis_limits(ax, scale1=.95,scale2=0.9):
    return ax.get_xlim()[1]*scale1, ax.get_ylim()[1]*scale2


#---------------Determine if I Want to Keep Any of the Variable Parameters Constant---------------#

const_params=[]
for xxx in it.combinations(['na'],0): #keep ht and t constant for constant vpd
    const_params+=[xxx]

#do this when I do not put any of the variable parameters as constant. instead I 
#vary each parameter one at a time while keeping the other parameters constant.
if const_params==[()]:
    const_params=[[-999999]]   

#---------------Begin Looping Through Photosynthesis Model---------------#

#each loop is for a constant value, or combinatin of constant values, of variable parameter as determined above
for ii in range(len(const_params)):

 
    #---------------Initialize Plots---------------#

    #extended growing season scenario
    #--figure 1--#

#    fig1,ax1= plt.subplots(figsize=(30,15))
#    ax1.set_xlabel('Time (days)',fontsize=40, fontname='Times New Roman')
#    ax1.set_ylabel('NUE (umol CO2/g N s)',fontsize=40, fontname='Times New Roman')
#    ax1.set_xlim([0,365])
#    ax1.set_title('NUE During an Extended Summer in a Dry Meadow', fontname='Times New Roman',fontsize=40,fontweight='bold')
#
#    #--figure 31--#
#
#    fig31,ax31= plt.subplots(figsize=(30,15))
#    ax31.set_xlabel('Time (days)',fontsize=40, fontname='Times New Roman')
#    ax31.set_ylabel('NUE (umol CO2/g N s)',fontsize=40, fontname='Times New Roman')
#    ax31.set_xlim([0,365])
#    ax31.set_title('NUE During an Extended Summer in a Wet Meadow', fontname='Times New Roman',fontsize=40,fontweight='bold')
#
#
#    #Extended growing season scenario just longer growing season
#    #--figure 22--#
#
#    fig22,ax22= plt.subplots(figsize=(30,15))
#    ax22.set_xlabel('Time (days)',fontsize=40, fontname='Times New Roman')
#    ax22.set_ylabel('NUE (umol CO2/g N s)',fontsize=40, fontname='Times New Roman')
#    ax22.set_xlim([0,365])
#    ax22.set_title('NUE During a Longer Growing Season in a Dry Meadow', fontname='Times New Roman',fontsize=40,fontweight='bold')
#
#    #--figure 30--#
#    fig30,ax30= plt.subplots(figsize=(30,15))
#    ax30.set_xlabel('Time (days)',fontsize=40, fontname='Times New Roman')
#    ax30.set_ylabel('NUE (umol CO2/g N s)',fontsize=40, fontname='Times New Roman')
#    ax30.set_xlim([0,365])
#    ax30.set_title('NUE During a Longer Growing Season in a Wet Meadow', fontname='Times New Roman',fontsize=40,fontweight='bold')
#
#
#    #Extended growing season scenario just higher temps
#    #--figure 2--#
#         
#    fig2,ax2= plt.subplots(figsize=(30,15))   
#    ax2.set_xlabel('Time (days)',fontsize=40, fontname='Times New Roman')
#    ax2.set_ylabel('NUE (umol CO2/g N s)', fontsize=40, fontname='Times New Roman')  
#    ax2.set_title('NUE During a Warmer Growing Season in a Dry Meadow', fontname='Times New Roman',fontsize=40,fontweight='bold')
#    ax2.set_xlim([0,365])
#
#
#    #--figure 32--#
#         
#    fig32,ax32= plt.subplots(figsize=(30,15))   
#    ax32.set_xlabel('Time (days)',fontsize=40, fontname='Times New Roman')
#    ax32.set_ylabel('NUE (umol CO2/g N s)', fontsize=40, fontname='Times New Roman')  
#    ax32.set_title('NUE During a Warmer Growing Season in a Wet Meadow', fontname='Times New Roman',fontsize=40,fontweight='bold')
#    ax32.set_xlim([0,365])
#  
    
    #extended growing season scenario

    #-----figure 36----#
    #put in correct ax value (e.g. axA, axB)
    fig36,ax36= plt.subplots(figsize=(12,7))
    ax36.set_xlabel('Time (days)',fontsize=20, fontname='Times New Roman')
    ax36.set_ylabel('Assimilation ($\mu$mol CO$_2$/m$^2$s)', fontsize=20, fontname='Times New Roman')  
#    ax36.set_title('Assimilation During an Extended Summer', fontname='Times New Roman',fontsize=30)
    ax36.set_xlim([0,365])


#    #-----figure 3----#
#                   
#    #put in correct ax value (e.g. axA, axB)
#    fig3,ax3= plt.subplots(figsize=(30,15))
#    ax3.set_xlabel('Time (days)',fontsize=40, fontname='Times New Roman')
#    ax3.set_ylabel('Assimilation (umol CO2/m2s)', fontsize=40, fontname='Times New Roman')  
#    ax3.set_title('Assimilation During an Extended Summer in a Dry Meadow', fontname='Times New Roman',fontsize=40,fontweight='bold')
#    ax3.set_xlim([0,365])
#
#
#    #-----figure 33---#
#                   
#    #put in correct ax value (e.g. axA, axB)
#    fig33,ax33= plt.subplots(figsize=(30,15))
#    ax33.set_xlabel('Time (days)',fontsize=40, fontname='Times New Roman')
#    ax33.set_ylabel('Assimilation (umol CO2/m2s)', fontsize=40, fontname='Times New Roman')  
#    ax33.set_title('Assimilation During an Extended Summer in a Wet Meadow', fontname='Times New Roman',fontsize=40,fontweight='bold')
#    ax33.set_xlim([0,365])
#
#
#
#    #Extended growing season scenario just longer growing season
#    #--figure 23--#
#         
#    fig23,ax23= plt.subplots(figsize=(30,15))   
#    ax23.set_xlabel('Time (days)',fontsize=40, fontname='Times New Roman')
#    ax23.set_ylabel('Assimilation (umol CO2/m2s)', fontsize=40, fontname='Times New Roman')  
#    ax23.set_title('Assimilation During a Longer Growing Season in a Dry Meadow', fontname='Times New Roman',fontsize=40,fontweight='bold')
#    ax23.set_xlim([0,365])
#
#
#    #--figure 34--#
#         
#    fig34,ax34= plt.subplots(figsize=(30,15))   
#    ax34.set_xlabel('Time (days)',fontsize=40, fontname='Times New Roman')
#    ax34.set_ylabel('Assimilation (umol CO2/m2s)', fontsize=40, fontname='Times New Roman')  
#    ax34.set_title('Assimilation During a Longer Growing Season in a Wet Meadow', fontname='Times New Roman',fontsize=40,fontweight='bold')
#    ax34.set_xlim([0,365])
#
#    #-----figure 24----#
#                   
#    #put in correct ax value (e.g. axA, axB)
#    fig24,ax24= plt.subplots(figsize=(30,15))
#    ax24.set_xlabel('Time (days)',fontsize=40, fontname='Times New Roman')
#    ax24.set_ylabel('Assimilation (umol CO2/m2s)', fontsize=40, fontname='Times New Roman')  
#    ax24.set_title('Assimilation During a Warmer Growing Season in a Dry Meadow', fontname='Times New Roman',fontsize=40,fontweight='bold')
#    ax24.set_xlim([0,365])
#
#
#    #-----figure 35----#
#                   
#    #put in correct ax value (e.g. axA, axB)
#    fig35,ax35= plt.subplots(figsize=(30,15))
#    ax35.set_xlabel('Time (days)',fontsize=40, fontname='Times New Roman')
#    ax35.set_ylabel('Assimilation (umol CO2/m2s)', fontsize=40, fontname='Times New Roman')  
#    ax35.set_title('Assimilation During a Warmer Growing Season in a Wet Meadow', fontname='Times New Roman',fontsize=40,fontweight='bold')
#    ax35.set_xlim([0,365])
#






    #Environmental Conditions During a Regular Growing Season
     #-----figure 4----#
          
    
    fig4,(ax4A,ax4B,ax4C) = plt.subplots(3,figsize=(9,20),sharex=True,sharey=False) 
    ax4C.set_xlabel('Time (days)',fontsize=20, fontname='Times New Roman')
    ax4A.set_ylabel('Temperature ($^\circ$C)', fontsize=20, fontname='Times New Roman')  
#    ax4A.set_title('Niwot Ridge Microclimates', fontname='Times New Roman',fontsize=30)
    ax4C.set_xlim([0,365])
    ax4B.set_ylabel('VWC (m$^3$/m$^3$)', fontsize=20, fontname='Times New Roman')  
    ax4C.set_ylabel('Foliar Nitrogen (g N/m$^2$)', fontsize=20, fontname='Times New Roman')  
    fig4.tight_layout()
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.2)  
 
    
    #Environmental Conditions in a Regular and Extended Summer
   #--figure 19--#
         
#    fig19,ax19= plt.subplots(figsize=(30,15))
#    ax19.set_xlabel('Time (days)',fontsize=40, fontname='Times New Roman')
#    ax19.set_ylabel('Temperature ($^\circ$C)', fontsize=40, fontname='Times New Roman')  
#    ax19.set_title('Surface Temperature', fontname='Times New Roman',fontsize=40,fontweight='bold')
#    ax19.set_xlim([0,365])
# 
#
#     #-----figure 20----#
#
#    fig20,ax20 = plt.subplots(figsize=(30,15))
#    ax20.set_xlabel('Time (days)',fontsize=40, fontname='Times New Roman')
#    ax20.set_ylabel('Volumetric Water Content (m3/m3)', fontsize=40, fontname='Times New Roman')  
#    ax20.set_title('Volumetric Water Content', fontname='Times New Roman',fontsize=40,fontweight='bold')
#    ax20.set_xlim([0,365])
#    
#     #-----figure 21----#
#          
#    fig21,ax21 = plt.subplots(figsize=(30,15))
#    ax21.set_xlabel('Time (days)',fontsize=40, fontname='Times New Roman')
#    ax21.set_ylabel('Leaf Nitrogen Content (g N/m2)', fontsize=40, fontname='Times New Roman')  
#    ax21.set_title('Leaf Nitrogen Content', fontname='Times New Roman',fontsize=40,fontweight='bold')
#    ax21.set_xlim([0,365])








    #Regular growing season scenario
    
    #with specific microclimates
    #----figure 7---#

    fig7,(ax7A,ax7B,ax7C) = plt.subplots(3,figsize=(9,20),sharex=True,sharey=False) 

    ax7C.set_xlabel('Time (days)',fontsize=18, fontname='Times New Roman')
    ax7A.set_ylabel('Assimilation ($\mu$mol CO$_2$/m$^2$ s)', fontsize=18, fontname='Times New Roman')  
#    ax7A.set_title('Assimilation, NUE, and WUE During the Growing Season', fontname='Times New Roman',fontsize=25)
    ax7C.set_xlim([0,365])   
    ax7B.set_ylabel('NUE ($\mu$mol CO$_2$/g N s)', fontsize=18, fontname='Times New Roman')  
    ax7C.set_ylabel('WUE ($\mu$mol CO$_2$/mmol H$_2$O)', fontsize=18, fontname='Times New Roman')  
    fig7.tight_layout()
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.2)  

   

    #in identical microclimates
    #----figure 8---#

    fig8,(ax8A,ax8B) = plt.subplots(2,figsize=(9,20),sharex=True,sharey=False) 
    ax8B.set_xlabel('Time (days)',fontsize=20, fontname='Times New Roman')
    ax8A.set_ylabel('Assimilation ($\mu$mol CO$_2$/m$^2$ s)', fontsize=20, fontname='Times New Roman')  
    ax8A.set_title('Dry Meadow', fontname='Times New Roman',fontsize=25)
    ax8B.set_title('Wet Meadow', fontname='Times New Roman',fontsize=25)
    ax8B.set_xlim([0,365])        
    ax8B.set_ylabel('Assimilation ($\mu$mol CO$_2$/m$^2$ s)', fontsize=20, fontname='Times New Roman')  
   


#    #----figure 10---#
#    fig10,ax10=plt.subplots(figsize=(30,15))
#    ax10.set_xlabel('Time (days)',fontsize=40, fontname='Times New Roman')
#    ax10.set_ylabel('NUE (umol CO2/g N s)', fontsize=40, fontname='Times New Roman')  
#    ax10.set_title('NUE During the Growing Season in a Dry Meadow', fontname='Times New Roman',fontsize=40,fontweight='bold')
#    ax10.set_xlim([0,365])    
#
#    #----figure 11---#
#    fig11,ax11=plt.subplots(figsize=(30,15))
#    ax11.set_xlabel('Time (days)',fontsize=40, fontname='Times New Roman')
#    ax11.set_ylabel('NUE (umol CO2/g N s)', fontsize=40, fontname='Times New Roman')  
#    ax11.set_title('NUE During the Growing Season in a Wet Meadow', fontname='Times New Roman',fontsize=40,fontweight='bold')
#    ax11.set_xlim([0,365])  

#
#   #----figure 12---#
#    fig12,ax12= plt.subplots(figsize=(30,15))
#    ax12.set_xlabel('Time (days)',fontsize=40, fontname='Times New Roman')
#    ax12.set_ylabel('WUE (umol CO2/mmol H2O)', fontsize=40, fontname='Times New Roman')  
#    ax12.set_title('WUE During the Growing Season in a Moist Meadow Climate', fontname='Times New Roman',fontsize=40,fontweight='bold')
#    ax12.set_xlim([0,365])    

  

#TESTING DIFFERENT TEMPS/MOISTURES/Nitrogen

    #----figure 13---#

    fig13,(ax13A,ax13B,ax13C) = plt.subplots(3,figsize=(9,20),sharey=True,sharex=False) 


#    ax13B.set_title('Assimilation as a Function of the Environment',fontsize=30, fontname='Times New Roman')
    ax13A.set_xlabel('Temperature ($^\circ$C)',fontsize=19, fontname='Times New Roman')
    ax13B.set_ylabel('Assimilation ($\mu$mol CO$_2$/m$^2$ s)', fontsize=19, fontname='Times New Roman')     
    ax13A.set_ylim([0,22])   
    ax13B.set_xlabel('Volumetric Water Content (m$^3$/m$^3$)',fontsize=19, fontname='Times New Roman')
    ax13B.set_xlim([0.05,0.4])
    ax13C.set_xlabel('Foliar Nitrogen Content (g N/m$^2$)',fontsize=19, fontname='Times New Roman')
    ax13A.tick_params(axis='x', labelsize=15)
    ax13B.tick_params(axis='x', labelsize=15)
    ax13C.tick_params(axis='x', labelsize=15)
    fig13.tight_layout()
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.3)  



#    #----figure 16---#
#    fig16,ax16= plt.subplots(figsize=(15,30))    
#    ax16.set_title('NUE as a Function of Temperature',fontsize=40, fontname='Times New Roman',fontweight='bold')
#    ax16.set_xlabel('Temperature ($^\circ$C)',fontsize=40, fontname='Times New Roman')
#    ax16.set_ylabel('NUE (umol CO2/g N s)', fontsize=40, fontname='Times New Roman')  
#    ax16.set_ylim([0,4])         
#
#    #----figure 17---# 
#    fig17,ax17= plt.subplots(figsize=(15,30))    
#    ax17.set_title('NUE as a Function of Soil Moisture',fontsize=40, fontname='Times New Roman',fontweight='bold')
#    ax17.set_xlabel('Volumetric Water Content (m3/m3)',fontsize=40, fontname='Times New Roman')
#    ax17.set_ylabel('NUE (umol CO2/g N s)', fontsize=40, fontname='Times New Roman')  
#    ax17.set_xlim([0.05,0.4])
#    ax17.set_ylim([0,4])     
#
#    #----figure 18---#     
#    fig18,ax18= plt.subplots(figsize=(15,30))    
#    ax18.set_title('NUE as a Function of Leaf Nitrogen Content',fontsize=40, fontname='Times New Roman',fontweight='bold')
#    ax18.set_xlabel('Nitrogen Content (gN/m2)',fontsize=40, fontname='Times New Roman')
#    ax18.set_ylabel('NUE (umol CO2/g N s)', fontsize=40, fontname='Times New Roman')  
#    ax18.set_ylim([0,4])     
#
#
#    #----figure 27---#
#    fig27,ax27= plt.subplots(figsize=(15,30))    
#    ax27.set_title('WUE as a Function of Temperature',fontsize=40, fontname='Times New Roman',fontweight='bold')
#    ax27.set_xlabel('Temperature ($^\circ$C)',fontsize=40, fontname='Times New Roman')
#    ax27.set_ylabel('WUE (umol CO2/mmol H2O)', fontsize=40, fontname='Times New Roman')  
#    ax27.set_ylim([0,9])         
#
#    #----figure 28---# 
#    fig28,ax28= plt.subplots(figsize=(15,30))    
#    ax28.set_title('WUE as a Function of Soil Moisture',fontsize=40, fontname='Times New Roman',fontweight='bold')
#    ax28.set_xlabel('Volumetric Water Content (m3/m3)',fontsize=40, fontname='Times New Roman')
#    ax28.set_ylabel('WUE (umol CO2/mmol H2O)', fontsize=40, fontname='Times New Roman')  
#    ax28.set_xlim([0.05,0.4])
#    ax28.set_ylim([0,9]) 
#
#    #----figure 29---#     
#    fig29,ax29= plt.subplots(figsize=(15,30))    
#    ax29.set_title('WUE as a Function of Leaf Nitrogen Content',fontsize=40, fontname='Times New Roman',fontweight='bold')
#    ax29.set_xlabel('Nitrogen Content (gN/m2)',fontsize=40, fontname='Times New Roman')
#    ax29.set_ylabel('WUE (umol CO2/mmol H2O)', fontsize=40, fontname='Times New Roman')  
#    ax29.set_ylim([0,9]) 


 





  
    #---------------Run through time series---------------#
    
    days=np.linspace(1,365,365)
    
    nue_d_inter_stats=[[],[],[]]
    wue_d_inter_stats=[[],[],[]]
    A_d_inter_stats=[[],[],[]]

    nue_w_inter_stats=[[],[],[]]
    wue_w_inter_stats=[[],[],[]]
    A_w_inter_stats=[[],[],[]]



    wue_d_dmclimate_gsl_stats=[[],[],[]]
    nue_d_dmclimate_gsl_stats=[[],[],[]]        
    A_d_dmclimate_gsl_stats=[[],[],[]]        

    wue_d_wmclimate_gsl_stats=[[],[],[]]
    nue_d_wmclimate_gsl_stats=[[],[],[]]      
    A_d_wmclimate_gsl_stats=[[],[],[]] 


    wue_d_dy_stats=[[],[],[]]
    nue_d_dy_stats=[[],[],[]]        
    A_d_dy_stats=[[],[],[]]   
    
    wue_w_dy_stats=[[],[],[]]
    nue_w_dy_stats=[[],[],[]]        
    A_w_dy_stats=[[],[],[]]   

    wue_d_dmclimate_dy_stats=[[],[],[]]
    nue_d_dmclimate_dy_stats=[[],[],[]]        
    A_d_dmclimate_dy_stats=[[],[],[]]        

    wue_d_wmclimate_dy_stats=[[],[],[]]
    nue_d_wmclimate_dy_stats=[[],[],[]]        
    A_d_wmclimate_dy_stats=[[],[],[]]  
    
    wue_d_dmclimate_temps_stats=[[],[],[]]
    nue_d_dmclimate_temps_stats=[[],[],[]]       
    A_d_dmclimate_temps_stats=[[],[],[]]        

    wue_d_wmclimate_temps_stats=[[],[],[]]
    nue_d_wmclimate_temps_stats=[[],[],[]]       
    A_d_wmclimate_temps_stats=[[],[],[]]  


    wue_w_dmclimate_gsl_stats=[[],[],[]]
    nue_w_dmclimate_gsl_stats=[[],[],[]]        
    A_w_dmclimate_gsl_stats=[[],[],[]]        

    wue_w_wmclimate_gsl_stats=[[],[],[]]
    nue_w_wmclimate_gsl_stats=[[],[],[]]      
    A_w_wmclimate_gsl_stats=[[],[],[]] 


    wue_w_dmclimate_dy_stats=[[],[],[]]
    nue_w_dmclimate_dy_stats=[[],[],[]]        
    A_w_dmclimate_dy_stats=[[],[],[]]        

    wue_w_wmclimate_dy_stats=[[],[],[]]
    nue_w_wmclimate_dy_stats=[[],[],[]]        
    A_w_wmclimate_dy_stats=[[],[],[]]  
    
    wue_w_dmclimate_temps_stats=[[],[],[]]
    nue_w_dmclimate_temps_stats=[[],[],[]]       
    A_w_dmclimate_temps_stats=[[],[],[]]        

    wue_w_wmclimate_temps_stats=[[],[],[]]
    nue_w_wmclimate_temps_stats=[[],[],[]]       
    A_w_wmclimate_temps_stats=[[],[],[]]  





    nue_d_inter_vartemp_stats=[[],[],[]]
    wue_d_inter_vartemp_stats=[[],[],[]]
    A_d_inter_vartemp_stats=[[],[],[]]

    nue_d_inter_varvwc_stats=[[],[],[]]
    wue_d_inter_varvwc_stats=[[],[],[]]
    A_d_inter_varvwc_stats=[[],[],[]]

    nue_d_inter_varna_stats=[[],[],[]]
    wue_d_inter_varna_stats=[[],[],[]]
    A_d_inter_varna_stats=[[],[],[]]
    
    nue_w_inter_vartemp_stats=[[],[],[]]
    wue_w_inter_vartemp_stats=[[],[],[]]
    A_w_inter_vartemp_stats=[[],[],[]]

    nue_w_inter_varvwc_stats=[[],[],[]]
    wue_w_inter_varvwc_stats=[[],[],[]]
    A_w_inter_varvwc_stats=[[],[],[]]

    nue_w_inter_varna_stats=[[],[],[]]
    wue_w_inter_varna_stats=[[],[],[]]
    A_w_inter_varna_stats=[[],[],[]]    

    wue_w_inter_mmclimate_stats=[[],[],[]]
    nue_w_inter_mmclimate_stats=[[],[],[]]
    A_w_inter_mmclimate_stats=[[],[],[]]
    wue_d_inter_mmclimate_stats=[[],[],[]]
    nue_d_inter_mmclimate_stats=[[],[],[]]
    A_d_inter_mmclimate_stats=[[],[],[]]    

    wue_w_inter_dmclimate_stats=[[],[],[]]
    nue_w_inter_dmclimate_stats=[[],[],[]]
    A_w_inter_dmclimate_stats=[[],[],[]]
    wue_d_inter_dmclimate_stats=[[],[],[]]
    nue_d_inter_dmclimate_stats=[[],[],[]]
    A_d_inter_dmclimate_stats=[[],[],[]]    

    wue_w_inter_wmclimate_stats=[[],[],[]]
    nue_w_inter_wmclimate_stats=[[],[],[]]
    A_w_inter_wmclimate_stats=[[],[],[]]
    wue_d_inter_wmclimate_stats=[[],[],[]]
    nue_d_inter_wmclimate_stats=[[],[],[]]
    A_d_inter_wmclimate_stats=[[],[],[]]   
    
    nue_d_inter_stats_smooth=[]
    wue_d_inter_stats_smooth=[]
    A_d_inter_stats_smooth=[]
    nue_w_inter_stats_smooth=[]
    wue_w_inter_stats_smooth=[]
    A_w_inter_stats_smooth=[]

           
     
    wue_d_dmclimate_gsl_stats_smooth=[]
    nue_d_dmclimate_gsl_stats_smooth=[]    
    A_d_dmclimate_gsl_stats_smooth=[]    

    wue_d_wmclimate_gsl_stats_smooth=[]
    nue_d_wmclimate_gsl_stats_smooth=[]    
    A_d_wmclimate_gsl_stats_smooth=[]


    wue_d_dmclimate_dy_stats_smooth=[]
    nue_d_dmclimate_dy_stats_smooth=[]      
    A_d_dmclimate_dy_stats_smooth=[]     

    wue_d_dy_stats_smooth=[]
    nue_d_dy_stats_smooth=[]      
    A_d_dy_stats_smooth=[]    

    wue_w_dy_stats_smooth=[]
    nue_w_dy_stats_smooth=[]      
    A_w_dy_stats_smooth=[]    

    wue_d_wmclimate_dy_stats_smooth=[]
    nue_d_wmclimate_dy_stats_smooth=[]      
    A_d_wmclimate_dy_stats_smooth=[]
    
    wue_d_dmclimate_temps_stats_smooth=[]
    nue_d_dmclimate_temps_stats_smooth=[]      
    A_d_dmclimate_temps_stats_smooth=[]     

    wue_d_wmclimate_temps_stats_smooth=[]
    nue_d_wmclimate_temps_stats_smooth=[]    
    A_d_wmclimate_temps_stats_smooth=[] 


    wue_w_dmclimate_gsl_stats_smooth=[]
    nue_w_dmclimate_gsl_stats_smooth=[]      
    A_w_dmclimate_gsl_stats_smooth=[]      

    wue_w_wmclimate_gsl_stats_smooth=[]
    nue_w_wmclimate_gsl_stats_smooth=[]   
    A_w_wmclimate_gsl_stats_smooth=[]


    wue_w_dmclimate_dy_stats_smooth=[]
    nue_w_dmclimate_dy_stats_smooth=[]       
    A_w_dmclimate_dy_stats_smooth=[]       

    wue_w_wmclimate_dy_stats_smooth=[]
    nue_w_wmclimate_dy_stats_smooth=[]       
    A_w_wmclimate_dy_stats_smooth=[]
    
    wue_w_dmclimate_temps_stats_smooth=[]
    nue_w_dmclimate_temps_stats_smooth=[]      
    A_w_dmclimate_temps_stats_smooth=[]   

    wue_w_wmclimate_temps_stats_smooth=[]
    nue_w_wmclimate_temps_stats_smooth=[]       
    A_w_wmclimate_temps_stats_smooth=[]
    
    
    
    
    nue_d_inter_vartemp_stats_smooth=[]
    wue_d_inter_vartemp_stats_smooth=[]
    A_d_inter_vartemp_stats_smooth=[]

    nue_d_inter_varvwc_stats_smooth=[]
    wue_d_inter_varvwc_stats_smooth=[]
    A_d_inter_varvwc_stats_smooth=[]

    nue_d_inter_varna_stats_smooth=[]
    wue_d_inter_varna_stats_smooth=[]
    A_d_inter_varna_stats_smooth=[]

    nue_w_inter_vartemp_stats_smooth=[]
    wue_w_inter_vartemp_stats_smooth=[]
    A_w_inter_vartemp_stats_smooth=[]

    nue_w_inter_varvwc_stats_smooth=[]
    wue_w_inter_varvwc_stats_smooth=[]
    A_w_inter_varvwc_stats_smooth=[]

    nue_w_inter_varna_stats_smooth=[]
    wue_w_inter_varna_stats_smooth=[]
    A_w_inter_varna_stats_smooth=[]
 
    wue_w_inter_mmclimate_stats_smooth=[]
    nue_w_inter_mmclimate_stats_smooth=[]
    A_w_inter_mmclimate_stats_smooth=[]
    wue_d_inter_mmclimate_stats_smooth=[]
    nue_d_inter_mmclimate_stats_smooth=[]
    A_d_inter_mmclimate_stats_smooth=[]    


    wue_w_inter_dmclimate_stats_smooth=[]
    nue_w_inter_dmclimate_stats_smooth=[]
    A_w_inter_dmclimate_stats_smooth=[]
    wue_d_inter_dmclimate_stats_smooth=[]
    nue_d_inter_dmclimate_stats_smooth=[]
    A_d_inter_dmclimate_stats_smooth=[]    


    wue_w_inter_wmclimate_stats_smooth=[]
    nue_w_inter_wmclimate_stats_smooth=[]
    A_w_inter_wmclimate_stats_smooth=[]
    wue_d_inter_wmclimate_stats_smooth=[]
    nue_d_inter_wmclimate_stats_smooth=[]
    A_d_inter_wmclimate_stats_smooth=[]            
    
    tot_temps_dm=[]
    tot_temps_wm=[]
    
    na_lst=[]


    for time in range(365):
        
       
        #arrays
#dry meadow        
        wue_d_inter=[]
        nue_d_inter=[]
        gsw_d_inter=[]
        A_d_inter=[]
        E_d_inter=[]
        vpd_d_inter=[]
    
        
        wue_d_dmclimate_gsl=[]
        nue_d_dmclimate_gsl=[]        
        A_d_dmclimate_gsl=[]        

        wue_d_wmclimate_gsl=[]
        nue_d_wmclimate_gsl=[]        
        A_d_wmclimate_gsl=[]  


        wue_d_dmclimate_dy=[]
        nue_d_dmclimate_dy=[]        
        A_d_dmclimate_dy=[]   
        
        wue_d_dy=[]
        nue_d_dy=[]        
        A_d_dy=[]           

        wue_w_dy=[]
        nue_w_dy=[]        
        A_w_dy=[]                   

        wue_d_wmclimate_dy=[]
        nue_d_wmclimate_dy=[]        
        A_d_wmclimate_dy=[]  
        
        wue_d_dmclimate_temps=[]
        nue_d_dmclimate_temps=[]        
        A_d_dmclimate_temps=[]        

        wue_d_wmclimate_temps=[]
        nue_d_wmclimate_temps=[]        
        A_d_wmclimate_temps=[]   

    
#wet meadow    
        wue_w_inter=[]
        nue_w_inter=[]
        gsw_w_inter=[]
        A_w_inter=[]
        E_w_inter=[]
        vpd_w_inter=[]
        

        
        wue_w_dmclimate_gsl=[]
        nue_w_dmclimate_gsl=[]        
        A_w_dmclimate_gsl=[]        

        wue_w_wmclimate_gsl=[]
        nue_w_wmclimate_gsl=[]        
        A_w_wmclimate_gsl=[]  


        wue_w_dmclimate_dy=[]
        nue_w_dmclimate_dy=[]        
        A_w_dmclimate_dy=[]        

        wue_w_wmclimate_dy=[]
        nue_w_wmclimate_dy=[]        
        A_w_wmclimate_dy=[]  
        
        wue_w_dmclimate_temps=[]
        nue_w_dmclimate_temps=[]        
        A_w_dmclimate_temps=[]        

        wue_w_wmclimate_temps=[]
        nue_w_wmclimate_temps=[]        
        A_w_wmclimate_temps=[]          
#constant temp/vwc  
        nue_d_inter_vartemp=[]
        wue_d_inter_vartemp=[]
        A_d_inter_vartemp=[]
    
        nue_d_inter_varvwc=[]
        wue_d_inter_varvwc=[]
        A_d_inter_varvwc=[]

        nue_d_inter_varna=[]
        wue_d_inter_varna=[]
        A_d_inter_varna=[]
        
        nue_w_inter_vartemp=[]
        wue_w_inter_vartemp=[]
        A_w_inter_vartemp=[]
    
        nue_w_inter_varvwc=[]
        wue_w_inter_varvwc=[]
        A_w_inter_varvwc=[]

        nue_w_inter_varna=[]
        wue_w_inter_varna=[]
        A_w_inter_varna=[]
         
        tempss_dm=[]
        tempss_wm=[]        
        na_lst0=[]
        
        
#constant climate between plant comms

        wue_d_inter_mmclimate=[]
        nue_d_inter_mmclimate=[]
        A_d_inter_mmclimate=[]
     
        wue_w_inter_mmclimate=[]
        nue_w_inter_mmclimate=[]
        A_w_inter_mmclimate=[]
     
        wue_d_inter_dmclimate=[]
        nue_d_inter_dmclimate=[]
        A_d_inter_dmclimate=[]
     
        wue_w_inter_dmclimate=[]
        nue_w_inter_dmclimate=[]
        A_w_inter_dmclimate=[]

        wue_d_inter_wmclimate=[]
        nue_d_inter_wmclimate=[]
        A_d_inter_wmclimate=[]
     
        wue_w_inter_wmclimate=[]
        nue_w_inter_wmclimate=[]
        A_w_inter_wmclimate=[]

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
   
        
        #------constant variable params-----#
        
        chl_c=np.zeros(shape=1)+(np.mean([396,465,476])) #Chlorophyll Content of the Leaf (umol chl/m2)
        ht_c=np.zeros(shape=1)+15.0 #Temperature of the Leaf (K)
        dia_c=np.zeros(shape=1)+(np.mean([1.4,2.3,2.6])/100.) #Mean diameter or size of leaf (m)
        na_c=np.zeros(shape=1)+4.0 #leaf nitrogen (g N/ m2)
        t_c=np.zeros(shape=1)+15.0 #temp (C)


        #----climate timeseries----#    
        
        #reg year (first list); dry/hot year (second list)
        temp_dm=[surtemp_dm, surtemp_dm_dy, np.linspace(0,20,365), [np.mean(surtemp_mm[gs0_m:gsf+1])]*365, surtemp_mm,surtemp_dm,surtemp_wm,[np.mean(surtemp_mm[gs0_m:gsf+1])]*365,surtemp_dm,surtemp_wm_dy,surtemp_wm, surtemp_dm_dy,surtemp_wm_dy,surtemp_dm_dy]
        temp_wm=[surtemp_wm, surtemp_dm_dy,np.linspace(0,20,365), [np.mean(surtemp_mm[gs0_m:gsf+1])]*365, surtemp_mm,surtemp_dm,surtemp_wm,[np.mean(surtemp_mm[gs0_m:gsf+1])]*365,surtemp_dm,surtemp_wm_dy,surtemp_wm, surtemp_dm_dy,surtemp_wm_dy,surtemp_wm_dy]
        sm_dm=[vwc_dm, vwc_dm_dy,[np.max(vwc_mm[gs0_m:gsf+1])]*365,np.linspace(0.08, 0.35, 365), vwc_mm,vwc_dm,vwc_wm,[np.mean(vwc_mm[gs0_m:gsf+1])]*365, vwc_dm,vwc_wm_dy,vwc_wm, vwc_dm,vwc_wm,vwc_dm_dy]
        sm_wm=[vwc_wm, vwc_dm_dy,[np.max(vwc_mm[gs0_m:gsf+1])]*365,np.linspace(0.08, 0.35, 365), vwc_mm,vwc_dm,vwc_wm,[np.mean(vwc_mm[gs0_m:gsf+1])]*365, vwc_dm,vwc_wm_dy,vwc_wm, vwc_dm,vwc_wm,vwc_wm_dy]        
        
        #----variables used in monte carlo simulations----#
#        #Inter Specific Variation (first list) & Intra Specific Variation (second list): Dominant Species (Kobresia, Carex)
#        chl_mean=[[395.7132,475.8913],[339.5,500.3]]
#        chl_sd=[[24.410199999999975, 29.185099999999977],[90.3,217.1]]
#        dia_mean=[[1.6/100.,3.0/100.],[(2*np.sqrt(0.48/np.pi))/100.,(2*np.sqrt(5.64/np.pi))/100.]]
#        dia_sd=[[0.9/100.0,1.2/100.0],[(2*np.sqrt(0.15/np.pi))/100.,(2*np.sqrt(1.03/np.pi))/100.]]
#        na_min=[[na_dm_min_inter[time],na_wm_min_inter[time]],[na_dm_min_intra[time],na_wm_min_intra[time]]]
#        na_max=[[na_dm_max_inter[time],na_wm_max_inter[time]],[na_dm_max_intra[time],na_wm_max_intra[time]]]
#        ht_mean=[[9.183549,19.98519],[7.7,27.9]]
#        ht_sd=[[1.5,3.1],[3.0,9.1]]

            
#no intraspecific variation
        chl_mean=[[395.7132,475.8913]]*8
        chl_sd=[[24.410199999999975,29.185099999999977]]*8
        dia_mean=[[1.6/100.,3.0/100.]]*8
        dia_sd=[[0.9/100.0,1.2/100.0]]*8
        na_min=[[na_dm_min_inter[time],na_wm_min_inter[time]],[np.mean(na_mm_min_inter[gs0_m:gsf+1]),np.mean(na_mm_min_inter[gs0_m:gsf+1])],[na_dm_min_inter_dy[time],na_dm_min_inter_dy[time]],[na_mm_min_inter[time],na_mm_min_inter[time]],[na_dm_min_inter[time],na_dm_min_inter[time]],[na_wm_min_inter[time],na_wm_min_inter[time]],[na_wm_min_inter_dy[time],na_wm_min_inter_dy[time]],[na_dm_min_inter_dy[time],na_wm_min_inter_dy[time]]]
        na_max=[[na_dm_max_inter[time],na_wm_max_inter[time]],[np.mean(na_mm_max_inter[gs0_m:gsf+1]),np.mean(na_mm_max_inter[gs0_m:gsf+1])],[na_dm_max_inter_dy[time],na_dm_max_inter_dy[time]],[na_mm_max_inter[time],na_mm_max_inter[time]],[na_dm_max_inter[time],na_dm_max_inter[time]],[na_wm_max_inter[time],na_wm_max_inter[time]],[na_wm_max_inter_dy[time],na_wm_max_inter_dy[time]],[na_dm_max_inter_dy[time],na_wm_max_inter_dy[time]]]
        ht_mean=[[9.183549,19.98519]]*8
        ht_sd=[[1.5,3.1]]*8
                   
    
        
        #Non-CWM Inter Specific Variation (first list) & Intra Specific Variation (second list): Dominant Species (Kobresia, Carex)
#        chl_mean=[[393.3512,452.7842],[339.5,500.3]]
#        chl_sd=[[44.54317, 31.35734],[90.3,217.1]]
#        dia_mean=[[(2*np.sqrt(3.192766/np.pi))/100.,(2*np.sqrt(5.769438/np.pi))/100.],[(2*np.sqrt(0.48/np.pi))/100.,(2*np.sqrt(5.64/np.pi))/100.]]
#        dia_sd=[[(2*np.sqrt(0.8345044/np.pi))/100.,(2*np.sqrt(1.186597/np.pi))/100.],[(2*np.sqrt(0.15/np.pi))/100.,(2*np.sqrt(1.03/np.pi))/100.]]
#        na_min=[[na_dm_min_inter[time],na_wm_min_inter[time]],[na_dm_min_intra[time],na_wm_min_intra[time]]]
#        na_max=[[na_dm_max_inter[time],na_wm_max_inter[time]],[na_dm_max_intra[time],na_wm_max_intra[time]]]
#        ht_mean=[[11.81,16.09],[7.7,27.9]]
#        ht_sd=[[3.16,2.23],[3.0,9.1]]
   
        for iiii in range(14):
            for iii in range(8):        
            
            #---------------Import Variable Parameter Arrays from Leaf Parameter File: Regular Year---------------#
                params=monte_carlo(chl_mean[iii], chl_sd[iii], dia_mean[iii], dia_sd[iii], na_min[iii], na_max[iii], ht_mean[iii], ht_sd[iii])        
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
                            
                            
            
                        
                        #number of simulations per meadow type:
                        m_sim=len(params) #meadow simulations 
                        
              
                        
                        if xx==0: 
                
                            #------calculate vapor pressure-----#
                            pa_v=611*np.exp((17.27*temp_dm[iiii][time])/(temp_dm[iiii][time]+237.3)) #saturation vapor pressure of air (Pa)
                            ea_str=pa_con_atmfrac(pa_v,3528) #saturation vapor pressure of air (Pa-->umol h20/mol air)
                            ea=rh*ea_str #vapor pressure (umol h2O/mol air)                
                            
        
                
                            #correct for leaf temperatures using leaf height
                   
                            t_diff=18-0.4*ht
                        
                            tl=temp_dm[iiii][time]+t_diff
                            
                            
        
                            #soil depth
                            z=0.2
                            
        
                            
                            #---------------Photosynthesis Function---------------#
                        
                            #alter this line of code for when implementing different photosynthesis functions
                            wue, nue, A, E, cs, ci, gsw, gs, gbw, gb, gm, cc,dd =photo(tk_25,ekc,eko,etau,ev,ej,toptv,toptj,na, qeff, PAR,tl,ea,chl,ij,kc25,ko25,o,ca,rh,m,a,frnr,flnr,ra,jm,g0,b,dia,u,q,vwc_min,vwc_max,sm_dm[iiii][time],z)
             
                            if isinstance(wue, np.ndarray):
                                wue=wue[0]
                    
                            if isinstance(nue, np.ndarray):
                                nue=nue[0]
                        
                            if isinstance(A, np.ndarray):
                                A=A[0]        
                        
                            if isinstance(wue, list):
                                wue=wue[0]
                        
                            if isinstance(nue, list):
                                nue=nue[0]
                        
                            if isinstance(A, list):
                                A=A[0]                        
                        
                            if iiii==0:
                                if iii==0:
                                    wue_d_inter+=[wue]
                                    nue_d_inter+=[nue]
                                    gsw_d_inter+=[gsw]
                                    A_d_inter+=[A]
                                    E_d_inter+=[E]
                                    vpd_d_inter+=[dd]
                                    tempss_dm+=[tl]

                            
                            if iiii==8:    
                                if iii==2:
                                    wue_d_dmclimate_gsl+=[wue]
                                    nue_d_dmclimate_gsl+=[nue]
                                    A_d_dmclimate_gsl+=[A]

                            if iiii==10:    
                                if iii==6:
                                    wue_d_wmclimate_gsl+=[wue]
                                    nue_d_wmclimate_gsl+=[nue]
                                    A_d_wmclimate_gsl+=[A]

                            if iiii==11:    
                                if iii==4:
                                    wue_d_dmclimate_temps+=[wue]
                                    nue_d_dmclimate_temps+=[nue]
                                    A_d_dmclimate_temps+=[A]

                            if iiii==12:    
                                if iii==5:
                                    wue_d_wmclimate_temps+=[wue]
                                    nue_d_wmclimate_temps+=[nue]
                                    A_d_wmclimate_temps+=[A]
                                    
                            if iiii==1:                                
                                if iii==2:
                                    wue_d_dmclimate_dy+=[wue]
                                    nue_d_dmclimate_dy+=[nue]
                                    A_d_dmclimate_dy+=[A]

                            if iiii==9:                                
                                if iii==6:
                                    wue_d_wmclimate_dy+=[wue]
                                    nue_d_wmclimate_dy+=[nue]
                                    A_d_wmclimate_dy+=[A]        

                            if iiii==2:                                
                                if iii==1:
                                    wue_d_inter_vartemp+=[wue]
                                    nue_d_inter_vartemp+=[nue]
                                    A_d_inter_vartemp+=[A]
        
                            if iiii==3:                                
                                if iii==1:
                                    wue_d_inter_varvwc+=[wue]
                                    nue_d_inter_varvwc+=[nue]
                                    A_d_inter_varvwc+=[A]
                                    
                            if iiii==7:                                
                                if iii==3:
                                    wue_d_inter_varna+=[wue]
                                    nue_d_inter_varna+=[nue]
                                    A_d_inter_varna+=[A]  
                                    na_lst0+=[na] 
                                   
 
                            if iiii==4:                                
                                if iii==3:
                                    wue_d_inter_mmclimate+=[wue]
                                    nue_d_inter_mmclimate+=[nue]
                                    A_d_inter_mmclimate+=[A]
                                    tempss_dm+=[tl]

                            if iiii==5:                                
                                if iii==4:
                                    wue_d_inter_dmclimate+=[wue]
                                    nue_d_inter_dmclimate+=[nue]
                                    A_d_inter_dmclimate+=[A]

                            if iiii==6:                                
                                if iii==5:
                                    wue_d_inter_wmclimate+=[wue]
                                    nue_d_inter_wmclimate+=[nue]
                                    A_d_inter_wmclimate+=[A]   
                                    
                            if iiii==13:                                
                                if iii==7:
                                    wue_d_dy+=[wue]
                                    nue_d_dy+=[nue]
                                    A_d_dy+=[A]                                        
                                    
                                    
                        elif xx==1:
            
                            #------calculate vapor pressure-----#
                            pa_v=611*np.exp((17.27*temp_wm[iiii][time])/(temp_wm[iiii][time]+237.3)) #saturation vapor pressure of air (Pa)
                            ea_str=pa_con_atmfrac(pa_v,3528) #saturation vapor pressure of air (Pa-->umol h20/mol air)
                            ea=rh*ea_str #vapor pressure (umol h2O/mol air)                
                
                
                            #correct for leaf temperatures using leaf height
                   
                            t_diff=18-0.4*ht
                            
                            if t_diff<0.0:
                                t_diff=0.0
                        
                            tl=temp_wm[iiii][time]+t_diff        
                            
        
        

                            
                            #soil depth
                            z=0.4
                            
        
                            #---------------Photosynthesis Function---------------#
        
                            #alter this line of code for when implementing different photosynthesis functions
                            wue, nue, A, E, cs, ci, gsw, gs, gbw, gb, gm, cc,dd =photo(tk_25,ekc,eko,etau,ev,ej,toptv,toptj,na, qeff, PAR,tl,ea,chl,ij,kc25,ko25,o,ca,rh,m,a,frnr,flnr,ra,jm,g0,b,dia,u,q,vwc_min,vwc_max,sm_wm[iiii][time],z)
             

                              
                            if isinstance(wue, np.ndarray):
                                wue=wue[0]
                    
                            if isinstance(nue, np.ndarray):
                                nue=nue[0]
                        
                            if isinstance(A, np.ndarray):
                                A=A[0]        
                        
                            if isinstance(wue, list):
                                wue=wue[0]
                        
                            if isinstance(nue, list):
                                nue=nue[0]
                        
                            if isinstance(A, list):
                                A=A[0]       
                         
                           
                            if iiii==0:
                                if iii==0:    
                                    wue_w_inter+=[wue]
                                    nue_w_inter+=[nue]
                                    gsw_w_inter+=[gsw]
                                    A_w_inter+=[A]
                                    E_w_inter+=[E]
                                    vpd_w_inter+=[dd]
                                    tempss_wm+=[tl]
                            
                            
                            if iiii==8:    
                                if iii==2:
                                    wue_w_dmclimate_gsl+=[wue]
                                    nue_w_dmclimate_gsl+=[nue]
                                    A_w_dmclimate_gsl+=[A]

                            if iiii==10:    
                                if iii==6:
                                    wue_w_wmclimate_gsl+=[wue]
                                    nue_w_wmclimate_gsl+=[nue]
                                    A_w_wmclimate_gsl+=[A]

                            if iiii==11:    
                                if iii==4:
                                    wue_w_dmclimate_temps+=[wue]
                                    nue_w_dmclimate_temps+=[nue]
                                    A_w_dmclimate_temps+=[A]

                            if iiii==12:    
                                if iii==5:
                                    wue_w_wmclimate_temps+=[wue]
                                    nue_w_wmclimate_temps+=[nue]
                                    A_w_wmclimate_temps+=[A]
                                    
                            if iiii==1:                                
                                if iii==2:
                                    wue_w_dmclimate_dy+=[wue]
                                    nue_w_dmclimate_dy+=[nue]
                                    A_w_dmclimate_dy+=[A]

                            if iiii==9:                                
                                if iii==6:
                                    wue_w_wmclimate_dy+=[wue]
                                    nue_w_wmclimate_dy+=[nue]
                                    A_w_wmclimate_dy+=[A]   
                                    
                            if iiii==2:                                
                                if iii==1:
                                    wue_w_inter_vartemp+=[wue]
                                    nue_w_inter_vartemp+=[nue]
                                    A_w_inter_vartemp+=[A]
                  
                            if iiii==3:                                
                                if iii==1:
                                    wue_w_inter_varvwc+=[wue]
                                    nue_w_inter_varvwc+=[nue]
                                    A_w_inter_varvwc+=[A]

                            if iiii==7:                                
                                if iii==3:
                                    wue_w_inter_varna+=[wue]
                                    nue_w_inter_varna+=[nue]
                                    A_w_inter_varna+=[A]         
                                    
                            if iiii==4:                                
                                if iii==3:
                                    wue_w_inter_mmclimate+=[wue]
                                    nue_w_inter_mmclimate+=[nue]
                                    A_w_inter_mmclimate+=[A]
                                    tempss_wm+=[tl]
 
                            if iiii==5:                                
                                if iii==4:
                                    wue_w_inter_dmclimate+=[wue]
                                    nue_w_inter_dmclimate+=[nue]
                                    A_w_inter_dmclimate+=[A]

                            if iiii==6:                                
                                if iii==5:
                                    wue_w_inter_wmclimate+=[wue]
                                    nue_w_inter_wmclimate+=[nue]
                                    A_w_inter_wmclimate+=[A]
                  
                            if iiii==13:                                
                                if iii==7:
                                    wue_w_dy+=[wue]
                                    nue_w_dy+=[nue]
                                    A_w_dy+=[A]        

        na_lst+=[np.mean(na_lst0)]
        
#        if time>140 and time<gsf+1:
#            print wue_d_inter
        for iii in range(3):
            
            if iii==0:
                nue_d_inter_stats[iii]+=[np.mean(nue_d_inter)]
                wue_d_inter_stats[iii]+=[np.mean(wue_d_inter)]
                A_d_inter_stats[iii]+=[np.mean(A_d_inter)]
                
                
                nue_w_inter_stats[iii]+=[np.mean(nue_w_inter)]
                wue_w_inter_stats[iii]+=[np.mean(wue_w_inter)]
                A_w_inter_stats[iii]+=[np.mean(A_w_inter)]
                
                
                wue_d_dmclimate_gsl_stats[iii]+=[np.mean(wue_d_dmclimate_gsl)]
                nue_d_dmclimate_gsl_stats[iii]+=[np.mean(nue_d_dmclimate_gsl)]      
                A_d_dmclimate_gsl_stats[iii]+=[np.mean(A_d_dmclimate_gsl)]       
            
                wue_d_wmclimate_gsl_stats[iii]+=[np.mean(wue_d_wmclimate_gsl)]
                nue_d_wmclimate_gsl_stats[iii]+=[np.mean(nue_d_wmclimate_gsl)]     
                A_d_wmclimate_gsl_stats[iii]+=[np.mean(A_d_wmclimate_gsl)]
            
            
                wue_d_dmclimate_dy_stats[iii]+=[np.mean(wue_d_dmclimate_dy)]
                nue_d_dmclimate_dy_stats[iii]+=[np.mean(nue_d_dmclimate_dy)]        
                A_d_dmclimate_dy_stats[iii]+=[np.mean(A_d_dmclimate_dy)]        
            
                wue_d_wmclimate_dy_stats[iii]+=[np.mean(wue_d_wmclimate_dy)]
                nue_d_wmclimate_dy_stats[iii]+=[np.mean(nue_d_wmclimate_dy)]      
                A_d_wmclimate_dy_stats[iii]+=[np.mean(A_d_wmclimate_dy)]
                
                wue_d_dmclimate_temps_stats[iii]+=[np.mean(wue_d_dmclimate_temps)]
                nue_d_dmclimate_temps_stats[iii]+=[np.mean(nue_d_dmclimate_temps)]     
                A_d_dmclimate_temps_stats[iii]+=[np.mean(A_d_dmclimate_temps)]    
            
                wue_d_wmclimate_temps_stats[iii]+=[np.mean(wue_d_wmclimate_temps)]
                nue_d_wmclimate_temps_stats[iii]+=[np.mean(nue_d_wmclimate_temps)]       
                A_d_wmclimate_temps_stats[iii]+=[np.mean(A_d_wmclimate_temps)] 
            
            
                wue_w_dmclimate_gsl_stats[iii]+=[np.mean(wue_w_dmclimate_gsl)]
                nue_w_dmclimate_gsl_stats[iii]+=[np.mean(nue_w_dmclimate_gsl)]       
                A_w_dmclimate_gsl_stats[iii]+=[np.mean(A_w_dmclimate_gsl)]     
            
                wue_w_wmclimate_gsl_stats[iii]+=[np.mean(wue_w_wmclimate_gsl)]
                nue_w_wmclimate_gsl_stats[iii]+=[np.mean(nue_w_wmclimate_gsl)]     
                A_w_wmclimate_gsl_stats[iii]+=[np.mean(A_w_wmclimate_gsl)]
            
            
                wue_w_dmclimate_dy_stats[iii]+=[np.mean(wue_w_dmclimate_dy)]
                nue_w_dmclimate_dy_stats[iii]+=[np.mean(nue_w_dmclimate_dy)]        
                A_w_dmclimate_dy_stats[iii]+=[np.mean(A_w_dmclimate_dy)]      
            
                wue_w_wmclimate_dy_stats[iii]+=[np.mean(wue_w_wmclimate_dy)]
                nue_w_wmclimate_dy_stats[iii]+=[np.mean(nue_w_wmclimate_dy)]      
                A_w_wmclimate_dy_stats[iii]+=[np.mean(A_w_wmclimate_dy)]
                
                wue_w_dmclimate_temps_stats[iii]+=[np.mean(wue_w_dmclimate_temps)]
                nue_w_dmclimate_temps_stats[iii]+=[np.mean(nue_w_dmclimate_temps)]      
                A_w_dmclimate_temps_stats[iii]+=[np.mean(A_w_dmclimate_temps)]    
            
                wue_w_wmclimate_temps_stats[iii]+=[np.mean(wue_w_wmclimate_temps)]
                nue_w_wmclimate_temps_stats[iii]+=[np.mean(nue_w_wmclimate_temps)]   
                A_w_wmclimate_temps_stats[iii]+=[np.mean(A_w_wmclimate_temps)]


                
                
                nue_d_inter_vartemp_stats[iii]+=[np.mean(nue_d_inter_vartemp)]
                wue_d_inter_vartemp_stats[iii]+=[np.mean(wue_d_inter_vartemp)]
                A_d_inter_vartemp_stats[iii]+=[np.mean(A_d_inter_vartemp)]
                nue_w_inter_vartemp_stats[iii]+=[np.mean(nue_w_inter_vartemp)]
                wue_w_inter_vartemp_stats[iii]+=[np.mean(wue_w_inter_vartemp)]
                A_w_inter_vartemp_stats[iii]+=[np.mean(A_w_inter_vartemp)]

                nue_d_inter_varvwc_stats[iii]+=[np.mean(nue_d_inter_varvwc)]
                wue_d_inter_varvwc_stats[iii]+=[np.mean(wue_d_inter_varvwc)]
                A_d_inter_varvwc_stats[iii]+=[np.mean(A_d_inter_varvwc)]
                nue_w_inter_varvwc_stats[iii]+=[np.mean(nue_w_inter_varvwc)]
                wue_w_inter_varvwc_stats[iii]+=[np.mean(wue_w_inter_varvwc)]
                A_w_inter_varvwc_stats[iii]+=[np.mean(A_w_inter_varvwc)]

                nue_d_inter_varna_stats[iii]+=[np.mean(nue_d_inter_varna)]
                wue_d_inter_varna_stats[iii]+=[np.mean(wue_d_inter_varna)]
                A_d_inter_varna_stats[iii]+=[np.mean(A_d_inter_varna)]
                nue_w_inter_varna_stats[iii]+=[np.mean(nue_w_inter_varna)]
                wue_w_inter_varna_stats[iii]+=[np.mean(wue_w_inter_varna)]
                A_w_inter_varna_stats[iii]+=[np.mean(A_w_inter_varna)]

                wue_w_inter_mmclimate_stats[iii]+=[np.mean(wue_w_inter_mmclimate)]
                nue_w_inter_mmclimate_stats[iii]+=[np.mean(nue_w_inter_mmclimate)]
                A_w_inter_mmclimate_stats[iii]+=[np.mean(A_w_inter_mmclimate)]
                wue_d_inter_mmclimate_stats[iii]+=[np.mean(wue_d_inter_mmclimate)]
                nue_d_inter_mmclimate_stats[iii]+=[np.mean(nue_d_inter_mmclimate)]
                A_d_inter_mmclimate_stats[iii]+=[np.mean(A_d_inter_mmclimate)]                


                wue_w_inter_dmclimate_stats[iii]+=[np.mean(wue_w_inter_dmclimate)]
                nue_w_inter_dmclimate_stats[iii]+=[np.mean(nue_w_inter_dmclimate)]
                A_w_inter_dmclimate_stats[iii]+=[np.mean(A_w_inter_dmclimate)]
                wue_d_inter_dmclimate_stats[iii]+=[np.mean(wue_d_inter_dmclimate)]
                nue_d_inter_dmclimate_stats[iii]+=[np.mean(nue_d_inter_dmclimate)]
                A_d_inter_dmclimate_stats[iii]+=[np.mean(A_d_inter_dmclimate)]    

                wue_w_inter_wmclimate_stats[iii]+=[np.mean(wue_w_inter_wmclimate)]
                nue_w_inter_wmclimate_stats[iii]+=[np.mean(nue_w_inter_wmclimate)]
                A_w_inter_wmclimate_stats[iii]+=[np.mean(A_w_inter_wmclimate)]
                wue_d_inter_wmclimate_stats[iii]+=[np.mean(wue_d_inter_wmclimate)]
                nue_d_inter_wmclimate_stats[iii]+=[np.mean(nue_d_inter_wmclimate)]
                A_d_inter_wmclimate_stats[iii]+=[np.mean(A_d_inter_wmclimate)]    

                
                
                wue_d_dy_stats[iii]+=[np.mean(wue_d_dy)]
                nue_d_dy_stats[iii]+=[np.mean(nue_d_dy)]        
                A_d_dy_stats[iii]+=[np.mean(A_d_dy)]                 
                
                wue_w_dy_stats[iii]+=[np.mean(wue_w_dy)]
                nue_w_dy_stats[iii]+=[np.mean(nue_w_dy)]        
                A_w_dy_stats[iii]+=[np.mean(A_w_dy)]          
         
                
                
            if iii==1:
                nue_d_inter_stats[iii]+=[np.min(nue_d_inter)]
                wue_d_inter_stats[iii]+=[np.min(wue_d_inter)]
                A_d_inter_stats[iii]+=[np.min(A_d_inter)]

                nue_w_inter_stats[iii]+=[np.min(nue_w_inter)]
                wue_w_inter_stats[iii]+=[np.min(wue_w_inter)]
                A_w_inter_stats[iii]+=[np.min(A_w_inter)]

                
                wue_d_dmclimate_gsl_stats[iii]+=[np.min(wue_d_dmclimate_gsl)]
                nue_d_dmclimate_gsl_stats[iii]+=[np.min(nue_d_dmclimate_gsl)]      
                A_d_dmclimate_gsl_stats[iii]+=[np.min(A_d_dmclimate_gsl)]       
            
                wue_d_wmclimate_gsl_stats[iii]+=[np.min(wue_d_wmclimate_gsl)]
                nue_d_wmclimate_gsl_stats[iii]+=[np.min(nue_d_wmclimate_gsl)]     
                A_d_wmclimate_gsl_stats[iii]+=[np.min(A_d_wmclimate_gsl)]
            
            
                wue_d_dmclimate_dy_stats[iii]+=[np.min(wue_d_dmclimate_dy)]
                nue_d_dmclimate_dy_stats[iii]+=[np.min(nue_d_dmclimate_dy)]        
                A_d_dmclimate_dy_stats[iii]+=[np.min(A_d_dmclimate_dy)]        
            
                wue_d_wmclimate_dy_stats[iii]+=[np.min(wue_d_wmclimate_dy)]
                nue_d_wmclimate_dy_stats[iii]+=[np.min(nue_d_wmclimate_dy)]      
                A_d_wmclimate_dy_stats[iii]+=[np.min(A_d_wmclimate_dy)]
                
                wue_d_dmclimate_temps_stats[iii]+=[np.min(wue_d_dmclimate_temps)]
                nue_d_dmclimate_temps_stats[iii]+=[np.min(nue_d_dmclimate_temps)]     
                A_d_dmclimate_temps_stats[iii]+=[np.min(A_d_dmclimate_temps)]    
            
                wue_d_wmclimate_temps_stats[iii]+=[np.min(wue_d_wmclimate_temps)]
                nue_d_wmclimate_temps_stats[iii]+=[np.min(nue_d_wmclimate_temps)]       
                A_d_wmclimate_temps_stats[iii]+=[np.min(A_d_wmclimate_temps)] 
            
            
                wue_w_dmclimate_gsl_stats[iii]+=[np.min(wue_w_dmclimate_gsl)]
                nue_w_dmclimate_gsl_stats[iii]+=[np.min(nue_w_dmclimate_gsl)]       
                A_w_dmclimate_gsl_stats[iii]+=[np.min(A_w_dmclimate_gsl)]     
            
                wue_w_wmclimate_gsl_stats[iii]+=[np.min(wue_w_wmclimate_gsl)]
                nue_w_wmclimate_gsl_stats[iii]+=[np.min(nue_w_wmclimate_gsl)]     
                A_w_wmclimate_gsl_stats[iii]+=[np.min(A_w_wmclimate_gsl)]
            
            
                wue_w_dmclimate_dy_stats[iii]+=[np.min(wue_w_dmclimate_dy)]
                nue_w_dmclimate_dy_stats[iii]+=[np.min(nue_w_dmclimate_dy)]        
                A_w_dmclimate_dy_stats[iii]+=[np.min(A_w_dmclimate_dy)]      
            
                wue_w_wmclimate_dy_stats[iii]+=[np.min(wue_w_wmclimate_dy)]
                nue_w_wmclimate_dy_stats[iii]+=[np.min(nue_w_wmclimate_dy)]      
                A_w_wmclimate_dy_stats[iii]+=[np.min(A_w_wmclimate_dy)]
                
                wue_w_dmclimate_temps_stats[iii]+=[np.min(wue_w_dmclimate_temps)]
                nue_w_dmclimate_temps_stats[iii]+=[np.min(nue_w_dmclimate_temps)]      
                A_w_dmclimate_temps_stats[iii]+=[np.min(A_w_dmclimate_temps)]    
            
                wue_w_wmclimate_temps_stats[iii]+=[np.min(wue_w_wmclimate_temps)]
                nue_w_wmclimate_temps_stats[iii]+=[np.min(nue_w_wmclimate_temps)]   
                A_w_wmclimate_temps_stats[iii]+=[np.min(A_w_wmclimate_temps)]
                
                
                

                nue_d_inter_vartemp_stats[iii]+=[np.min(nue_d_inter_vartemp)]
                wue_d_inter_vartemp_stats[iii]+=[np.min(wue_d_inter_vartemp)]
                A_d_inter_vartemp_stats[iii]+=[np.min(A_d_inter_vartemp)]
                nue_w_inter_vartemp_stats[iii]+=[np.min(nue_w_inter_vartemp)]
                wue_w_inter_vartemp_stats[iii]+=[np.min(wue_w_inter_vartemp)]
                A_w_inter_vartemp_stats[iii]+=[np.min(A_w_inter_vartemp)]

                nue_d_inter_varvwc_stats[iii]+=[np.min(nue_d_inter_varvwc)]
                wue_d_inter_varvwc_stats[iii]+=[np.min(wue_d_inter_varvwc)]
                A_d_inter_varvwc_stats[iii]+=[np.min(A_d_inter_varvwc)]
                nue_w_inter_varvwc_stats[iii]+=[np.min(nue_w_inter_varvwc)]
                wue_w_inter_varvwc_stats[iii]+=[np.min(wue_w_inter_varvwc)]
                A_w_inter_varvwc_stats[iii]+=[np.min(A_w_inter_varvwc)]

                nue_d_inter_varna_stats[iii]+=[np.min(nue_d_inter_varna)]
                wue_d_inter_varna_stats[iii]+=[np.min(wue_d_inter_varna)]
                A_d_inter_varna_stats[iii]+=[np.min(A_d_inter_varna)]
                nue_w_inter_varna_stats[iii]+=[np.min(nue_w_inter_varna)]
                wue_w_inter_varna_stats[iii]+=[np.min(wue_w_inter_varna)]
                A_w_inter_varna_stats[iii]+=[np.min(A_w_inter_varna)]

                wue_w_inter_mmclimate_stats[iii]+=[np.min(wue_w_inter_mmclimate)]
                nue_w_inter_mmclimate_stats[iii]+=[np.min(nue_w_inter_mmclimate)]
                A_w_inter_mmclimate_stats[iii]+=[np.min(A_w_inter_mmclimate)]
                wue_d_inter_mmclimate_stats[iii]+=[np.min(wue_d_inter_mmclimate)]
                nue_d_inter_mmclimate_stats[iii]+=[np.min(nue_d_inter_mmclimate)]
                A_d_inter_mmclimate_stats[iii]+=[np.min(A_d_inter_mmclimate)]     

                wue_w_inter_dmclimate_stats[iii]+=[np.min(wue_w_inter_dmclimate)]
                nue_w_inter_dmclimate_stats[iii]+=[np.min(nue_w_inter_dmclimate)]
                A_w_inter_dmclimate_stats[iii]+=[np.min(A_w_inter_dmclimate)]
                wue_d_inter_dmclimate_stats[iii]+=[np.min(wue_d_inter_dmclimate)]
                nue_d_inter_dmclimate_stats[iii]+=[np.min(nue_d_inter_dmclimate)]
                A_d_inter_dmclimate_stats[iii]+=[np.min(A_d_inter_dmclimate)]  

                wue_w_inter_wmclimate_stats[iii]+=[np.min(wue_w_inter_wmclimate)]
                nue_w_inter_wmclimate_stats[iii]+=[np.min(nue_w_inter_wmclimate)]
                A_w_inter_wmclimate_stats[iii]+=[np.min(A_w_inter_wmclimate)]
                wue_d_inter_wmclimate_stats[iii]+=[np.min(wue_d_inter_wmclimate)]
                nue_d_inter_wmclimate_stats[iii]+=[np.min(nue_d_inter_wmclimate)]
                A_d_inter_wmclimate_stats[iii]+=[np.min(A_d_inter_wmclimate)]  



                wue_d_dy_stats[iii]+=[np.min(wue_d_dy)]
                nue_d_dy_stats[iii]+=[np.min(nue_d_dy)]        
                A_d_dy_stats[iii]+=[np.min(A_d_dy)]                 
                
                wue_w_dy_stats[iii]+=[np.min(wue_w_dy)]
                nue_w_dy_stats[iii]+=[np.min(nue_w_dy)]        
                A_w_dy_stats[iii]+=[np.min(A_w_dy)]          
         

            if iii==2:
                nue_d_inter_stats[iii]+=[np.max(nue_d_inter)]
                wue_d_inter_stats[iii]+=[np.max(wue_d_inter)]
                A_d_inter_stats[iii]+=[np.max(A_d_inter)]
    
                
                nue_w_inter_stats[iii]+=[np.max(nue_w_inter)]
                wue_w_inter_stats[iii]+=[np.max(wue_w_inter)]
                A_w_inter_stats[iii]+=[np.max(A_w_inter)]        
                
                

                wue_d_dmclimate_gsl_stats[iii]+=[np.max(wue_d_dmclimate_gsl)]
                nue_d_dmclimate_gsl_stats[iii]+=[np.max(nue_d_dmclimate_gsl)]      
                A_d_dmclimate_gsl_stats[iii]+=[np.max(A_d_dmclimate_gsl)]       
            
                wue_d_wmclimate_gsl_stats[iii]+=[np.max(wue_d_wmclimate_gsl)]
                nue_d_wmclimate_gsl_stats[iii]+=[np.max(nue_d_wmclimate_gsl)]     
                A_d_wmclimate_gsl_stats[iii]+=[np.max(A_d_wmclimate_gsl)]
            
            
                wue_d_dmclimate_dy_stats[iii]+=[np.max(wue_d_dmclimate_dy)]
                nue_d_dmclimate_dy_stats[iii]+=[np.max(nue_d_dmclimate_dy)]        
                A_d_dmclimate_dy_stats[iii]+=[np.max(A_d_dmclimate_dy)]        
            
                wue_d_wmclimate_dy_stats[iii]+=[np.max(wue_d_wmclimate_dy)]
                nue_d_wmclimate_dy_stats[iii]+=[np.max(nue_d_wmclimate_dy)]      
                A_d_wmclimate_dy_stats[iii]+=[np.max(A_d_wmclimate_dy)]
                
                wue_d_dmclimate_temps_stats[iii]+=[np.max(wue_d_dmclimate_temps)]
                nue_d_dmclimate_temps_stats[iii]+=[np.max(nue_d_dmclimate_temps)]     
                A_d_dmclimate_temps_stats[iii]+=[np.max(A_d_dmclimate_temps)]    
            
                wue_d_wmclimate_temps_stats[iii]+=[np.max(wue_d_wmclimate_temps)]
                nue_d_wmclimate_temps_stats[iii]+=[np.max(nue_d_wmclimate_temps)]       
                A_d_wmclimate_temps_stats[iii]+=[np.max(A_d_wmclimate_temps)] 
            
            
                wue_w_dmclimate_gsl_stats[iii]+=[np.max(wue_w_dmclimate_gsl)]
                nue_w_dmclimate_gsl_stats[iii]+=[np.max(nue_w_dmclimate_gsl)]       
                A_w_dmclimate_gsl_stats[iii]+=[np.max(A_w_dmclimate_gsl)]     
            
                wue_w_wmclimate_gsl_stats[iii]+=[np.max(wue_w_wmclimate_gsl)]
                nue_w_wmclimate_gsl_stats[iii]+=[np.max(nue_w_wmclimate_gsl)]     
                A_w_wmclimate_gsl_stats[iii]+=[np.max(A_w_wmclimate_gsl)]
            
            
                wue_w_dmclimate_dy_stats[iii]+=[np.max(wue_w_dmclimate_dy)]
                nue_w_dmclimate_dy_stats[iii]+=[np.max(nue_w_dmclimate_dy)]        
                A_w_dmclimate_dy_stats[iii]+=[np.max(A_w_dmclimate_dy)]      
            
                wue_w_wmclimate_dy_stats[iii]+=[np.max(wue_w_wmclimate_dy)]
                nue_w_wmclimate_dy_stats[iii]+=[np.max(nue_w_wmclimate_dy)]      
                A_w_wmclimate_dy_stats[iii]+=[np.max(A_w_wmclimate_dy)]
                
                wue_w_dmclimate_temps_stats[iii]+=[np.max(wue_w_dmclimate_temps)]
                nue_w_dmclimate_temps_stats[iii]+=[np.max(nue_w_dmclimate_temps)]      
                A_w_dmclimate_temps_stats[iii]+=[np.max(A_w_dmclimate_temps)]    
            
                wue_w_wmclimate_temps_stats[iii]+=[np.max(wue_w_wmclimate_temps)]
                nue_w_wmclimate_temps_stats[iii]+=[np.max(nue_w_wmclimate_temps)]   
                A_w_wmclimate_temps_stats[iii]+=[np.max(A_w_wmclimate_temps)]
                
                
              

                nue_d_inter_vartemp_stats[iii]+=[np.max(nue_d_inter_vartemp)]
                wue_d_inter_vartemp_stats[iii]+=[np.max(wue_d_inter_vartemp)]
                A_d_inter_vartemp_stats[iii]+=[np.max(A_d_inter_vartemp)]
                nue_w_inter_vartemp_stats[iii]+=[np.max(nue_w_inter_vartemp)]
                wue_w_inter_vartemp_stats[iii]+=[np.max(wue_w_inter_vartemp)]
                A_w_inter_vartemp_stats[iii]+=[np.max(A_w_inter_vartemp)]

                nue_d_inter_varvwc_stats[iii]+=[np.max(nue_d_inter_varvwc)]
                wue_d_inter_varvwc_stats[iii]+=[np.max(wue_d_inter_varvwc)]
                A_d_inter_varvwc_stats[iii]+=[np.max(A_d_inter_varvwc)]
                nue_w_inter_varvwc_stats[iii]+=[np.max(nue_w_inter_varvwc)]
                wue_w_inter_varvwc_stats[iii]+=[np.max(wue_w_inter_varvwc)]
                A_w_inter_varvwc_stats[iii]+=[np.max(A_w_inter_varvwc)]        

                nue_d_inter_varna_stats[iii]+=[np.max(nue_d_inter_varna)]
                wue_d_inter_varna_stats[iii]+=[np.max(wue_d_inter_varna)]
                A_d_inter_varna_stats[iii]+=[np.max(A_d_inter_varna)]
                nue_w_inter_varna_stats[iii]+=[np.max(nue_w_inter_varna)]
                wue_w_inter_varna_stats[iii]+=[np.max(wue_w_inter_varna)]
                A_w_inter_varna_stats[iii]+=[np.max(A_w_inter_varna)] 

                wue_w_inter_mmclimate_stats[iii]+=[np.max(wue_w_inter_mmclimate)]
                nue_w_inter_mmclimate_stats[iii]+=[np.max(nue_w_inter_mmclimate)]
                A_w_inter_mmclimate_stats[iii]+=[np.max(A_w_inter_mmclimate)]
                wue_d_inter_mmclimate_stats[iii]+=[np.max(wue_d_inter_mmclimate)]
                nue_d_inter_mmclimate_stats[iii]+=[np.max(nue_d_inter_mmclimate)]
                A_d_inter_mmclimate_stats[iii]+=[np.max(A_d_inter_mmclimate)]     


                wue_w_inter_dmclimate_stats[iii]+=[np.max(wue_w_inter_dmclimate)]
                nue_w_inter_dmclimate_stats[iii]+=[np.max(nue_w_inter_dmclimate)]
                A_w_inter_dmclimate_stats[iii]+=[np.max(A_w_inter_dmclimate)]
                wue_d_inter_dmclimate_stats[iii]+=[np.max(wue_d_inter_dmclimate)]
                nue_d_inter_dmclimate_stats[iii]+=[np.max(nue_d_inter_dmclimate)]
                A_d_inter_dmclimate_stats[iii]+=[np.max(A_d_inter_dmclimate)]  

                wue_w_inter_wmclimate_stats[iii]+=[np.max(wue_w_inter_wmclimate)]
                nue_w_inter_wmclimate_stats[iii]+=[np.max(nue_w_inter_wmclimate)]
                A_w_inter_wmclimate_stats[iii]+=[np.max(A_w_inter_wmclimate)]
                wue_d_inter_wmclimate_stats[iii]+=[np.max(wue_d_inter_wmclimate)]
                nue_d_inter_wmclimate_stats[iii]+=[np.max(nue_d_inter_wmclimate)]
                A_d_inter_wmclimate_stats[iii]+=[np.max(A_d_inter_wmclimate)]  


                wue_d_dy_stats[iii]+=[np.max(wue_d_dy)]
                nue_d_dy_stats[iii]+=[np.max(nue_d_dy)]        
                A_d_dy_stats[iii]+=[np.max(A_d_dy)]                 
                
                wue_w_dy_stats[iii]+=[np.max(wue_w_dy)]
                nue_w_dy_stats[iii]+=[np.max(nue_w_dy)]        
                A_w_dy_stats[iii]+=[np.max(A_w_dy)]    
               
        tot_temps_dm+=[np.mean(tempss_dm)]
        tot_temps_wm+=[np.mean(tempss_wm)]    


#---------------Smooth Model Output---------------#      
    

    for iii in range(3):
        
        nue_d_inter_stats_smooth+=[gs_smooth(nue_d_inter_stats[iii],gs0_d,gsf)]
        wue_d_inter_stats_smooth+=[gs_smooth(wue_d_inter_stats[iii],gs0_d,gsf)]
        A_d_inter_stats_smooth+=[gs_smooth(A_d_inter_stats[iii],gs0_d,gsf)]


        nue_w_inter_stats_smooth+=[gs_smooth(nue_w_inter_stats[iii],gs0_w,gsf)]
        wue_w_inter_stats_smooth+=[gs_smooth(wue_w_inter_stats[iii],gs0_w,gsf)]
        A_w_inter_stats_smooth+=[gs_smooth(A_w_inter_stats[iii],gs0_w,gsf)]
        
        
        
        wue_d_dmclimate_gsl_stats_smooth+=[gs_smooth(wue_d_dmclimate_gsl_stats[iii],gs0_ddy,gsf)]
        nue_d_dmclimate_gsl_stats_smooth+=[gs_smooth(nue_d_dmclimate_gsl_stats[iii],gs0_ddy,gsf)]   
        A_d_dmclimate_gsl_stats_smooth+=[gs_smooth( A_d_dmclimate_gsl_stats[iii],gs0_ddy,gsf)]     
    
        wue_d_wmclimate_gsl_stats_smooth+=[gs_smooth(wue_d_wmclimate_gsl_stats[iii],gs0_wdy,gsf)]
        nue_d_wmclimate_gsl_stats_smooth+=[gs_smooth(nue_d_wmclimate_gsl_stats[iii],gs0_wdy,gsf)]    
        A_d_wmclimate_gsl_stats_smooth+=[gs_smooth(A_d_wmclimate_gsl_stats[iii],gs0_wdy,gsf)]
    
    
        wue_d_dmclimate_dy_stats_smooth+=[gs_smooth(wue_d_dmclimate_dy_stats[iii],gs0_ddy,gsf)]
        nue_d_dmclimate_dy_stats_smooth+=[gs_smooth(nue_d_dmclimate_dy_stats[iii],gs0_ddy,gsf)]       
        A_d_dmclimate_dy_stats_smooth+=[gs_smooth(A_d_dmclimate_dy_stats[iii],gs0_ddy,gsf)]        

        if iii==0:
            for i in range(len(wue_d_dmclimate_dy_stats_smooth[0])):
                if wue_d_dmclimate_dy_stats_smooth[0][i]<0.0:
                    wue_d_dmclimate_dy_stats_smooth[0][i]=0.0
        
            for i in range(len(nue_d_dmclimate_dy_stats_smooth[0])):
                if nue_d_dmclimate_dy_stats_smooth[0][i]<0.0:
                    nue_d_dmclimate_dy_stats_smooth[0][i]=0.0        
            
            for i in range(len(A_d_dmclimate_dy_stats_smooth[0])):
                if A_d_dmclimate_dy_stats_smooth[0][i]<0.0:
                    A_d_dmclimate_dy_stats_smooth[0][i]=0.0    



        wue_d_dy_stats_smooth+=[gs_smooth(wue_d_dy_stats[iii],gs0_ddy,gsf)]
        nue_d_dy_stats_smooth+=[gs_smooth(nue_d_dy_stats[iii],gs0_ddy,gsf)]       
        A_d_dy_stats_smooth+=[gs_smooth(A_d_dy_stats[iii],gs0_ddy,gsf)]        

        if iii==0:
            for i in range(len(wue_d_dy_stats_smooth[0])):
                if wue_d_dy_stats_smooth[0][i]<0.0:
                    wue_d_dy_stats_smooth[0][i]=0.0
        
            for i in range(len(nue_d_dy_stats_smooth[0])):
                if nue_d_dy_stats_smooth[0][i]<0.0:
                    nue_d_dy_stats_smooth[0][i]=0.0        
            
            for i in range(len(A_d_dy_stats_smooth[0])):
                if A_d_dy_stats_smooth[0][i]<0.0:
                    A_d_dy_stats_smooth[0][i]=0.0    
    
    
        wue_d_wmclimate_dy_stats_smooth+=[gs_smooth(wue_d_wmclimate_dy_stats[iii],gs0_wdy,gsf)]
        nue_d_wmclimate_dy_stats_smooth+=[gs_smooth(nue_d_wmclimate_dy_stats[iii],gs0_wdy,gsf)]     
        A_d_wmclimate_dy_stats_smooth+=[gs_smooth(A_d_wmclimate_dy_stats[iii],gs0_wdy,gsf)]
        
        wue_d_dmclimate_temps_stats_smooth+=[gs_smooth(wue_d_dmclimate_temps_stats[iii],gs0_d,gsf)]
        nue_d_dmclimate_temps_stats_smooth+=[gs_smooth(nue_d_dmclimate_temps_stats[iii],gs0_d,gsf)]  
        A_d_dmclimate_temps_stats_smooth+=[gs_smooth(A_d_dmclimate_temps_stats[iii],gs0_d,gsf)]    
    
        wue_d_wmclimate_temps_stats_smooth+=[gs_smooth(wue_d_wmclimate_temps_stats[iii],gs0_w,gsf)]
        nue_d_wmclimate_temps_stats_smooth+=[gs_smooth(nue_d_wmclimate_temps_stats[iii],gs0_w,gsf)]       
        A_d_wmclimate_temps_stats_smooth+=[gs_smooth(A_d_wmclimate_temps_stats[iii],gs0_w,gsf)]
    
    
        wue_w_dmclimate_gsl_stats_smooth+=[gs_smooth(wue_w_dmclimate_gsl_stats[iii],gs0_ddy,gsf)]
        nue_w_dmclimate_gsl_stats_smooth+=[gs_smooth(nue_w_dmclimate_gsl_stats[iii],gs0_ddy,gsf)]   
        A_w_dmclimate_gsl_stats_smooth+=[gs_smooth(A_w_dmclimate_gsl_stats[iii],gs0_ddy,gsf)]    
    
        wue_w_wmclimate_gsl_stats_smooth+=[gs_smooth(wue_w_wmclimate_gsl_stats[iii],gs0_wdy,gsf)]
        nue_w_wmclimate_gsl_stats_smooth+=[gs_smooth(nue_w_wmclimate_gsl_stats[iii],gs0_wdy,gsf)]    
        A_w_wmclimate_gsl_stats_smooth+=[gs_smooth(A_w_wmclimate_gsl_stats[iii],gs0_wdy,gsf)]
    
    
        wue_w_dmclimate_dy_stats_smooth+=[gs_smooth(wue_w_dmclimate_dy_stats[iii],gs0_ddy,gsf)]
        nue_w_dmclimate_dy_stats_smooth+=[gs_smooth(nue_w_dmclimate_dy_stats[iii],gs0_ddy,gsf)]     
        A_w_dmclimate_dy_stats_smooth+=[gs_smooth(A_w_dmclimate_dy_stats[iii],gs0_ddy,gsf)]  
        
        if iii==0:
            for i in range(len(wue_w_dmclimate_dy_stats_smooth[0])):
                if wue_w_dmclimate_dy_stats_smooth[0][i]<0.0:
                    wue_w_dmclimate_dy_stats_smooth[0][i]=0.0
        
            for i in range(len(nue_w_dmclimate_dy_stats_smooth[0])):
                if nue_w_dmclimate_dy_stats_smooth[0][i]<0.0:
                    nue_w_dmclimate_dy_stats_smooth[0][i]=0.0        
            
            for i in range(len(A_w_dmclimate_dy_stats_smooth[0])):
                if A_w_dmclimate_dy_stats_smooth[0][i]<0.0:
                    A_w_dmclimate_dy_stats_smooth[0][i]=0.0         
                    
                    
        wue_w_dy_stats_smooth+=[gs_smooth(wue_w_dy_stats[iii],gs0_wdy,gsf)]
        nue_w_dy_stats_smooth+=[gs_smooth(nue_w_dy_stats[iii],gs0_wdy,gsf)]       
        A_w_dy_stats_smooth+=[gs_smooth(A_w_dy_stats[iii],gs0_wdy,gsf)]        

        if iii==0:
            for i in range(len(wue_w_dy_stats_smooth[0])):
                if wue_w_dy_stats_smooth[0][i]<0.0:
                    wue_w_dy_stats_smooth[0][i]=0.0
        
            for i in range(len(nue_w_dy_stats_smooth[0])):
                if nue_w_dy_stats_smooth[0][i]<0.0:
                    nue_w_dy_stats_smooth[0][i]=0.0        
            
            for i in range(len(A_w_dy_stats_smooth[0])):
                if A_w_dy_stats_smooth[0][i]<0.0:
                    A_w_dy_stats_smooth[0][i]=0.0  
                    
                    
    
        wue_w_wmclimate_dy_stats_smooth+=[gs_smooth(wue_w_wmclimate_dy_stats[iii],gs0_wdy,gsf)]
        nue_w_wmclimate_dy_stats_smooth+=[gs_smooth(nue_w_wmclimate_dy_stats[iii],gs0_wdy,gsf)]     
        A_w_wmclimate_dy_stats_smooth+=[gs_smooth(A_w_wmclimate_dy_stats[iii],gs0_wdy,gsf)]
        
        wue_w_dmclimate_temps_stats_smooth+=[gs_smooth(wue_w_dmclimate_temps_stats[iii],gs0_d,gsf)]
        nue_w_dmclimate_temps_stats_smooth+=[gs_smooth(nue_w_dmclimate_temps_stats[iii],gs0_d,gsf)]  
        A_w_dmclimate_temps_stats_smooth+=[gs_smooth(A_w_dmclimate_temps_stats[iii],gs0_d,gsf)]   
    
        wue_w_wmclimate_temps_stats_smooth+=[gs_smooth(wue_w_wmclimate_temps_stats[iii],gs0_w,gsf)]
        nue_w_wmclimate_temps_stats_smooth+=[gs_smooth(nue_w_wmclimate_temps_stats[iii],gs0_w,gsf)] 
        A_w_wmclimate_temps_stats_smooth+=[gs_smooth(A_w_wmclimate_temps_stats[iii],gs0_w,gsf)]
                
                

        nue_d_inter_vartemp_stats_smooth+=[gs_smooth(nue_d_inter_vartemp_stats[iii],0,364)]
        wue_d_inter_vartemp_stats_smooth+=[gs_smooth(wue_d_inter_vartemp_stats[iii],0,364)]
        A_d_inter_vartemp_stats_smooth+=[gs_smooth(A_d_inter_vartemp_stats[iii],0,364)]
        nue_w_inter_vartemp_stats_smooth+=[gs_smooth(nue_w_inter_vartemp_stats[iii],0,364)]
        wue_w_inter_vartemp_stats_smooth+=[gs_smooth(wue_w_inter_vartemp_stats[iii],0,364)]
        A_w_inter_vartemp_stats_smooth+=[gs_smooth(A_w_inter_vartemp_stats[iii],0,364)]

        nue_d_inter_varvwc_stats_smooth+=[gs_smooth(nue_d_inter_varvwc_stats[iii],0,364)]
        wue_d_inter_varvwc_stats_smooth+=[gs_smooth(wue_d_inter_varvwc_stats[iii],0,364)]
        A_d_inter_varvwc_stats_smooth+=[gs_smooth(A_d_inter_varvwc_stats[iii],0,364)]
        nue_w_inter_varvwc_stats_smooth+=[gs_smooth(nue_w_inter_varvwc_stats[iii],0,364)]
        wue_w_inter_varvwc_stats_smooth+=[gs_smooth(wue_w_inter_varvwc_stats[iii],0,364)]
        A_w_inter_varvwc_stats_smooth+=[gs_smooth(A_w_inter_varvwc_stats[iii],0,364)]   

        nue_d_inter_varna_stats_smooth+=[gs_smooth(nue_d_inter_varna_stats[iii],gs0_m,gsf)]
        wue_d_inter_varna_stats_smooth+=[gs_smooth(wue_d_inter_varna_stats[iii],gs0_m,gsf)]
        A_d_inter_varna_stats_smooth+=[gs_smooth(A_d_inter_varna_stats[iii],gs0_m,gsf)]
        nue_w_inter_varna_stats_smooth+=[gs_smooth(nue_w_inter_varna_stats[iii],gs0_m,gsf)]
        wue_w_inter_varna_stats_smooth+=[gs_smooth(wue_w_inter_varna_stats[iii],gs0_m,gsf)]
        A_w_inter_varna_stats_smooth+=[gs_smooth(A_w_inter_varna_stats[iii],gs0_m,gsf)]   


        wue_w_inter_mmclimate_stats_smooth+=[gs_smooth(wue_w_inter_mmclimate_stats[iii],gs0_m,gsf)]
        nue_w_inter_mmclimate_stats_smooth+=[gs_smooth(nue_w_inter_mmclimate_stats[iii],gs0_m,gsf)]
        A_w_inter_mmclimate_stats_smooth+=[gs_smooth(A_w_inter_mmclimate_stats[iii],gs0_m,gsf)]
        wue_d_inter_mmclimate_stats_smooth+=[gs_smooth(wue_d_inter_mmclimate_stats[iii],gs0_m,gsf)]
        nue_d_inter_mmclimate_stats_smooth+=[gs_smooth(nue_d_inter_mmclimate_stats[iii],gs0_m,gsf)]
        A_d_inter_mmclimate_stats_smooth+=[gs_smooth(A_d_inter_mmclimate_stats[iii],gs0_m,gsf)]                

        wue_w_inter_dmclimate_stats_smooth+=[gs_smooth(wue_w_inter_dmclimate_stats[iii],gs0_d,gsf)]
        nue_w_inter_dmclimate_stats_smooth+=[gs_smooth(nue_w_inter_dmclimate_stats[iii],gs0_d,gsf)]
        A_w_inter_dmclimate_stats_smooth+=[gs_smooth(A_w_inter_dmclimate_stats[iii],gs0_d,gsf)]
        wue_d_inter_dmclimate_stats_smooth+=[gs_smooth(wue_d_inter_dmclimate_stats[iii],gs0_d,gsf)]
        nue_d_inter_dmclimate_stats_smooth+=[gs_smooth(nue_d_inter_dmclimate_stats[iii],gs0_d,gsf)]
        A_d_inter_dmclimate_stats_smooth+=[gs_smooth(A_d_inter_dmclimate_stats[iii],gs0_d,gsf)]      

        wue_w_inter_wmclimate_stats_smooth+=[gs_smooth(wue_w_inter_wmclimate_stats[iii],gs0_w,gsf)]
        nue_w_inter_wmclimate_stats_smooth+=[gs_smooth(nue_w_inter_wmclimate_stats[iii],gs0_w,gsf)]
        A_w_inter_wmclimate_stats_smooth+=[gs_smooth(A_w_inter_wmclimate_stats[iii],gs0_w,gsf)]
        wue_d_inter_wmclimate_stats_smooth+=[gs_smooth(wue_d_inter_wmclimate_stats[iii],gs0_w,gsf)]
        nue_d_inter_wmclimate_stats_smooth+=[gs_smooth(nue_d_inter_wmclimate_stats[iii],gs0_w,gsf)]
        A_d_inter_wmclimate_stats_smooth+=[gs_smooth(A_d_inter_wmclimate_stats[iii],gs0_w,gsf)]  

###------------------------------------------------------######    
   #dry and wet meadow nue (high temps)

    #extended growing season scenario
#    
#    #dry meadow climate
#    ax1.plot(days, nue_d_dmclimate_dy_stats_smooth[0], 'r-',linewidth=3) 
##    ax3.plot(days, A_d_inter_dy_stats_smooth[2], 'r-',linewidth=3)  
##    ax3.fill_between(days, A_d_inter_stats_smooth[1], A_d_inter_stats_smooth[2],alpha=0.0,color='red')
##    ax3.fill_between(days, A_d_inter_dy_stats_smooth[1], A_d_inter_dy_stats_smooth[2],alpha=0.3,color='red') 
#    ax1.plot(days, nue_w_dmclimate_dy_stats_smooth[0], 'b-',linewidth=3) 
##    ax3.plot(days, A_w_inter_dy_stats_smooth[2], 'b-',linewidth=3)      
#    ax1.plot(days, nue_d_inter_dmclimate_stats_smooth[0], 'r--',linewidth=3) 
##    ax3.plot(days, A_d_inter_stats_smooth[2], 'r--',linewidth=3)    
#    ax1.plot(days, nue_w_inter_dmclimate_stats_smooth[0], 'b--',linewidth=3) 
##    ax3.plot(days, A_w_inter_stats_smooth[2], 'b--',linewidth=3)    
##    ax3.fill_between(days, A_w_inter_stats_smooth[1], A_w_inter_stats_smooth[2],alpha=0.0,color='blue')
##    ax3.fill_between(days, A_w_inter_dy_stats_smooth[1], A_w_inter_dy_stats_smooth[2],alpha=0.3,color='blue')    
##    ax1.legend(loc="upper left",fontsize=20)
#    ax1.annotate('B', xy=get_axis_limits(ax1),fontsize=40)
#
##Wet meadow climate
#    ax31.plot(days, nue_d_wmclimate_dy_stats_smooth[0], 'r-',linewidth=3) 
##    ax3.plot(days, A_d_inter_dy_stats_smooth[2], 'r-',linewidth=3)  
##    ax3.fill_between(days, A_d_inter_stats_smooth[1], A_d_inter_stats_smooth[2],alpha=0.0,color='red')
##    ax3.fill_between(days, A_d_inter_dy_stats_smooth[1], A_d_inter_dy_stats_smooth[2],alpha=0.3,color='red') 
#    ax31.plot(days, nue_w_wmclimate_dy_stats_smooth[0], 'b-',linewidth=3) 
##    ax3.plot(days, A_w_inter_dy_stats_smooth[2], 'b-',linewidth=3)      
#    ax31.plot(days, nue_d_inter_wmclimate_stats_smooth[0], 'r--',linewidth=3) 
##    ax3.plot(days, A_d_inter_stats_smooth[2], 'r--',linewidth=3)    
#    ax31.plot(days, nue_w_inter_wmclimate_stats_smooth[0], 'b--',linewidth=3) 
##    ax3.plot(days, A_w_inter_stats_smooth[2], 'b--',linewidth=3)    
##    ax3.fill_between(days, A_w_inter_stats_smooth[1], A_w_inter_stats_smooth[2],alpha=0.0,color='blue')
##    ax3.fill_between(days, A_w_inter_dy_stats_smooth[1], A_w_inter_dy_stats_smooth[2],alpha=0.3,color='blue')    
##    ax31.legend(loc="upper left",fontsize=20)
#    ax31.annotate('B', xy=get_axis_limits(ax31),fontsize=40)
#    
#    
#    #longer growing season
#    #dry meadow climate
#    ax22.plot(days, nue_d_dmclimate_gsl_stats_smooth[0], 'r-',linewidth=3) 
##    ax24.plot(days, A_d_inter_dy_gsl_stats_smooth[2], 'r-',linewidth=3)  
##    ax24.fill_between(days, A_d_inter_stats_smooth[1], A_d_inter_stats_smooth[2],alpha=0.0,color='red')
##    ax24.fill_between(days, A_d_inter_dy_gsl_stats_smooth[1], A_d_inter_dy_gsl_stats_smooth[2],alpha=0.3,color='red') 
#    ax22.plot(days, nue_w_dmclimate_gsl_stats_smooth[0], 'b-',linewidth=3) 
##    ax24.plot(days, A_w_inter_dy_gsl_stats_smooth[2], 'b-',linewidth=3)      
#    ax22.plot(days, nue_d_inter_dmclimate_stats_smooth[0], 'r--',linewidth=3) 
##    ax24.plot(days, A_d_inter_stats_smooth[2], 'r--',linewidth=3)    
#    ax22.plot(days, nue_w_inter_dmclimate_stats_smooth[0], 'b--',linewidth=3) 
##    ax24.plot(days, A_w_inter_stats_smooth[2], 'b--',linewidth=3)    
##    ax24.fill_between(days, A_w_inter_stats_smooth[1], A_w_inter_stats_smooth[2],alpha=0.0,color='blue')
##    ax24.fill_between(days, A_w_inter_dy_gsl_stats_smooth[1], A_w_inter_dy_gsl_stats_smooth[2],alpha=0.3,color='blue')    
##    ax22.legend(loc="upper left",fontsize=20)
#    ax22.annotate('B', xy=get_axis_limits(ax22),fontsize=40)    
#
##wet meadow climate
#    ax30.plot(days, nue_d_wmclimate_gsl_stats_smooth[0], 'r-',linewidth=3) 
##    ax24.plot(days, A_d_inter_dy_gsl_stats_smooth[2], 'r-',linewidth=3)  
##    ax24.fill_between(days, A_d_inter_stats_smooth[1], A_d_inter_stats_smooth[2],alpha=0.0,color='red')
##    ax24.fill_between(days, A_d_inter_dy_gsl_stats_smooth[1], A_d_inter_dy_gsl_stats_smooth[2],alpha=0.3,color='red') 
#    ax30.plot(days, nue_w_wmclimate_gsl_stats_smooth[0], 'b-',linewidth=3) 
##    ax24.plot(days, A_w_inter_dy_gsl_stats_smooth[2], 'b-',linewidth=3)      
#    ax30.plot(days, nue_d_inter_wmclimate_stats_smooth[0], 'r--',linewidth=3) 
##    ax24.plot(days, A_d_inter_stats_smooth[2], 'r--',linewidth=3)    
#    ax30.plot(days, nue_w_inter_wmclimate_stats_smooth[0], 'b--',linewidth=3) 
##    ax24.plot(days, A_w_inter_stats_smooth[2], 'b--',linewidth=3)    
##    ax24.fill_between(days, A_w_inter_stats_smooth[1], A_w_inter_stats_smooth[2],alpha=0.0,color='blue')
##    ax24.fill_between(days, A_w_inter_dy_gsl_stats_smooth[1], A_w_inter_dy_gsl_stats_smooth[2],alpha=0.3,color='blue')    
##    ax30.legend(loc="upper left",fontsize=20)
#    ax30.annotate('B', xy=get_axis_limits(ax30),fontsize=40) 
#    
#    #hotter temperatures
#    #dry meadow climate
#    ax2.plot(days, nue_d_dmclimate_temps_stats_smooth[0], 'r-',linewidth=3) 
##    ax24.plot(days, A_d_inter_dy_gsl_stats_smooth[2], 'r-',linewidth=3)  
##    ax24.fill_between(days, A_d_inter_stats_smooth[1], A_d_inter_stats_smooth[2],alpha=0.0,color='red')
##    ax24.fill_between(days, A_d_inter_dy_gsl_stats_smooth[1], A_d_inter_dy_gsl_stats_smooth[2],alpha=0.3,color='red') 
#    ax2.plot(days, nue_w_dmclimate_temps_stats_smooth[0], 'b-',linewidth=3) 
##    ax24.plot(days, A_w_inter_dy_gsl_stats_smooth[2], 'b-',linewidth=3)      
#    ax2.plot(days, nue_d_inter_dmclimate_stats_smooth[0], 'r--',linewidth=3) 
##    ax24.plot(days, A_d_inter_stats_smooth[2], 'r--',linewidth=3)    
#    ax2.plot(days, nue_w_inter_dmclimate_stats_smooth[0], 'b--',linewidth=3) 
##    ax24.plot(days, A_w_inter_stats_smooth[2], 'b--',linewidth=3)    
##    ax24.fill_between(days, A_w_inter_stats_smooth[1], A_w_inter_stats_smooth[2],alpha=0.0,color='blue')
##    ax24.fill_between(days, A_w_inter_dy_gsl_stats_smooth[1], A_w_inter_dy_gsl_stats_smooth[2],alpha=0.3,color='blue')    
##    ax2.legend(loc="upper left",fontsize=20)
#    ax2.annotate('B', xy=get_axis_limits(ax2),fontsize=40)
#
##wet meadow climate
#    ax32.plot(days, nue_d_wmclimate_temps_stats_smooth[0], 'r-',linewidth=3) 
##    ax24.plot(days, A_d_inter_dy_gsl_stats_smooth[2], 'r-',linewidth=3)  
##    ax24.fill_between(days, A_d_inter_stats_smooth[1], A_d_inter_stats_smooth[2],alpha=0.0,color='red')
##    ax24.fill_between(days, A_d_inter_dy_gsl_stats_smooth[1], A_d_inter_dy_gsl_stats_smooth[2],alpha=0.3,color='red') 
#    ax32.plot(days, nue_w_wmclimate_temps_stats_smooth[0], 'b-',linewidth=3) 
##    ax24.plot(days, A_w_inter_dy_gsl_stats_smooth[2], 'b-',linewidth=3)      
#    ax32.plot(days, nue_d_inter_wmclimate_stats_smooth[0], 'r--',linewidth=3) 
##    ax24.plot(days, A_d_inter_stats_smooth[2], 'r--',linewidth=3)    
#    ax32.plot(days, nue_w_inter_wmclimate_stats_smooth[0], 'b--',linewidth=3) 
##    ax24.plot(days, A_w_inter_stats_smooth[2], 'b--',linewidth=3)    
##    ax24.fill_between(days, A_w_inter_stats_smooth[1], A_w_inter_stats_smooth[2],alpha=0.0,color='blue')
##    ax24.fill_between(days, A_w_inter_dy_gsl_stats_smooth[1], A_w_inter_dy_gsl_stats_smooth[2],alpha=0.3,color='blue')    
##    ax32.legend(loc="upper left",fontsize=20)
#    ax32.annotate('B', xy=get_axis_limits(ax32),fontsize=40)
 

###------------------------------------------------------######    
    

    #dry and wet meadow A

    #extended growing season scenario

    #conservative and acquisitive species in their respective climates
    l1, =ax36.plot(days, A_d_inter_stats_smooth[0], 'k-',linewidth=3,alpha=0.5,label="Conservative Species\nAverage Growing Season\n") 
    l2, =ax36.plot(days, A_d_dy_stats_smooth[0], 'k--',linewidth=3,alpha=0.5,label="Extended Summer\n$\Delta$ Assimilation=%d%%"%((sum((np.array(A_d_dy_stats_smooth[0])*3600*6)/1000000.*44.)-sum((np.array(A_d_inter_stats_smooth[0])*3600*6)/1000000.*44.))/sum((np.array(A_d_inter_stats_smooth[0])*3600*6)/1000000.*44.)*100.0)) 
    l3, =ax36.plot(days, A_w_inter_stats_smooth[0], 'k-',linewidth=3,alpha=1.0,label="Acquisitive Species\nAverage Growing Season\n") 
    l4, =ax36.plot(days, A_w_dy_stats_smooth[0], 'k--',linewidth=3,alpha=1.0, label="Extended Summer\n$\Delta$ Assimilation=%d%%"%((sum((np.array(A_w_dy_stats_smooth[0])*3600*6)/1000000.*44.)-sum((np.array(A_w_inter_stats_smooth[0])*3600*6)/1000000.*44.))/sum((np.array(A_w_inter_stats_smooth[0])*3600*6)/1000000.*44.)*100.0)) 
    plot_lines=[]
    plot_lines.append([l1,l2])
    plot_lines.append([l3,l4])    
    legend1=ax36.legend(plot_lines[0],["Conservative Species\nAverage Growing Season\n","Extended Summer\n$\Delta$ Assimilation=%d%%"%((sum((np.array(A_d_dy_stats_smooth[0])*3600*6)/1000000.*44.)-sum((np.array(A_d_inter_stats_smooth[0])*3600*6)/1000000.*44.))/sum((np.array(A_d_inter_stats_smooth[0])*3600*6)/1000000.*44.)*100.0)],loc=2,fontsize=15)
    ax36.legend(plot_lines[1],["Acquisitive Species\nAverage Growing Season\n","Extended Summer\n$\Delta$ Assimilation=%d%%"%((sum((np.array(A_w_dy_stats_smooth[0])*3600*6)/1000000.*44.)-sum((np.array(A_w_inter_stats_smooth[0])*3600*6)/1000000.*44.))/sum((np.array(A_w_inter_stats_smooth[0])*3600*6)/1000000.*44.)*100.0)],loc=1,fontsize=15)
    ax36.add_artist(legend1)

    
#    #dry meadow climate
#    ax3.plot(days, A_d_dmclimate_dy_stats_smooth[0], 'r-',linewidth=3,label="Conservative Species Scenario\n$\Delta$ Average Cumulative Assimilation=%d%%"%((sum((np.array(A_d_dmclimate_dy_stats_smooth[0])*3600*6)/1000000.*44.)-sum((np.array(A_d_inter_dmclimate_stats_smooth[0])*3600*6)/1000000.*44.))/sum((np.array(A_d_inter_dmclimate_stats_smooth[0])*3600*6)/1000000.*44.)*100.0)) 
##    ax3.plot(days, A_d_inter_dy_stats_smooth[2], 'r-',linewidth=3)  
##    ax3.fill_between(days, A_d_inter_stats_smooth[1], A_d_inter_stats_smooth[2],alpha=0.0,color='red')
##    ax3.fill_between(days, A_d_inter_dy_stats_smooth[1], A_d_inter_dy_stats_smooth[2],alpha=0.3,color='red') 
#    ax3.plot(days, A_w_dmclimate_dy_stats_smooth[0], 'b-',linewidth=3,label="Acquisitive Species Scenario\n$\Delta$ Average Cumulative Assimilation=%d%%"%((sum((np.array(A_w_dmclimate_dy_stats_smooth[0])*3600*6)/1000000.*44.)-sum((np.array(A_w_inter_dmclimate_stats_smooth[0])*3600*6)/1000000.*44.))/sum((np.array(A_w_inter_dmclimate_stats_smooth[0])*3600*6)/1000000.*44.)*100.0)) 
##    ax3.plot(days, A_w_inter_dy_stats_smooth[2], 'b-',linewidth=3)      
#    ax3.plot(days, A_d_inter_dmclimate_stats_smooth[0], 'r--',linewidth=3,label="Typical Growing Season") 
##    ax3.plot(days, A_d_inter_stats_smooth[2], 'r--',linewidth=3)    
#    ax3.plot(days, A_w_inter_dmclimate_stats_smooth[0], 'b--',linewidth=3,label="Typical Growing Season") 
##    ax3.plot(days, A_w_inter_stats_smooth[2], 'b--',linewidth=3)    
##    ax3.fill_between(days, A_w_inter_stats_smooth[1], A_w_inter_stats_smooth[2],alpha=0.0,color='blue')
##    ax3.fill_between(days, A_w_inter_dy_stats_smooth[1], A_w_inter_dy_stats_smooth[2],alpha=0.3,color='blue')    
#    ax3.legend(loc="upper left",fontsize=20)
#    ax3.annotate('A', xy=get_axis_limits(ax3),fontsize=40)
#
##Wet meadow climate
#    ax33.plot(days, A_d_wmclimate_dy_stats_smooth[0], 'r-',linewidth=3,label="Conservative Species Scenario\n$\Delta$ Average Cumulative Assimilation=%d%%"%((sum((np.array(A_d_wmclimate_dy_stats_smooth[0])*3600*6)/1000000.*44.)-sum((np.array(A_d_inter_wmclimate_stats_smooth[0])*3600*6)/1000000.*44.))/sum((np.array(A_d_inter_wmclimate_stats_smooth[0])*3600*6)/1000000.*44.)*100.0)) 
##    ax3.plot(days, A_d_inter_dy_stats_smooth[2], 'r-',linewidth=3)  
##    ax3.fill_between(days, A_d_inter_stats_smooth[1], A_d_inter_stats_smooth[2],alpha=0.0,color='red')
##    ax3.fill_between(days, A_d_inter_dy_stats_smooth[1], A_d_inter_dy_stats_smooth[2],alpha=0.3,color='red') 
#    ax33.plot(days, A_w_wmclimate_dy_stats_smooth[0], 'b-',linewidth=3,label="Acquisitive Species Scenario\n$\Delta$ Average Cumulative Assimilation=%d%%"%((sum((np.array(A_w_wmclimate_dy_stats_smooth[0])*3600*6)/1000000.*44.)-sum((np.array(A_w_inter_wmclimate_stats_smooth[0])*3600*6)/1000000.*44.))/sum((np.array(A_w_inter_wmclimate_stats_smooth[0])*3600*6)/1000000.*44.)*100.0)) 
##    ax3.plot(days, A_w_inter_dy_stats_smooth[2], 'b-',linewidth=3)      
#    ax33.plot(days, A_d_inter_wmclimate_stats_smooth[0], 'r--',linewidth=3,label="Typical Growing Season") 
##    ax3.plot(days, A_d_inter_stats_smooth[2], 'r--',linewidth=3)    
#    ax33.plot(days, A_w_inter_wmclimate_stats_smooth[0], 'b--',linewidth=3,label="Typical Growing Season") 
##    ax3.plot(days, A_w_inter_stats_smooth[2], 'b--',linewidth=3)    
##    ax3.fill_between(days, A_w_inter_stats_smooth[1], A_w_inter_stats_smooth[2],alpha=0.0,color='blue')
##    ax3.fill_between(days, A_w_inter_dy_stats_smooth[1], A_w_inter_dy_stats_smooth[2],alpha=0.3,color='blue')    
#    ax33.legend(loc="upper left",fontsize=15)
#    ax33.annotate('A', xy=get_axis_limits(ax33),fontsize=40)
#    
#    
#    #longer growing season
#    #dry meadow climate
#    ax23.plot(days, A_d_dmclimate_gsl_stats_smooth[0], 'r-',linewidth=3,label="Conservative Species Scenario\n$\Delta$ Average Cumulative Assimilation=%d%%"%((sum((np.array(A_d_dmclimate_gsl_stats_smooth[0])*3600*6)/1000000.*44.)-sum((np.array(A_d_inter_dmclimate_stats_smooth[0])*3600*6)/1000000.*44.))/sum((np.array(A_d_inter_dmclimate_stats_smooth[0])*3600*6)/1000000.*44.)*100.0)) 
##    ax24.plot(days, A_d_inter_dy_gsl_stats_smooth[2], 'r-',linewidth=3)  
##    ax24.fill_between(days, A_d_inter_stats_smooth[1], A_d_inter_stats_smooth[2],alpha=0.0,color='red')
##    ax24.fill_between(days, A_d_inter_dy_gsl_stats_smooth[1], A_d_inter_dy_gsl_stats_smooth[2],alpha=0.3,color='red') 
#    ax23.plot(days, A_w_dmclimate_gsl_stats_smooth[0], 'b-',linewidth=3,label="Acquisitive Species Scenario\n$\Delta$ Average Cumulative Assimilation=%d%%"%((sum((np.array(A_w_dmclimate_gsl_stats_smooth[0])*3600*6)/1000000.*44.)-sum((np.array(A_w_inter_dmclimate_stats_smooth[0])*3600*6)/1000000.*44.))/sum((np.array(A_w_inter_dmclimate_stats_smooth[0])*3600*6)/1000000.*44.)*100.0)) 
##    ax24.plot(days, A_w_inter_dy_gsl_stats_smooth[2], 'b-',linewidth=3)      
#    ax23.plot(days,A_d_inter_dmclimate_stats_smooth[0], 'r--',linewidth=3,label="Typical Growing Season") 
##    ax24.plot(days, A_d_inter_stats_smooth[2], 'r--',linewidth=3)    
#    ax23.plot(days, A_w_inter_dmclimate_stats_smooth[0], 'b--',linewidth=3,label="Typical Growing Season") 
##    ax24.plot(days, A_w_inter_stats_smooth[2], 'b--',linewidth=3)    
##    ax24.fill_between(days, A_w_inter_stats_smooth[1], A_w_inter_stats_smooth[2],alpha=0.0,color='blue')
##    ax24.fill_between(days, A_w_inter_dy_gsl_stats_smooth[1], A_w_inter_dy_gsl_stats_smooth[2],alpha=0.3,color='blue')    
#    ax23.legend(loc="upper left",fontsize=15)
#    ax23.annotate('A', xy=get_axis_limits(ax23),fontsize=40)    
#
##wet meadow climate
#    ax34.plot(days, A_d_wmclimate_gsl_stats_smooth[0], 'r-',linewidth=3,label="Conservative Species Scenario\n$\Delta$ Average Cumulative Assimilation=%d%%"%((sum((np.array(A_d_wmclimate_gsl_stats_smooth[0])*3600*6)/1000000.*44.)-sum((np.array(A_d_inter_wmclimate_stats_smooth[0])*3600*6)/1000000.*44.))/sum((np.array(A_d_inter_wmclimate_stats_smooth[0])*3600*6)/1000000.*44.)*100.0)) 
##    ax24.plot(days, A_d_inter_dy_gsl_stats_smooth[2], 'r-',linewidth=3)  
##    ax24.fill_between(days, A_d_inter_stats_smooth[1], A_d_inter_stats_smooth[2],alpha=0.0,color='red')
##    ax24.fill_between(days, A_d_inter_dy_gsl_stats_smooth[1], A_d_inter_dy_gsl_stats_smooth[2],alpha=0.3,color='red') 
#    ax34.plot(days, A_w_wmclimate_gsl_stats_smooth[0], 'b-',linewidth=3,label="Acquisitive Species Scenario\n$\Delta$ Average Cumulative Assimilation=%d%%"%((sum((np.array(A_w_wmclimate_gsl_stats_smooth[0])*3600*6)/1000000.*44.)-sum((np.array(A_w_inter_wmclimate_stats_smooth[0])*3600*6)/1000000.*44.))/sum((np.array(A_w_inter_wmclimate_stats_smooth[0])*3600*6)/1000000.*44.)*100.0)) 
##    ax24.plot(days, A_w_inter_dy_gsl_stats_smooth[2], 'b-',linewidth=3)      
#    ax34.plot(days, A_d_inter_wmclimate_stats_smooth[0], 'r--',linewidth=3,label="Typical Growing Season") 
##    ax24.plot(days, A_d_inter_stats_smooth[2], 'r--',linewidth=3)    
#    ax34.plot(days, A_w_inter_wmclimate_stats_smooth[0], 'b--',linewidth=3,label="Typical Growing Season") 
##    ax24.plot(days, A_w_inter_stats_smooth[2], 'b--',linewidth=3)    
##    ax24.fill_between(days, A_w_inter_stats_smooth[1], A_w_inter_stats_smooth[2],alpha=0.0,color='blue')
##    ax24.fill_between(days, A_w_inter_dy_gsl_stats_smooth[1], A_w_inter_dy_gsl_stats_smooth[2],alpha=0.3,color='blue')    
#    ax34.legend(loc="upper left",fontsize=15)
#    ax34.annotate('A', xy=get_axis_limits(ax34),fontsize=40) 
#    
#    #hotter temperatures
#    ax24.plot(days, A_d_dmclimate_temps_stats_smooth[0], 'r-',linewidth=3,label="Conservative Species Scenario\n$\Delta$ Average Cumulative Assimilation=%d%%"%((sum((np.array(A_d_dmclimate_temps_stats_smooth[0])*3600*6)/1000000.*44.)-sum((np.array(A_d_inter_dmclimate_stats_smooth[0])*3600*6)/1000000.*44.))/sum((np.array(A_d_inter_dmclimate_stats_smooth[0])*3600*6)/1000000.*44.)*100.0)) 
##    ax24.plot(days, A_d_inter_dy_gsl_stats_smooth[2], 'r-',linewidth=3)  
##    ax24.fill_between(days, A_d_inter_stats_smooth[1], A_d_inter_stats_smooth[2],alpha=0.0,color='red')
##    ax24.fill_between(days, A_d_inter_dy_gsl_stats_smooth[1], A_d_inter_dy_gsl_stats_smooth[2],alpha=0.3,color='red') 
#    ax24.plot(days, A_w_dmclimate_temps_stats_smooth[0], 'b-',linewidth=3,label="Acquisitive Species Scnario\n$\Delta$ Average Cumulative Assimilation=%d%%"%((sum((np.array(A_w_dmclimate_temps_stats_smooth[0])*3600*6)/1000000.*44.)-sum((np.array(A_w_inter_dmclimate_stats_smooth[0])*3600*6)/1000000.*44.))/sum((np.array(A_w_inter_dmclimate_stats_smooth[0])*3600*6)/1000000.*44.)*100.0)) 
##    ax24.plot(days, A_w_inter_dy_gsl_stats_smooth[2], 'b-',linewidth=3)      
#    ax24.plot(days, A_d_inter_dmclimate_stats_smooth[0], 'r--',linewidth=3,label="Typical Growing Season") 
##    ax24.plot(days, A_d_inter_stats_smooth[2], 'r--',linewidth=3)    
#    ax24.plot(days, A_w_inter_dmclimate_stats_smooth[0], 'b--',linewidth=3,label="Typical Growing Season") 
##    ax24.plot(days, A_w_inter_stats_smooth[2], 'b--',linewidth=3)    
##    ax24.fill_between(days, A_w_inter_stats_smooth[1], A_w_inter_stats_smooth[2],alpha=0.0,color='blue')
##    ax24.fill_between(days, A_w_inter_dy_gsl_stats_smooth[1], A_w_inter_dy_gsl_stats_smooth[2],alpha=0.3,color='blue')    
#    ax24.legend(loc="upper left",fontsize=15)
#    ax24.annotate('A', xy=get_axis_limits(ax24),fontsize=40)
#
#    ax35.plot(days, A_d_wmclimate_temps_stats_smooth[0], 'r-',linewidth=3,label="Conservative Species Scenario\n$\Delta$ Average Cumulative Assimilation=%d%%"%((sum((np.array(A_d_wmclimate_temps_stats_smooth[0])*3600*6)/1000000.*44.)-sum((np.array(A_d_inter_wmclimate_stats_smooth[0])*3600*6)/1000000.*44.))/sum((np.array(A_d_inter_wmclimate_stats_smooth[0])*3600*6)/1000000.*44.)*100.0)) 
##    ax24.plot(days, A_d_inter_dy_gsl_stats_smooth[2], 'r-',linewidth=3)  
##    ax24.fill_between(days, A_d_inter_stats_smooth[1], A_d_inter_stats_smooth[2],alpha=0.0,color='red')
##    ax24.fill_between(days, A_d_inter_dy_gsl_stats_smooth[1], A_d_inter_dy_gsl_stats_smooth[2],alpha=0.3,color='red') 
#    ax35.plot(days, A_w_wmclimate_temps_stats_smooth[0], 'b-',linewidth=3,label="Acquisitive Species Scenario\n$\Delta$ Average Cumulative Assimilation=%d%%"%((sum((np.array(A_w_wmclimate_temps_stats_smooth[0])*3600*6)/1000000.*44.)-sum((np.array(A_w_inter_wmclimate_stats_smooth[0])*3600*6)/1000000.*44.))/sum((np.array(A_w_inter_wmclimate_stats_smooth[0])*3600*6)/1000000.*44.)*100.0)) 
##    ax24.plot(days, A_w_inter_dy_gsl_stats_smooth[2], 'b-',linewidth=3)      
#    ax35.plot(days, A_d_inter_wmclimate_stats_smooth[0], 'r--',linewidth=3,label="Typical Growing Season") 
##    ax24.plot(days, A_d_inter_stats_smooth[2], 'r--',linewidth=3)    
#    ax35.plot(days, A_w_inter_wmclimate_stats_smooth[0], 'b--',linewidth=3,label="Typical Growing Season") 
##    ax24.plot(days, A_w_inter_stats_smooth[2], 'b--',linewidth=3)    
##    ax24.fill_between(days, A_w_inter_stats_smooth[1], A_w_inter_stats_smooth[2],alpha=0.0,color='blue')
##    ax24.fill_between(days, A_w_inter_dy_gsl_stats_smooth[1], A_w_inter_dy_gsl_stats_smooth[2],alpha=0.3,color='blue')    
#    ax35.legend(loc="upper left",fontsize=15)
#    ax35.annotate('A', xy=get_axis_limits(ax35),fontsize=40)
   
###------------------------------------------------------######   


#---Figure 4---#
    
#temperatures in wm and dm
    ax4A.plot(days,surtemp_dm,'k-',linewidth=5,alpha=0.5,label="Dry Meadow")
    ax4A.plot(days,surtemp_wm,'k-',linewidth=5,alpha=1.0,label="Wet Meadow")
    ax4A.legend(fontsize=15,loc='upper left')
    ax4A.annotate('A', xy=get_axis_limits(ax4A,scale2=0.80),fontsize=20,fontname="Times New Roman")
#vwc in wm and dm    
    ax4B.plot(days,vwc_dm,'k-',linewidth=5,alpha=0.5,label="Dry Meadow Environment")
    ax4B.plot(days,vwc_wm,'k-',linewidth=5,alpha=1.0,label="Wet Meadow Environment")
    ax4B.annotate('B', xy=get_axis_limits(ax4B, scale2=0.85),fontsize=20,fontname="Times New Roman")
#nitrogen in wm and dm
    ax4C.plot(days,na_dm_min_inter,'k-',alpha=0.5,linewidth=5,label="Dry Meadow Environment")
    ax4C.plot(days,na_dm_max_inter,'k-',alpha=0.5,linewidth=5)
    ax4C.plot(days,na_wm_min_inter,'k-',alpha=1.0,linewidth=5,label="Wet Meadow Environment")
    ax4C.plot(days,na_wm_max_inter,'k-',alpha=1.0,linewidth=5)
    ax4C.fill_between(days, na_dm_min_inter, na_dm_max_inter,alpha=0.2,color='black')
    ax4C.fill_between(days, na_wm_min_inter, na_wm_max_inter,alpha=0.5,color='black')
    ax4C.annotate('C', xy=get_axis_limits(ax4C,scale2=0.85),fontsize=20,fontname="Times New Roman")





##temperatures in normal and extended summer
#    ax19.plot(days,gs_smooth(np.average(np.array([surtemp_dm,surtemp_wm]),axis=0).tolist(),0,364),'k-',linewidth=5,label="Typical Growing Season")
#    ax19.plot(days,gs_smooth(np.average(np.array([surtemp_dm_dy,surtemp_wm_dy]),axis=0).tolist(),0,364),'r-',linewidth=5,label="Extended Summer Growing Season")
#    ax19.legend(fontsize=20,loc='upper left')
#    ax19.annotate('A', xy=get_axis_limits(ax19),fontsize=40)
#
##vwc in normal and extended summer    
#    ax20.plot(days,gs_smooth(np.average(np.array([vwc_dm,vwc_wm]),axis=0).tolist(),0,364,poly=5),'k-',linewidth=5,label="Typical Growing Season")
#    ax20.plot(days,gs_smooth(np.average(np.array([vwc_dm_dy,vwc_wm_dy]),axis=0).tolist(),0,364,poly=5),'r-',linewidth=5,label="Extended Summer Growing Season")
##    axT.legend(fontsize=20,loc='upper left')
#    ax20.annotate('B', xy=get_axis_limits(ax20),fontsize=40)
#
##na in a normal and extended summer
#    ax21.plot(days,na_dm_min_inter,'k-',linewidth=5,label="Typical Growing Season Dry Meadow")
#    ax21.plot(days,na_wm_min_inter,'b-',linewidth=5,label="Typical Growing Season Wet Meadow")
#    ax21.plot(days,na_dm_max_inter,'k-',linewidth=5)
#    ax21.plot(days,na_wm_max_inter,'b-',linewidth=5)
#    ax21.plot(days,na_dm_min_inter_dy,'k--',linewidth=5,label="Extended Summer Growing Season Dry Meadow")
#    ax21.plot(days,na_wm_min_inter_dy,'b--',linewidth=5,label="Extended Summer Growing Season Wet Meadow")
#    ax21.plot(days,na_dm_max_inter_dy,'k--',linewidth=5)
#    ax21.plot(days,na_wm_max_inter_dy,'b--',linewidth=5)
##    ax21.fill_between(days, na_dm_min_inter, na_dm_max_inter,alpha=0.3,color='black')
##    ax21.fill_between(days, na_wm_min_inter,na_wm_max_inter,alpha=0.3,color='blue')
#    ax21.legend(fontsize=20,loc='upper left')
#    ax21.annotate('C', xy=get_axis_limits(ax21),fontsize=40)


###------------------------------------------------------######  

#####------different climates A time series------######
    ax7A.plot(days, A_d_inter_stats_smooth[1], 'k-',alpha=0.5,linewidth=5,label="Conservative Species") 
    ax7A.plot(days, A_d_inter_stats_smooth[2], 'k-',alpha=0.5,linewidth=5)    
    ax7A.fill_between(days, A_d_inter_stats_smooth[1], A_d_inter_stats_smooth[2],alpha=0.2,color='black')
    ax7A.plot(days, A_w_inter_stats_smooth[1], 'k-',alpha=1.0,linewidth=5,label="Acquisitive Species") 
    ax7A.plot(days, A_w_inter_stats_smooth[2], 'k-',alpha=1.0,linewidth=5)    
    ax7A.fill_between(days, A_w_inter_stats_smooth[1], A_w_inter_stats_smooth[2],alpha=0.5,color='black')
    ax7A.legend(loc="upper left",fontsize=15)
    ax7A.annotate('A', xy=get_axis_limits(ax7A,scale2=0.88),fontsize=20,fontname="Times New Roman")


#####--------diff climates nue time series----###    

    ax7B.plot(days, nue_d_inter_stats_smooth[1], 'k-',alpha=0.5,linewidth=5,label="Conservative Species") 
    ax7B.plot(days, nue_d_inter_stats_smooth[2], 'k-',alpha=0.5,linewidth=5)    
    ax7B.fill_between(days, nue_d_inter_stats_smooth[1], nue_d_inter_stats_smooth[2],alpha=0.2,color='black')
    ax7B.plot(days, nue_w_inter_stats_smooth[1], 'k-',alpha=1.0,linewidth=5,label="Acquisitive Species") 
    ax7B.plot(days, nue_w_inter_stats_smooth[2], 'k-',alpha=1.0,linewidth=5)    
    ax7B.fill_between(days, nue_w_inter_stats_smooth[1], nue_w_inter_stats_smooth[2],alpha=0.5,color='black')    
    ax7B.annotate('B', xy=get_axis_limits(ax7B,scale2=0.88),fontsize=20,fontname="Times New Roman")


######-------dif climates wue time series-----#### 
    ax7C.plot(days, wue_d_inter_stats_smooth[1], 'k-',alpha=0.5,linewidth=5,label="Dry Meadow Conservative Species") 
    ax7C.plot(days, wue_d_inter_stats_smooth[2], 'k-',alpha=0.5,linewidth=5)    
    ax7C.fill_between(days, wue_d_inter_stats_smooth[1], wue_d_inter_stats_smooth[2],alpha=0.2,color='black')
    ax7C.plot(days, wue_w_inter_stats_smooth[1], 'k-',alpha=1.0,linewidth=5,label="Wet Meadow Acquisitive Species") 
    ax7C.plot(days, wue_w_inter_stats_smooth[2], 'k-',alpha=1.0,linewidth=5)    
    ax7C.fill_between(days, wue_w_inter_stats_smooth[1], wue_w_inter_stats_smooth[2],alpha=0.5,color='black')      
    ax7C.annotate('C', xy=get_axis_limits(ax7C,scale2=0.88),fontsize=20,fontname="Times New Roman")



#####-----dry meadow climate nue time series-----###
#    ax10.plot(days, nue_d_inter_dmclimate_stats_smooth[1], 'r-',linewidth=5,label="Dry Meadow Conservative Species") 
#    ax10.plot(days, nue_d_inter_dmclimate_stats_smooth[2], 'r-',linewidth=5)    
#    ax10.fill_between(days, nue_d_inter_dmclimate_stats_smooth[1], nue_d_inter_dmclimate_stats_smooth[2],alpha=0.3,color='red')
#    ax10.plot(days, nue_w_inter_dmclimate_stats_smooth[1], 'b-',linewidth=5,label="Wet Meadow Acquisitive Species") 
#    ax10.plot(days, nue_w_inter_dmclimate_stats_smooth[2], 'b-',linewidth=5)    
#    ax10.fill_between(days, nue_w_inter_dmclimate_stats_smooth[1], nue_w_inter_dmclimate_stats_smooth[2],alpha=0.3,color='blue')    
##    axN2.legend(loc="upper left",fontsize=20)
#    ax10.annotate('B', xy=get_axis_limits(ax10),fontsize=40)
#
#####-----wet meadow climate nue time series-----###
#    ax11.plot(days, nue_d_inter_wmclimate_stats_smooth[1], 'r-',linewidth=5,label="Dry Meadow Conservative Species") 
#    ax11.plot(days, nue_d_inter_wmclimate_stats_smooth[2], 'r-',linewidth=5)    
#    ax11.fill_between(days, nue_d_inter_wmclimate_stats_smooth[1], nue_d_inter_wmclimate_stats_smooth[2],alpha=0.3,color='red')
#    ax11.plot(days, nue_w_inter_wmclimate_stats_smooth[1], 'b-',linewidth=5,label="Wet Meadow Acquisitive Species") 
#    ax11.plot(days, nue_w_inter_wmclimate_stats_smooth[2], 'b-',linewidth=5)    
#    ax11.fill_between(days, nue_w_inter_wmclimate_stats_smooth[1], nue_w_inter_wmclimate_stats_smooth[2],alpha=0.3,color='blue')    
##    axN2.legend(loc="upper left",fontsize=20)
#    ax11.annotate('B', xy=get_axis_limits(ax11),fontsize=40)  


#
#
##------moist meadow climate wue time series-----###
#    ax12.plot(days, wue_d_inter_mmclimate_stats_smooth[1], 'r-',linewidth=5,label="Dry Meadow Conservative Species") 
#    ax12.plot(days, wue_d_inter_mmclimate_stats_smooth[2], 'r-',linewidth=5)    
#    ax12.fill_between(days, wue_d_inter_mmclimate_stats_smooth[1], wue_d_inter_mmclimate_stats_smooth[2],alpha=0.3,color='red')
#    ax12.plot(days, wue_w_inter_mmclimate_stats_smooth[1], 'b-',linewidth=5,label="Wet Meadow Acquisitive Species") 
#    ax12.plot(days, wue_w_inter_mmclimate_stats_smooth[2], 'b-',linewidth=5)    
#    ax12.fill_between(days, wue_w_inter_mmclimate_stats_smooth[1], wue_w_inter_mmclimate_stats_smooth[2],alpha=0.3,color='blue')      
##    axO2.legend(loc="upper left",fontsize=20)
#    ax12.annotate('C', xy=get_axis_limits(ax12),fontsize=40)



#####-----wet meadow climate A time series------###   
    ax8B.plot(days, A_d_inter_wmclimate_stats_smooth[1], 'k-',alpha=0.5,linewidth=5,label="Conservative Species") 
    ax8B.plot(days, A_d_inter_wmclimate_stats_smooth[2], 'k-',alpha=0.5,linewidth=5)    
    ax8B.fill_between(days, A_d_inter_wmclimate_stats_smooth[1], A_d_inter_wmclimate_stats_smooth[2],alpha=0.2,color='black')
    ax8B.plot(days, A_w_inter_wmclimate_stats_smooth[1], 'k-',alpha=1.0,linewidth=5,label="Acquisitive Species") 
    ax8B.plot(days, A_w_inter_wmclimate_stats_smooth[2], 'k-',alpha=1.0,linewidth=5)    
    ax8B.fill_between(days, A_w_inter_wmclimate_stats_smooth[1], A_w_inter_wmclimate_stats_smooth[2],alpha=0.5,color='black')
    ax8B.annotate('B', xy=get_axis_limits(ax8B),fontsize=20,fontname="Times New Roman")


#####-----dry meadow climate A time series------###   
    ax8A.plot(days, A_d_inter_dmclimate_stats_smooth[1], 'k-',alpha=0.5,linewidth=5,label="Conservative Species") 
    ax8A.plot(days, A_d_inter_dmclimate_stats_smooth[2], 'k-',alpha=0.5,linewidth=5)    
    ax8A.fill_between(days, A_d_inter_dmclimate_stats_smooth[1], A_d_inter_dmclimate_stats_smooth[2],alpha=0.2,color='black')
    ax8A.plot(days, A_w_inter_dmclimate_stats_smooth[1], 'k-',alpha=1.0,linewidth=5,label="Acquisitive Species") 
    ax8A.plot(days, A_w_inter_dmclimate_stats_smooth[2], 'k-',alpha=1.0,linewidth=5)    
    ax8A.fill_between(days, A_w_inter_dmclimate_stats_smooth[1], A_w_inter_dmclimate_stats_smooth[2],alpha=0.5,color='black')
    ax8A.legend(loc="upper left", fontsize=15)
    ax8A.annotate('A', xy=get_axis_limits(ax8A),fontsize=20,fontname="Times New Roman")




###-------variation in abiotic environment (vwc and temp)-----##### 
    list1_min_d, list2_min_d = (list(t) for t in zip(*sorted(zip(temp_dm[2], A_d_inter_vartemp_stats_smooth[1]))))
    list1_max_d, list2_max_d = (list(t) for t in zip(*sorted(zip(temp_dm[2], A_d_inter_vartemp_stats_smooth[2]))))
    list1_min_w, list2_min_w = (list(t) for t in zip(*sorted(zip(temp_wm[2], A_w_inter_vartemp_stats_smooth[1]))))
    list1_max_w, list2_max_w = (list(t) for t in zip(*sorted(zip(temp_wm[2], A_w_inter_vartemp_stats_smooth[2]))))
    ax13A.plot(list1_min_d,list2_min_d, 'k-',alpha=0.5,linewidth=5,label="Conservative Species") 
    ax13A.plot(list1_max_d,list2_max_d, 'k-',alpha=0.5,linewidth=5) 
    ax13A.fill_between(list1_min_d, list2_min_d, list2_max_d,alpha=0.2,color='black')    
    ax13A.plot(list1_min_w,list2_min_w, 'k-',alpha=1.0,linewidth=5,label="Acquisitive Species") 
    ax13A.plot(list1_max_w,list2_max_w, 'k-',alpha=1.0,linewidth=5) 
    ax13A.fill_between(list1_min_w, list2_min_w, list2_max_w,alpha=0.5,color='black')        
    #add lines to show where temperatures fall during growing season
    ax13A.text(0.8,ax13A.get_ylim()[1]*0.1,'Growing Season\n   Days 1-2',fontsize=15,fontname="Times New Roman") 
    ax13A.plot([5]*100,np.linspace(0,25,100),'k-',linewidth=4)
    ax13A.text(5.8,ax13A.get_ylim()[1]*0.1,'Growing Season\n   Days 3-34',fontsize=15,fontname="Times New Roman") 
    ax13A.plot([10]*100,np.linspace(0,25,100),'k-',linewidth=4)
    ax13A.text(10.8,ax13A.get_ylim()[1]*0.1,'Growing Season\n  Days 35-112',fontsize=15,fontname="Times New Roman")         
    ax13A.plot([15]*100,np.linspace(0,25,100),'k-',linewidth=4) 
    ax13A.legend(loc='upper left',ncol=2,fontsize=15)
    ax13A.annotate('A', xy=get_axis_limits(ax13A,scale2=0.88),fontsize=20,fontname="Times New Roman")

    
    list1_min_d, list2_min_d = (list(t) for t in zip(*sorted(zip(sm_dm[3], A_d_inter_varvwc_stats[1][0:50]+A_d_inter_varvwc_stats_smooth[1][50:364]))))
    list1_max_d, list2_max_d = (list(t) for t in zip(*sorted(zip(sm_dm[3], A_d_inter_varvwc_stats[2][0:50]+A_d_inter_varvwc_stats_smooth[2][50:364]))))
    list1_min_w, list2_min_w = (list(t) for t in zip(*sorted(zip(sm_wm[3], A_w_inter_varvwc_stats[1][0:30]+A_w_inter_varvwc_stats_smooth[1][30:364]))))
    list1_max_w, list2_max_w = (list(t) for t in zip(*sorted(zip(sm_wm[3], A_w_inter_varvwc_stats[2][0:30]+A_w_inter_varvwc_stats_smooth[2][30:364]))))
    ax13B.plot(list1_min_d,list2_min_d, 'k-',alpha=0.5,linewidth=5,label="Conservative Species") 
    ax13B.plot(list1_max_d,list2_max_d, 'k-',alpha=0.5,linewidth=5) 
    ax13B.fill_between(list1_min_d, list2_min_d, list2_max_d,alpha=0.2,color='black')    
    ax13B.plot(list1_min_w,list2_min_w, 'k-',alpha=1.0,linewidth=5,label="Acquisitive Species") 
    ax13B.plot(list1_max_w,list2_max_w, 'k-',alpha=1.0,linewidth=5) 
    ax13B.fill_between(list1_min_w, list2_min_w, list2_max_w,alpha=0.5,color='black')    
    ax13B.annotate('B', xy=get_axis_limits(ax13B,scale2=0.88),fontsize=20,fontname="Times New Roman")

    
    
    list1_min_d, list2_min_d = (list(t) for t in zip(*sorted(zip(na_lst[gs0_m:gsf+1], A_d_inter_varna_stats[1][gs0_m:gsf+1]))))
    list1_max_d, list2_max_d = (list(t) for t in zip(*sorted(zip(na_lst[gs0_m:gsf+1], A_d_inter_varna_stats[2][gs0_m:gsf+1]))))
    list1_min_w, list2_min_w = (list(t) for t in zip(*sorted(zip(na_lst[gs0_m:gsf+1], A_w_inter_varna_stats[1][gs0_m:gsf+1]))))
    list1_max_w, list2_max_w = (list(t) for t in zip(*sorted(zip(na_lst[gs0_m:gsf+1], A_w_inter_varna_stats[2][gs0_m:gsf+1]))))
    ax13C.plot(list1_min_d,gs_smooth(list2_min_d,0,len(list2_min_d)-1), 'k-',alpha=0.5,linewidth=5,label="Conservative Species") 
    ax13C.plot(list1_max_d,gs_smooth(list2_max_d,0,len(list2_max_d)-1), 'k-',alpha=0.5,linewidth=5) 
    ax13C.fill_between(list1_min_d, gs_smooth(list2_min_d,0,len(list2_min_d)-1), gs_smooth(list2_max_d,0,len(list2_max_d)-1),alpha=0.2,color='black')    
    ax13C.plot(list1_min_w,gs_smooth(list2_min_w,0,len(list2_min_w)-1), 'k-',alpha=1.0,linewidth=5,label="Acquisitive Species") 
    ax13C.plot(list1_max_w,gs_smooth(list2_max_w,0,len(list2_max_w)-1), 'k-',alpha=1.0,linewidth=5) 
    ax13C.fill_between(list1_min_w, gs_smooth(list2_min_w,0,len(list2_min_w)-1), gs_smooth(list2_max_w,0,len(list2_max_w)-1),alpha=0.5,color='black')    
    ax13C.annotate('C', xy=get_axis_limits(ax13C,scale1=0.97,scale2=0.88),fontsize=20,fontname="Times New Roman")
    
#    list1_min_d, list2_min_d = (list(t) for t in zip(*sorted(zip(temp_dm[2], nue_d_inter_vartemp_stats_smooth[1]))))
#    list1_max_d, list2_max_d = (list(t) for t in zip(*sorted(zip(temp_dm[2], nue_d_inter_vartemp_stats_smooth[2]))))
#    list1_min_w, list2_min_w = (list(t) for t in zip(*sorted(zip(temp_wm[2], nue_w_inter_vartemp_stats_smooth[1]))))
#    list1_max_w, list2_max_w = (list(t) for t in zip(*sorted(zip(temp_wm[2], nue_w_inter_vartemp_stats_smooth[2]))))
#    ax16.plot(list1_min_d,list2_min_d, 'r-',linewidth=5,label="Conservative Speciees") 
#    ax16.plot(list1_max_d,list2_max_d, 'r-',linewidth=5) 
#    ax16.fill_between(list1_min_d, list2_min_d, list2_max_d,alpha=0.3,color='red')    
#    ax16.plot(list1_min_w,list2_min_w, 'b-',linewidth=5,label="Acquisitive Species") 
#    ax16.plot(list1_max_w,list2_max_w, 'b-',linewidth=5) 
#    ax16.fill_between(list1_min_w, list2_min_w, list2_max_w,alpha=0.3,color='blue')  
#    #add lines to show where temperatures fall during growing season
#    ax16.text(0.2,ax16.get_ylim()[1]*0.2,'Growing Season Days\n             1-2',fontsize=25,fontname="Times New Roman") 
#    ax16.plot([5]*100,np.linspace(0,25,100),'k-',linewidth=8)
#    ax16.text(5.2,ax16.get_ylim()[1]*0.2,'Growing Season Days\n             3-34',fontsize=25,fontname="Times New Roman") 
#    ax16.plot([10]*100,np.linspace(0,25,100),'k-',linewidth=8)
#    ax16.text(10.2,ax16.get_ylim()[1]*0.2,'Growing Season Days\n            35-112',fontsize=25,fontname="Times New Roman")     
#    ax16.plot([15]*100,np.linspace(0,25,100),'k-',linewidth=8)     
#    ax16.annotate('B', xy=get_axis_limits(ax16, scale1=0.1),fontsize=40)
#
#
#    list1_min_d, list2_min_d = (list(t) for t in zip(*sorted(zip(sm_dm[3], nue_d_inter_varvwc_stats[1][0:35]+nue_d_inter_varvwc_stats_smooth[1][35:364]))))
#    list1_max_d, list2_max_d = (list(t) for t in zip(*sorted(zip(sm_dm[3], nue_d_inter_varvwc_stats[2][0:35]+nue_d_inter_varvwc_stats_smooth[2][35:364]))))
#    list1_min_w, list2_min_w = (list(t) for t in zip(*sorted(zip(sm_wm[3], nue_w_inter_varvwc_stats[1][0:35]+nue_w_inter_varvwc_stats_smooth[1][35:364]))))
#    list1_max_w, list2_max_w = (list(t) for t in zip(*sorted(zip(sm_wm[3], nue_w_inter_varvwc_stats[2][0:35]+nue_w_inter_varvwc_stats_smooth[2][35:364]))))
#    ax17.plot(list1_min_d,list2_min_d, 'r-',linewidth=5,label="Conservative Species") 
#    ax17.plot(list1_max_d,list2_max_d, 'r-',linewidth=5,label="Conservative") 
#    ax17.fill_between(list1_min_d, list2_min_d, list2_max_d,alpha=0.3,color='red')    
#    ax17.plot(list1_min_w,list2_min_w, 'b-',linewidth=5,label="Acquisitive") 
#    ax17.plot(list1_max_w,list2_max_w, 'b-',linewidth=5,label="Acquisitive") 
#    ax17.fill_between(list1_min_w, list2_min_w, list2_max_w,alpha=0.3,color='blue')    
#    ax17.annotate('B', xy=get_axis_limits(ax17, scale1=0.2,scale2=0.9),fontsize=40)
#
#    list1_min_d, list2_min_d = (list(t) for t in zip(*sorted(zip(na_lst[gs0_m:gsf+1], nue_d_inter_varna_stats[1][gs0_m:gsf+1]))))
#    list1_max_d, list2_max_d = (list(t) for t in zip(*sorted(zip(na_lst[gs0_m:gsf+1], nue_d_inter_varna_stats[2][gs0_m:gsf+1]))))
#    list1_min_w, list2_min_w = (list(t) for t in zip(*sorted(zip(na_lst[gs0_m:gsf+1], nue_w_inter_varna_stats[1][gs0_m:gsf+1]))))
#    list1_max_w, list2_max_w = (list(t) for t in zip(*sorted(zip(na_lst[gs0_m:gsf+1], nue_w_inter_varna_stats[2][gs0_m:gsf+1]))))
#    ax18.plot(list1_min_d,gs_smooth(list2_min_d,0,len(list2_min_d)-1), 'r-',linewidth=5,label="Conservative") 
#    ax18.plot(list1_max_d,gs_smooth(list2_max_d,0,len(list2_max_d)-1), 'r-',linewidth=5,label="Conservative") 
#    ax18.fill_between(list1_min_d, gs_smooth(list2_min_d,0,len(list2_min_d)-1), gs_smooth(list2_max_d,0,len(list2_max_d)-1),alpha=0.3,color='red')    
#    ax18.plot(list1_min_w,gs_smooth(list2_min_w,0,len(list2_min_w)-1), 'b-',linewidth=5,label="Acquisitive") 
#    ax18.plot(list1_max_w,gs_smooth(list2_max_w,0,len(list2_max_w)-1), 'b-',linewidth=5,label="Acquisitive") 
#    ax18.fill_between(list1_min_w, gs_smooth(list2_min_w,0,len(list2_min_w)-1), gs_smooth(list2_max_w,0,len(list2_max_w)-1),alpha=0.3,color='blue')    
#    ax18.annotate('B', xy=get_axis_limits(ax18, scale1=0.9,scale2=0.98),fontsize=40)
#
#    list1_min_d, list2_min_d = (list(t) for t in zip(*sorted(zip(sm_dm[3], wue_d_inter_varvwc_stats_smooth[1]))))
#    list1_max_d, list2_max_d = (list(t) for t in zip(*sorted(zip(sm_dm[3], wue_d_inter_varvwc_stats_smooth[2]))))
#    list1_min_w, list2_min_w = (list(t) for t in zip(*sorted(zip(sm_dm[3], wue_w_inter_varvwc_stats_smooth[1]))))
#    list1_max_w, list2_max_w = (list(t) for t in zip(*sorted(zip(sm_dm[3], wue_w_inter_varvwc_stats_smooth[2]))))
#    ax28.plot(list1_min_d,list2_min_d, 'r-',linewidth=5,label="Conservative") 
#    ax28.plot(list1_max_d,list2_max_d, 'r-',linewidth=5) 
#    ax28.fill_between(list1_min_d, list2_min_d, list2_max_d,alpha=0.3,color='red')    
#    ax28.plot(list1_min_w,list2_min_w, 'b-',linewidth=5,label="Acquisitive") 
#    ax28.plot(list1_max_w,list2_max_w, 'b-',linewidth=5) 
#    ax28.fill_between(list1_min_w, list2_min_w, list2_max_w,alpha=0.3,color='blue')    
#    ax28.annotate('C', xy=get_axis_limits(ax28, scale1=0.2,scale2=0.9),fontsize=40)
#
#
#    list1_min_d, list2_min_d = (list(t) for t in zip(*sorted(zip(temp_dm[2], wue_d_inter_vartemp_stats_smooth[1]))))
#    list1_max_d, list2_max_d = (list(t) for t in zip(*sorted(zip(temp_dm[2], wue_d_inter_vartemp_stats_smooth[2]))))
#    list1_min_w, list2_min_w = (list(t) for t in zip(*sorted(zip(temp_dm[2], wue_w_inter_vartemp_stats_smooth[1]))))
#    list1_max_w, list2_max_w = (list(t) for t in zip(*sorted(zip(temp_dm[2], wue_w_inter_vartemp_stats_smooth[2]))))
#    ax27.plot(list1_min_d,list2_min_d, 'r-',linewidth=5,label="Conservative") 
#    ax27.plot(list1_max_d,list2_max_d, 'r-',linewidth=5,label="Conservative") 
#    ax27.fill_between(list1_min_d, list2_min_d, list2_max_d,alpha=0.3,color='red')    
#    ax27.plot(list1_min_w,list2_min_w, 'b-',linewidth=5,label="Acquisitive") 
#    ax27.plot(list1_max_w,list2_max_w, 'b-',linewidth=5,label="Acquisitive") 
#    ax27.fill_between(list1_min_w, list2_min_w, list2_max_w,alpha=0.3,color='blue')
#
#    ax27.text(0.2,ax27.get_ylim()[1]*0.1,'Growing Season Days\n             1-2',fontsize=25,fontname="Times New Roman") 
#    ax27.plot([5]*100,np.linspace(0,25,100),'k-',linewidth=8)
#    ax27.text(5.2,ax27.get_ylim()[1]*0.1,'Growing Season Days\n             3-34',fontsize=25,fontname="Times New Roman") 
#    ax27.plot([10]*100,np.linspace(0,25,100),'k-',linewidth=8)
#    ax27.text(10.2,ax27.get_ylim()[1]*0.1,'Growing Season Days\n            35-112',fontsize=25,fontname="Times New Roman")         
#    ax27.plot([15]*100,np.linspace(0,25,100),'k-',linewidth=8)    
#    ax27.annotate('C', xy=get_axis_limits(ax27, scale1=0.1),fontsize=40)
#
#
#    list1_min_d, list2_min_d = (list(t) for t in zip(*sorted(zip(na_lst[gs0_m:gsf+1], wue_d_inter_varna_stats[1][gs0_m:gsf+1]))))
#    list1_max_d, list2_max_d = (list(t) for t in zip(*sorted(zip(na_lst[gs0_m:gsf+1], wue_d_inter_varna_stats[2][gs0_m:gsf+1]))))
#    list1_min_w, list2_min_w = (list(t) for t in zip(*sorted(zip(na_lst[gs0_m:gsf+1], wue_w_inter_varna_stats[1][gs0_m:gsf+1]))))
#    list1_max_w, list2_max_w = (list(t) for t in zip(*sorted(zip(na_lst[gs0_m:gsf+1], wue_w_inter_varna_stats[2][gs0_m:gsf+1]))))
#    ax29.plot(list1_min_d,gs_smooth(list2_min_d,0,len(list2_min_d)-1), 'r-',linewidth=5,label="Conservative") 
#    ax29.plot(list1_max_d,gs_smooth(list2_max_d,0,len(list2_max_d)-1), 'r-',linewidth=5,label="Conservative") 
#    ax29.fill_between(list1_min_d, gs_smooth(list2_min_d,0,len(list2_min_d)-1), gs_smooth(list2_max_d,0,len(list2_max_d)-1),alpha=0.3,color='red')    
#    ax29.plot(list1_min_w,gs_smooth(list2_min_w,0,len(list2_min_w)-1), 'b-',linewidth=5,label="Acquisitive") 
#    ax29.plot(list1_max_w,gs_smooth(list2_max_w,0,len(list2_max_w)-1), 'b-',linewidth=5,label="Acquisitive") 
#    ax29.fill_between(list1_min_w, gs_smooth(list2_min_w,0,len(list2_min_w)-1), gs_smooth(list2_max_w,0,len(list2_max_w)-1),alpha=0.3,color='blue')    
#    ax29.annotate('C', xy=get_axis_limits(ax29, scale1=0.95,scale2=0.9),fontsize=40)
#    
    
#    ll=[nue_w_inter_dmclimate_stats_smooth[0],nue_w_inter_wmclimate_stats_smooth[0],nue_w_dmclimate_dy_stats_smooth[0],nue_w_wmclimate_dy_stats_smooth[0],nue_w_dmclimate_gsl_stats_smooth[0],nue_w_wmclimate_gsl_stats_smooth[0],nue_w_dmclimate_temps_stats_smooth[0],nue_w_wmclimate_temps_stats_smooth[0]]
#    ll2=[nue_d_inter_dmclimate_stats_smooth[0],nue_d_inter_wmclimate_stats_smooth[0],nue_d_dmclimate_dy_stats_smooth[0],nue_d_wmclimate_dy_stats_smooth[0],nue_d_dmclimate_gsl_stats_smooth[0],nue_d_wmclimate_gsl_stats_smooth[0],nue_d_dmclimate_temps_stats_smooth[0],nue_d_wmclimate_temps_stats_smooth[0]]


#first 30 day changes


    ll=[A_w_inter_dmclimate_stats_smooth[0],A_w_inter_wmclimate_stats_smooth[0],A_w_dmclimate_dy_stats_smooth[0],A_w_wmclimate_dy_stats_smooth[0],A_w_dmclimate_gsl_stats_smooth[0],A_w_wmclimate_gsl_stats_smooth[0],A_w_dmclimate_temps_stats_smooth[0],A_w_wmclimate_temps_stats_smooth[0]]
    ll2=[A_d_inter_dmclimate_stats_smooth[0],A_d_inter_wmclimate_stats_smooth[0],A_d_dmclimate_dy_stats_smooth[0],A_d_wmclimate_dy_stats_smooth[0],A_d_dmclimate_gsl_stats_smooth[0],A_d_wmclimate_gsl_stats_smooth[0],A_d_dmclimate_temps_stats_smooth[0],A_d_wmclimate_temps_stats_smooth[0]]

    newlist=[]
    newlist2=[]
    for ii in range(len(ll)):
        i=1
        i2=1
        mylist=[]
        mylist2=[]
        for y in ll[ii]:
            if y>0.0:
                i+=1
                mylist+=[y]                
                if i>30:
                    break
        for y in ll2[ii]:
            if y>0.0:
                i2+=1
                mylist2+=[y]                
                if i2>30:
                    break

        newlist+=[np.mean(mylist)]
        newlist2+=[np.mean(mylist2)]
        
    print np.array(newlist2)/np.array(newlist)



#what is the change in assimilation for different scenarios?

    #dm es, wm es, dm lgs, wm lgs, dm ht, wm ht
    
    dm_es_con=int((sum((np.array(A_d_dmclimate_dy_stats_smooth[0])*3600*6)/1000000.*44.)-sum((np.array(A_d_inter_dmclimate_stats_smooth[0])*3600*6)/1000000.*44.))/sum((np.array(A_d_inter_dmclimate_stats_smooth[0])*3600*6)/1000000.*44.)*100.0)
    dm_es_acq=int((sum((np.array(A_w_dmclimate_dy_stats_smooth[0])*3600*6)/1000000.*44.)-sum((np.array(A_w_inter_dmclimate_stats_smooth[0])*3600*6)/1000000.*44.))/sum((np.array(A_w_inter_dmclimate_stats_smooth[0])*3600*6)/1000000.*44.)*100.0)
    wm_es_con=int((sum((np.array(A_d_wmclimate_dy_stats_smooth[0])*3600*6)/1000000.*44.)-sum((np.array(A_d_inter_wmclimate_stats_smooth[0])*3600*6)/1000000.*44.))/sum((np.array(A_d_inter_wmclimate_stats_smooth[0])*3600*6)/1000000.*44.)*100.0)
    wm_es_acq=int((sum((np.array(A_w_wmclimate_dy_stats_smooth[0])*3600*6)/1000000.*44.)-sum((np.array(A_w_inter_wmclimate_stats_smooth[0])*3600*6)/1000000.*44.))/sum((np.array(A_w_inter_wmclimate_stats_smooth[0])*3600*6)/1000000.*44.)*100.0)
    
    dm_lgs_con=int((sum((np.array(A_d_dmclimate_gsl_stats_smooth[0])*3600*6)/1000000.*44.)-sum((np.array(A_d_inter_dmclimate_stats_smooth[0])*3600*6)/1000000.*44.))/sum((np.array(A_d_inter_dmclimate_stats_smooth[0])*3600*6)/1000000.*44.)*100.0)
    dm_lgs_acq=int((sum((np.array(A_w_dmclimate_gsl_stats_smooth[0])*3600*6)/1000000.*44.)-sum((np.array(A_w_inter_dmclimate_stats_smooth[0])*3600*6)/1000000.*44.))/sum((np.array(A_w_inter_dmclimate_stats_smooth[0])*3600*6)/1000000.*44.)*100.0)
    wm_lgs_con=int((sum((np.array(A_d_wmclimate_gsl_stats_smooth[0])*3600*6)/1000000.*44.)-sum((np.array(A_d_inter_wmclimate_stats_smooth[0])*3600*6)/1000000.*44.))/sum((np.array(A_d_inter_wmclimate_stats_smooth[0])*3600*6)/1000000.*44.)*100.0)
    wm_lgs_acq=int((sum((np.array(A_w_wmclimate_gsl_stats_smooth[0])*3600*6)/1000000.*44.)-sum((np.array(A_w_inter_wmclimate_stats_smooth[0])*3600*6)/1000000.*44.))/sum((np.array(A_w_inter_wmclimate_stats_smooth[0])*3600*6)/1000000.*44.)*100.0)
    
    dm_ht_con=int((sum((np.array(A_d_dmclimate_temps_stats_smooth[0])*3600*6)/1000000.*44.)-sum((np.array(A_d_inter_dmclimate_stats_smooth[0])*3600*6)/1000000.*44.))/sum((np.array(A_d_inter_dmclimate_stats_smooth[0])*3600*6)/1000000.*44.)*100.0)
    dm_ht_acq=int((sum((np.array(A_w_dmclimate_temps_stats_smooth[0])*3600*6)/1000000.*44.)-sum((np.array(A_w_inter_dmclimate_stats_smooth[0])*3600*6)/1000000.*44.))/sum((np.array(A_w_inter_dmclimate_stats_smooth[0])*3600*6)/1000000.*44.)*100.0)
    wm_ht_con=int((sum((np.array(A_d_wmclimate_temps_stats_smooth[0])*3600*6)/1000000.*44.)-sum((np.array(A_d_inter_wmclimate_stats_smooth[0])*3600*6)/1000000.*44.))/sum((np.array(A_d_inter_wmclimate_stats_smooth[0])*3600*6)/1000000.*44.)*100.0)
    wm_ht_acq=int((sum((np.array(A_w_wmclimate_temps_stats_smooth[0])*3600*6)/1000000.*44.)-sum((np.array(A_w_inter_wmclimate_stats_smooth[0])*3600*6)/1000000.*44.))/sum((np.array(A_w_inter_wmclimate_stats_smooth[0])*3600*6)/1000000.*44.)*100.0)

    #CUSTOMIZE TICK PARAMETERS ON FIGURES

    axes=[ax36,ax4A,ax4B,ax4C,ax7A,ax7B,ax7C,ax8A,ax8B,ax13A,ax13B,ax13C]
    
    for i in range(len(axes)):
        axes[i].tick_params(axis='y', labelsize=15)
        axes[i].tick_params(axis='x', labelsize=15)
        
        for tick in axes[i].get_xticklabels():
            tick.set_fontname("Times New Roman")
        for tick in axes[i].get_yticklabels():
            tick.set_fontname("Times New Roman")
#        figs[i].tight_layout()


    
#---------------Finalize Figure---------------#    
    
#    fig1.savefig('NUE_DM_ES.png')
#    fig2.savefig('NUE_DM_TEMP.png')
#    fig3.savefig('A_DM_ES.png')
#    fig4.savefig('Temp_TimeSeries.png')    
#    fig5.savefig('vwc_TimeSeries.png')
#    fig6.savefig('na_TimeSeries.png')
#    fig7.savefig('A_TimeSeries.png')
#    fig8.savefig('A_DM.png')
#    fig9.savefig('NUE_TimeSeries.png')
#    fig10.savefig('NUE_DM.png') 
#    fig11.savefig('NUE_WM.png') 
#    fig12.savefig('A_WM.png')     
#    fig13.savefig('A_temp_var.png') 
#    fig14.savefig('A_vwc_var.png')        
#    fig15.savefig('A_na_var.png')   
#    fig16.savefig('NUE_temp_var.png')       
#    fig17.savefig('NUE_vwc_var.png')
#    fig18.savefig('NUE_na_var.png')    
#    fig19.savefig('Temp_TimeSeries_DiffClim.png')   
#    fig20.savefig('vwc_TimeSeries_DiffClim.png')   
#    fig21.savefig('na_TimeSeries_DiffClim.png') 
#    fig22.savefig('NUE_DM_GSL.png')
#    fig23.savefig('A_DM_GSL.png')
#    fig24.savefig('A_DM_TEMP.png')
#    fig25.savefig('WUE_TimeSeries.png')
#    fig27.savefig('WUE_temp_var.png')
#    fig28.savefig('WUE_vwc_var.png')
#    fig29.savefig('WUE_na_var.png')
#    fig30.savefig('NUE_WM_GSL.png')
#    fig31.savefig('NUE_WM_ES.png')
#    fig32.savefig('NUE_WM_TEMP.png')
#    fig33.savefig('A_WM_ES.png')
#    fig34.savefig('A_WM_GSL.png')
#    fig35.savefig('A_WM_TEMP.png')
    

fig36.savefig('A_Extended_Summer.png')
fig4.savefig('temp_vwc_na.png')
fig7.savefig('A_NUE_WUE.png')
fig8.savefig('A_DM_WM.png')
fig13.savefig('A_temp_vwc_na.png')