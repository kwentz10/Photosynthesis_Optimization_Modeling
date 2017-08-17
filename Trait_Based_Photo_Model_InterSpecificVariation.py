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
from uncertain_params import monte_carlo_all

#Import photosynthesis model
from Photosynthesis_Model import photo_bound_meso_eqstom as photo

#Import functions to switch between Pa and umol/mol at sea level
from photo_functions import pa_con_atmfrac


#for doing letters on graphs for multiple plots in same figure
def get_axis_limits(ax, scale1=.95,scale2=0.9):
    return ax.get_xlim()[1]*scale1, ax.get_ylim()[1]*scale2


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
    
    #simulated and empirical NUE in plant communities 
    fig1, (ax1A,ax1B) = plt.subplots(2,figsize=(9,20),sharex=True)   
    ax1A.set_ylabel('NUE ($\mu$mol CO$_2$/g N s)',fontsize=20, fontname='Times New Roman')
    ax1A.set_ylim([0,5])
#    ax1A.set_title('NUE in Plant Communities', fontname='Times New Roman',fontsize=30)
    
    ax1B.set_ylabel('NUE (g biomass/g N)',fontsize=20, fontname='Times New Roman')
    ax1B.set_ylim([60,100])


     #-----figure 2----#
    
    #simulated and empirical WUE in plant communities
    fig2,(ax2A,ax2B,ax2C) = plt.subplots(3,figsize=(9,20),sharex=True,sharey=True) 
    ax2B.set_ylabel('WUE ($\mu$mol CO$_2$/mmol H$_2$O)',fontsize=20, fontname='Times New Roman')
    ax2A.set_ylim([0,5])
#    ax2A.set_title('WUE in Plant Communities', fontname='Times New Roman',fontsize=30)
   

    #--figure 3--#
    
    #simulated and empirical Assimilation in plant communities 
    fig3, (ax3A,ax3B) = plt.subplots(2,figsize=(9,20),sharex=True)   
    ax3A.set_ylabel('Assimilation ($\mu$mol CO$_2$/m$^2$s)',fontsize=20, fontname='Times New Roman')
    ax3A.set_ylim([0,30])
#    ax3A.set_title('Assimilation in Plant Communities', fontname='Times New Roman',fontsize=30)
    
    ax3B.set_ylabel('Growth Rate (g C/m$^2$ day)',fontsize=20, fontname='Times New Roman')
    ax3B.set_ylim([0,5])


    #--figure 4--#
    
    #leaf height vs. temperature
    fig4,ax4 = plt.subplots(figsize=(11,11))
    ax4.set_xlabel('Leaf Height (cm)',fontsize=25, fontname='Times New Roman')
    ax4.set_ylabel('Difference Between Leaf and Air Temperature ($^\circ$C)',fontsize=25, fontname='Times New Roman')
#    ax4.set_title('Leaf Height vs. Leaf & Air Temperature Difference', fontname='Times New Roman',fontsize=30)


    #-----figure 5----#

    #assimilation vs. stomatal conductance
    fig5,ax5 = plt.subplots(figsize=(11,11))

    ax5.set_xlabel('Assimilation ($\mu$mol CO$_2$/m$^2$s)',fontsize=25, fontname='Times New Roman')
    ax5.set_ylabel('Stomatal Conductance (mol CO$_2$/m$^2$s)',fontsize=25, fontname='Times New Roman')
#    ax5.set_title('Simulated Assimilation vs. Stomatal Conductance', fontname='Times New Roman',fontsize=30)
    ax5.set_ylim([0,0.55])
    ax5.set_xlim([0,25])
  
    
    #--figure 6--#
    
    #vpd vs. wue
    fig6, (ax6A,ax6B) = plt.subplots(2,figsize=(9,20),sharex=False)   
 #   ax6A.set_xlabel('VPD (cmol H$_2$O/mol air)',fontsize=25, fontname='Times New Roman')
    ax6A.set_ylabel('WUE ($\mu$mol CO$_2$/mmol H$_2$O)',fontsize=20, fontname='Times New Roman')
#    ax6A.set_title('Vapor Pressure Deficit vs. WUE', fontname='Times New Roman',fontsize=30)
    ax6B.set_xlabel('VPD (cmol H$_2$O/mol air)',fontsize=20, fontname='Times New Roman')
    ax6B.set_ylabel('WUE ($\mu$mol CO$_2$/mmol H$_2$O)',fontsize=20, fontname='Times New Roman')
    
    




    #---------------Initialize Arrays for Each Meadow---------------#
        
    #total nue and wue
    nue_tot=[]
    wue_tot=[]
   
    #wue and nue arrays
    wue_d=[]
    nue_d=[]
    wue_d_const=[]
    
    wue_m_const=[]

    wue_w=[]
    wue_w_const=[]
    nue_w=[]    
     

    
    
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
    ht_c=np.zeros(shape=1)+10.0 #Temperature of the Leaf (K)
    dia_c=np.zeros(shape=1)+(np.mean([1.4,2.3,2.6])/100.) #Mean diameter or size of leaf (m)
    na_c=np.zeros(shape=1)+(np.mean([2.5,5.6,6.3])) #leaf nitrogen (g N/ m2)
    t_c=np.zeros(shape=1)+15.0 #temp (C)
    

#---------------Import Variable Parameter Arrays from Leaf Parameter File---------------#
    params=monte_carlo_all()        
    
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
    
                z=0.2
                #---------------Photosynthesis Function---------------#
            
                #alter this line of code for when implementing different photosynthesis functions
                wue, nue, A, E, cs, ci, gsw, gs, gbw, gb, gm, cc,dd =photo(tk_25,ekc,eko,etau,ev,ej,toptv,toptj,na, qeff, PAR,tl,ea,chl,ij,kc25,ko25,o,ca,rh,m,a,frnr,flnr,ra,jm,g0,b,dia,u,q,vwc_min,vwc_max,vwc,z)
     
           
                if isinstance(wue, np.ndarray):
                    wue=wue[0]
        
                if isinstance(nue, np.ndarray):
                    nue=nue[0]
            
                if isinstance(A, np.ndarray):
                    A=A[0]        

                if isinstance(gs, np.ndarray):
                    gs=gs[0]      

                if isinstance(gsw, np.ndarray):
                    gsw=gsw[0]   

                if isinstance(E, np.ndarray):
                    E=E[0]   

                if isinstance(dd, np.ndarray):
                    dd=dd[0] 
            
                if isinstance(wue, list):
                    wue=wue[0]
            
                if isinstance(nue, list):
                    nue=nue[0]
            
                if isinstance(A, list):
                    A=A[0]      

                if isinstance(gs, list):
                    gs=gs[0]                     

                if isinstance(gsw, list):
                    gsw=gsw[0]     
                    
                if isinstance(E, list):
                    E=E[0]                         
                    
                if isinstance(dd, list):
                    dd=dd[0]  
                    
                wue_d+=[wue]
                nue_d+=[nue]
                gsw_d+=[gsw]
                A_d+=[A]
                E_d+=[E]
                vpd_d+=[dd]
                tl_d+=[tl]
                gs_d+=[gs]

                    

            elif xx==1:
                
                z=0.4
                #---------------Photosynthesis Function---------------#
            
                #alter this line of code for when implementing different photosynthesis functions
                wue, nue, A, E, cs, ci, gsw, gs, gbw, gb, gm, cc,dd =photo(tk_25,ekc,eko,etau,ev,ej,toptv,toptj,na, qeff, PAR,tl,ea,chl,ij,kc25,ko25,o,ca,rh,m,a,frnr,flnr,ra,jm,g0,b,dia,u,q,vwc_min,vwc_max,vwc,z)
     
                if isinstance(wue, np.ndarray):
                    wue=wue[0]
        
                if isinstance(nue, np.ndarray):
                    nue=nue[0]
            
                if isinstance(A, np.ndarray):
                    A=A[0]        

                if isinstance(gs, np.ndarray):
                    gs=gs[0]      

                if isinstance(gsw, np.ndarray):
                    gsw=gsw[0]   

                if isinstance(E, np.ndarray):
                    E=E[0]   

                if isinstance(dd, np.ndarray):
                    dd=dd[0] 
            
                if isinstance(wue, list):
                    wue=wue[0]
            
                if isinstance(nue, list):
                    nue=nue[0]
            
                if isinstance(A, list):
                    A=A[0]      

                if isinstance(gs, list):
                    gs=gs[0]                     

                if isinstance(gsw, list):
                    gsw=gsw[0]     
                    
                if isinstance(E, list):
                    E=E[0]                         
                    
                if isinstance(dd, list):
                    dd=dd[0] 
                    
                wue_m+=[wue]

                
     
            
            elif xx==2:
                
                z=0.4
                #---------------Photosynthesis Function---------------#
            
                #alter this line of code for when implementing different photosynthesis functions
                wue, nue, A, E, cs, ci, gsw, gs, gbw, gb, gm, cc,dd =photo(tk_25,ekc,eko,etau,ev,ej,toptv,toptj,na, qeff, PAR,tl,ea,chl,ij,kc25,ko25,o,ca,rh,m,a,frnr,flnr,ra,jm,g0,b,dia,u,q,vwc_min,vwc_max,vwc,z)
     
                if isinstance(wue, np.ndarray):
                    wue=wue[0]
        
                if isinstance(nue, np.ndarray):
                    nue=nue[0]
            
                if isinstance(A, np.ndarray):
                    A=A[0]        

                if isinstance(gs, np.ndarray):
                    gs=gs[0]      

                if isinstance(gsw, np.ndarray):
                    gsw=gsw[0]   

                if isinstance(E, np.ndarray):
                    E=E[0]   

                if isinstance(dd, np.ndarray):
                    dd=dd[0] 
            
                if isinstance(wue, list):
                    wue=wue[0]
            
                if isinstance(nue, list):
                    nue=nue[0]
            
                if isinstance(A, list):
                    A=A[0]      

                if isinstance(gs, list):
                    gs=gs[0]                     

                if isinstance(gsw, list):
                    gsw=gsw[0]     
                    
                if isinstance(E, list):
                    E=E[0]                         
                    
                if isinstance(dd, list):
                    dd=dd[0] 
                    
                wue_w+=[wue]
                nue_w+=[nue]
                gsw_w+=[gsw]
                A_w+=[A]
                E_w+=[E]
                vpd_w+=[dd]
                tl_w+=[tl]
                gs_w+=[gs]
                




#---------------Constant WUE---------------#
    params=monte_carlo_all()        
    
    for xx in range(len(params)):
        for yy in range(len(params[xx])):
            for key,val in params[xx][yy].items():
                exec(key + '=val')
        
        
           
            ht=ht_c

            temp=t_c

                
            
            
            #------calculate vapor pressure-----#
            pa_v=611*np.exp((17.27*temp)/(temp+237.3)) #saturation vapor pressure of air (Pa)
            ea_str=pa_con_atmfrac(pa_v,3528) #saturation vapor pressure of air (Pa-->umol h20/mol air)
            ea=rh*ea_str #vapor pressure (umol h2O/mol air)                
    
    
            #correct for leaf temperatures using leaf height
       
            t_diff=18-0.4*ht
        
            tl=temp+t_diff
    
            

            
            if xx==0: 
    
                z=0.2
                #---------------Photosynthesis Function---------------#
            
                #alter this line of code for when implementing different photosynthesis functions
                wue, nue, A, E, cs, ci, gsw, gs, gbw, gb, gm, cc,dd =photo(tk_25,ekc,eko,etau,ev,ej,toptv,toptj,na, qeff, PAR,tl,ea,chl,ij,kc25,ko25,o,ca,rh,m,a,frnr,flnr,ra,jm,g0,b,dia,u,q,vwc_min,vwc_max,vwc,z)
     
           
                if isinstance(wue, np.ndarray):
                    wue=wue[0]
        
                if isinstance(nue, np.ndarray):
                    nue=nue[0]
            
                if isinstance(A, np.ndarray):
                    A=A[0]        

                if isinstance(gs, np.ndarray):
                    gs=gs[0]      

                if isinstance(gsw, np.ndarray):
                    gsw=gsw[0]   

                if isinstance(E, np.ndarray):
                    E=E[0]   

                if isinstance(dd, np.ndarray):
                    dd=dd[0] 
            
                if isinstance(wue, list):
                    wue=wue[0]
            
                if isinstance(nue, list):
                    nue=nue[0]
            
                if isinstance(A, list):
                    A=A[0]      

                if isinstance(gs, list):
                    gs=gs[0]                     

                if isinstance(gsw, list):
                    gsw=gsw[0]     
                    
                if isinstance(E, list):
                    E=E[0]                         
                    
                if isinstance(dd, list):
                    dd=dd[0]  
                    
                wue_d_const+=[wue]
                
            if xx==1: 
    
                z=0.2
                #---------------Photosynthesis Function---------------#
            
                #alter this line of code for when implementing different photosynthesis functions
                wue, nue, A, E, cs, ci, gsw, gs, gbw, gb, gm, cc,dd =photo(tk_25,ekc,eko,etau,ev,ej,toptv,toptj,na, qeff, PAR,tl,ea,chl,ij,kc25,ko25,o,ca,rh,m,a,frnr,flnr,ra,jm,g0,b,dia,u,q,vwc_min,vwc_max,vwc,z)
     
           
                if isinstance(wue, np.ndarray):
                    wue=wue[0]
        
                if isinstance(nue, np.ndarray):
                    nue=nue[0]
            
                if isinstance(A, np.ndarray):
                    A=A[0]        

                if isinstance(gs, np.ndarray):
                    gs=gs[0]      

                if isinstance(gsw, np.ndarray):
                    gsw=gsw[0]   

                if isinstance(E, np.ndarray):
                    E=E[0]   

                if isinstance(dd, np.ndarray):
                    dd=dd[0] 
            
                if isinstance(wue, list):
                    wue=wue[0]
            
                if isinstance(nue, list):
                    nue=nue[0]
            
                if isinstance(A, list):
                    A=A[0]      

                if isinstance(gs, list):
                    gs=gs[0]                     

                if isinstance(gsw, list):
                    gsw=gsw[0]     
                    
                if isinstance(E, list):
                    E=E[0]                         
                    
                if isinstance(dd, list):
                    dd=dd[0]  
                    
                wue_m_const+=[wue]
    
            
            elif xx==2:
                
                z=0.4
                #---------------Photosynthesis Function---------------#
            
                #alter this line of code for when implementing different photosynthesis functions
                wue, nue, A, E, cs, ci, gsw, gs, gbw, gb, gm, cc,dd =photo(tk_25,ekc,eko,etau,ev,ej,toptv,toptj,na, qeff, PAR,tl,ea,chl,ij,kc25,ko25,o,ca,rh,m,a,frnr,flnr,ra,jm,g0,b,dia,u,q,vwc_min,vwc_max,vwc,z)
     
                if isinstance(wue, np.ndarray):
                    wue=wue[0]
        
                if isinstance(nue, np.ndarray):
                    nue=nue[0]
            
                if isinstance(A, np.ndarray):
                    A=A[0]        

                if isinstance(gs, np.ndarray):
                    gs=gs[0]      

                if isinstance(gsw, np.ndarray):
                    gsw=gsw[0]   

                if isinstance(E, np.ndarray):
                    E=E[0]   

                if isinstance(dd, np.ndarray):
                    dd=dd[0] 
            
                if isinstance(wue, list):
                    wue=wue[0]
            
                if isinstance(nue, list):
                    nue=nue[0]
            
                if isinstance(A, list):
                    A=A[0]      

                if isinstance(gs, list):
                    gs=gs[0]                     

                if isinstance(gsw, list):
                    gsw=gsw[0]     
                    
                if isinstance(E, list):
                    E=E[0]                         
                    
                if isinstance(dd, list):
                    dd=dd[0] 
                    
                wue_w_const+=[wue]


    
#---------------Figure 1: Plant Communities vs. NUE ---------------#      


    #model simulations

    nue_bp=ax1A.boxplot([nue_d,nue_w], patch_artist=True, showmeans=True, showfliers=False)

    ax1A.set_xticks([1.0,2.0])
    ax1A.set_xticklabels(['Dry Meadow','Wet Meadow'],fontname='Times New Roman',fontsize=20)

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

    ax1A.annotate('A', xy=get_axis_limits(ax1A,scale1=0.95,scale2=0.88),fontsize=20,fontname="Times New Roman")

    #model validation
    
    ax1B.bar(ax1A.get_xticks(),[72,88], yerr=[2.08,3.81], edgecolor='black', align="center",width=0.2, color='black',alpha=0.5,error_kw={'ecolor':'black', 'lw':2, 'capsize':5, 'capthick':2})
    ax1B.text(0.97,76,"a",fontsize=20)
    ax1B.text(1.98,93.5,"b",fontsize=20)
    ax1B.annotate('B', xy=get_axis_limits(ax1B,scale2=0.95),fontsize=20,fontname="Times New Roman")  
    ax1B.set_xticklabels(['Dry Meadow','Wet Meadow'],fontname='Times New Roman',fontsize=25)



#---------------Figure 2: Plant Communities vs. WUE ---------------#      

    #model simulations

    wue_bp=ax2A.boxplot([wue_d,wue_w], patch_artist=True, showmeans=True, showfliers=False)

    ax2A.set_xticks([1.0,2.0])
    ax2A.set_xticklabels(['Dry Meadow','Moist Meadow'],fontname='Times New Roman',fontsize=25)

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

    ax2A.annotate('A', xy=get_axis_limits(ax2A,scale1=0.95,scale2=0.85),fontsize=20,fontname="Times New Roman")

    #model validation
        
    ax2B.bar([1.0,2.0],[1.49,1.57], yerr=[0.06,0.04], edgecolor='black', align="center",width=0.2, color='black',alpha=0.5,error_kw={'ecolor':'black', 'lw':2, 'capsize':5, 'capthick':2})    
    ax2B.text(0.97,2,"a",fontsize=20)
    ax2B.text(1.98,2,"a",fontsize=20)
    ax2B.annotate('B', xy=get_axis_limits(ax2B,scale1=0.95,scale2=0.85),fontsize=20,fontname="Times New Roman")  
    ax2B.set_xticklabels(['Dry Meadow','Moist Meadow'],fontname='Times New Roman',fontsize=20)


    #model simulations when temp and leaf height are constant
    wue_bp_2=ax2C.boxplot([wue_d_const,wue_m_const], patch_artist=True, showmeans=True, showfliers=False)

    ax2C.set_xticks([1.0,2.0])
    ax2C.set_xticklabels(['Dry Meadow','Moist Meadow'],fontname='Times New Roman',fontsize=30)

    for box in wue_bp_2['boxes']:
        #change outline color
        box.set(color='black',linewidth=2)
        #change fill color
        box.set(facecolor='black',alpha=0.2)

    for whisker in wue_bp_2['whiskers']:
        whisker.set(color='black',linewidth=2)
    
    for cap in wue_bp_2['caps']:
        cap.set(color='black',linewidth=2)

    for median in wue_bp_2['medians']:
        median.set(color='black', linewidth=2)

    for flier in wue_bp_2['fliers']:
        flier.set(marker='*',color='black',alpha=0.5)

    for means in wue_bp_2['means']:
        means.set(marker='o',markerfacecolor='black')    

    ax2C.annotate('C', xy=get_axis_limits(ax2C,scale1=0.95,scale2=0.85),fontsize=20,fontname="Times New Roman")
           

#---------------Figure 3: Plant Communities vs. Assimilation--------------- #     

    A_bp=ax3A.boxplot([A_d,A_w], patch_artist=True, showmeans=True, showfliers=False)

    ax3A.set_xticks([1.0,2.0])
    ax3A.set_xticklabels(['Dry Meadow','Wet Meadow'],fontname='Times New Roman')

    
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
        
    ax3A.annotate('A', xy=get_axis_limits(ax3A),fontsize=20,fontname="Times New Roman")
        
      
    ax3B.bar([1.0,2.0],[0.91,1.92], yerr=[0.44, 0.81], edgecolor='black', align="center",width=0.2, color='black',alpha=0.5,error_kw={'ecolor':'black', 'lw':2, 'capsize':5, 'capthick':2})
    ax3B.text(0.97,1.6,"a",fontsize=20)
    ax3B.text(1.98,3.0,"b",fontsize=20)
    ax3B.annotate('B', xy=get_axis_limits(ax3B),fontsize=20,fontname="Times New Roman")  
    ax3B.set_xticklabels(['Dry Meadow','Wet Meadow'],fontname='Times New Roman',fontsize=30)



#---------------Figure 4: leaf height vs. air temperature plot        
    ax4.scatter([25,25,25,9,9,5,5,5,5,3,3,3,3,10,10,10,2],[3.6,5.8,12.5,17.2,20.8,13.1,7.2,5.8,11.2,14.,17.7,15.5,24.,14.5,6.3,24.,30.],edgecolors="black",facecolors="black",marker='o',s=30)
    ax4.plot(np.unique([25,25,25,9,9,5,5,5,5,3,3,3,3,10,10,10,2]), np.poly1d(np.polyfit([25,25,25,9,9,5,5,5,5,3,3,3,3,10,10,10,2], [3.6,5.8,12.5,17.2,20.8,13.1,7.2,5.8,11.2,14.,17.7,15.5,24.,14.5,6.3,24.,30.], 1))(np.unique([25,25,25,9,9,5,5,5,5,3,3,3,3,10,10,10,2])),color="black",linewidth=3)


    

#---------------Figure 5: Regression Plot Assimilation vs. Stomatal Conductance---------------#   
    ax5.scatter(A_d+A_w,gs_d+gs_w,edgecolors='black',facecolors='black',marker='o',s=30)
    ax5.plot(np.unique(A_d+A_w), np.poly1d(np.polyfit(A_d+A_w, gs_d+gs_w, 1))(np.unique(A_d+A_w)),color='black',linewidth=3)


    
#---------------Figure 6: Plot VPD vs. WUE for validation (use all points rather than mean of points)-------------#
    def func(x, a, b, c):
        return a * np.exp(-b * x) + c
    xdata=np.array(vpd_d+vpd_w)/10000.
    ydata=np.array(wue_d+wue_w)
    L=sorted(zip(xdata,ydata),key=operator.itemgetter(0))
    new_x,new_y=zip(*L)
    popt, pcov = curve_fit(func, np.array(new_x), np.array(new_y))
    
    ax6B.plot(np.array(new_x), func(np.array(new_x), *popt), 'k-', linewidth=3,label='Model Simulation')
    ax6B.scatter(xdata,ydata,edgecolors='black',facecolors='black',marker='o',s=30)
    ax6B.annotate('B', xy=get_axis_limits(ax6B),fontsize=20,fontname="Times New Roman")



    def func(x, a, b, c):
        return a * np.exp(-b * x) + c
    xdata=np.array([13734.61196533664,13734.61196533664,13734.61196533664,13734.61196533664,21364.95194607922,21364.95194607922,21364.95194607922,21364.95194607922,28995.2919268218,28995.2919268218,28995.2919268218,28995.2919268218,13734.61196533664,13734.61196533664,13734.61196533664,13734.61196533664,21364.95194607922,21364.95194607922,21364.95194607922,21364.95194607922,28995.2919268218,28995.2919268218,28995.2919268218,28995.2919268218])/10000.
    ydata=np.array([11,12,11,3,5,4.5,1,5,4,2,1,.5,9,9,10,11,5,5,5,5,5,5,5,2])*.4091
    L=sorted(zip(xdata,ydata),key=operator.itemgetter(0))
    new_x,new_y=zip(*L)
    popt, pcov = curve_fit(func, np.array(new_x), np.array(new_y))
    
    ax6A.plot(np.array(new_x), func(np.array(new_x), *popt), 'k-', linewidth=3,label='fit')
    ax6A.scatter(xdata,ydata,edgecolors="black",facecolors="black",marker='o',s=30,label='Empirical Data')
    ax6A.annotate('A', xy=get_axis_limits(ax6A),fontsize=20,fontname="Times New Roman")
    

    
#-----------FINALIZE PLOT------#


    figs=[fig1,fig2,fig3,fig4,fig5,fig6]
    axes=[ax1A,ax2A,ax3A,ax1B,ax2B,ax3B,ax2C,ax4,ax5,ax6A,ax6B]
    
    for i in range(len(axes)):
        axes[i].tick_params(axis='y', labelsize=15)
        axes[i].tick_params(axis='x', labelsize=15)
        for tick in axes[i].get_xticklabels():
            tick.set_fontname("Times New Roman")
        for tick in axes[i].get_yticklabels():
            tick.set_fontname("Times New Roman")


    ax1B.set_xticklabels(['Dry Meadow','Wet Meadow'],fontname='Times New Roman',fontsize=20)
    ax2C.set_xticklabels(['Dry Meadow','Moist Meadow'],fontname='Times New Roman',fontsize=20)
    ax3B.set_xticklabels(['Dry Meadow','Wet Meadow'],fontname='Times New Roman',fontsize=20)

        
#    for i in range(len(figs)):        
#        figs[i].tight_layout()

#---------------Finalize Figure---------------#    

    ##---Save Figure--##
    fig1.savefig('NUE_val.png') 
    fig2.savefig('WUE_val.png')
    fig3.savefig('Assimilation_val.png')
    fig4.savefig('Leaf_Ht_Temp.png')
    fig5.savefig('Assimilation_vs_Conductance.png')
    fig6.savefig('VPD_vs_WUE.png')

