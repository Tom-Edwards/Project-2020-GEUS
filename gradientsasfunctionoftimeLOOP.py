#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  9 19:03:33 2020

@author: ufbu
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 25 13:56:17 2020

@author: ub
"""

import pandas as pd
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import pathlib
import datetime
from gfunct import even_array
from gfunct import NGN
from numba import jit

number_of_layers_int = 100
start_date = datetime.date(1989,1,1)

    
list_of_Delta_rhos=list(np.zeros(11))
list_of_Dates =list(np.zeros(11))
list_of_total_mean = list(np.zeros(11))
list_of_Delta_rhos_minus_total_mean = list(np.zeros(11))

count=0

#for thick perturbations
#for l in range(10,32,2): 

# for homogenous pertubationsx
for l in range(10,40,10):

#%%    
    #Define dictionary with all the information needed to specify which animation 
    #was analysed
    
    
    #for thick perturbations
    #simulation_inf = {'matlab_output_folder': 'PerturbationsWD2012_1989/PThick'+str(l), 'initial_density_profile': 'Thick'+str(l), 'number_of_layers': '100layers', 'time_period': '1989-2019'}
    
    #for thick perturbations    
    simulation_inf = {'matlab_output_folder': 'PerturbationsWD2012_1989/PHomo'+str(l), 'initial_density_profile': 'Homo'+str(l), 'number_of_layers': '100layers', 'time_period': '1989-2019'}
    
    
    outputfolder_matlabsim = simulation_inf['matlab_output_folder']
    initial_density_profile=simulation_inf['initial_density_profile']
    number_of_layers=simulation_inf['number_of_layers']
    time_period=simulation_inf['time_period']
    
    
    target_hr = 204525
    start_date = datetime.date(1989,1,1)
    target_hr_delta = datetime.timedelta(hours=target_hr)
    
    
    #Define path to where the output will be placed
    path_to_outputs='python_outputs/'+initial_density_profile+'/'+time_period+'/'+number_of_layers+'/'
    
    
    #Create folder if they dont exstis for the plots and the animations
    pathlib.Path(path_to_outputs+'plots/').mkdir(parents=True, exist_ok=True)
    pathlib.Path(path_to_outputs+'Animations/').mkdir(parents=True, exist_ok=True)
    
    
    
    
    
    
    #Adds this name to to output files (Animations and plots):
    name='Kan_U_'+initial_density_profile+'_'+time_period+'_'+number_of_layers
    
    
    #Load dataset rho_bin
    fn_rho = 'Output/'+outputfolder_matlabsim+'/rho_bin_1.nc'
    ds_rho = nc.Dataset(fn_rho)
    
    #Put rho into a dataframe
    rho = ds_rho['rho'][:]
    rho_dataframe = pd.DataFrame(rho)
    
    
    
    
    #Load dataset snowc
    fn_snowc = 'Output/'+outputfolder_matlabsim+'/snowc_bin_1.nc'
    ds_snowc = nc.Dataset(fn_snowc)
    
    
    #Put snowc into a dataframe
    snowc = ds_snowc['snowc'][:]
    snowc_dataframe = pd.DataFrame(snowc)
    
    
    
    #Load dataset snic
    fn_snic = 'Output/'+outputfolder_matlabsim+'/snic_bin_1.nc'
    ds_snic = nc.Dataset(fn_snic)
    
    
    #Put snowc into a dataframe
    snic = ds_snic['snic'][:]
    snic_dataframe = pd.DataFrame(snic)
    
    rho_all = (snowc_dataframe + snic_dataframe)/ (snowc_dataframe/rho_dataframe + snic_dataframe/917)
    rho_all_np = rho_all.to_numpy()
    
    #%%
    #Load dataset T_ice_bin_1
    fn_T_ice = 'Output/'+outputfolder_matlabsim+'/T_ice_bin_1.nc'
    ds_T_ice = nc.Dataset(fn_T_ice)
    
    #Put T_ice into a dataframe
    T_ice = ds_T_ice['T_ice'][:]
    T_ice_dataframe = pd.DataFrame(T_ice)
    
    
    #Put Depth into a dataframe
    Depth = ds_rho['Depth'][:]
    Depth_dataframe = pd.DataFrame(Depth)
    Depth_np= Depth_dataframe.to_numpy()
    
    #Load initial_density profile from data
    init_density = pd.read_csv('Input/Initial state/density/'+initial_density_profile+'.csv', sep=',', header='infer')
    
    
    #Load 2012_density profile from data
    density_2012 = pd.read_csv('Input/Initial state/density/RetMIP_density_KAN-U_2012(deleted_NaN).csv', sep=',', header=None)
    
    
    
    #Load 2013_density profile from data
    density_2013 = pd.read_csv('Input/Initial state/density/DensityProfile_KAN_U_2013_Machguth_et_al.csv', sep=';')
    
    
    #%%
    
    
    # =============================================================================
    # @jit
    # def even_array(rho_array, depth_array, output_depth=2000):
    #         #This block of code generates an array, with differences between the simulation
    #     #And the experimental data. The density profiles from simulation and data does not 
    #     #have the same dimensions, and the simkulation have density in terms of layers
    #     #instead of depth. This Block should take take of this by, creating an-
    #     #array, by adding data points with same density as the layers it is in.
    #     alt_rho_all = np.zeros(output_depth)
    #     j=0
    #     for i in range(output_depth):
    #         cond=0
    #         while cond==0:
    #             if j>number_of_layers_int:
    #                 cond=1
    #             elif i < depth_array[j]*100:
    #                 #diff_dens_2013[i] = density_2013.iloc[:,1][i] - rho_all.iloc[k][j]
    #                 alt_rho_all[i] = rho_array[j]
    #                 #sqr_diffs_2013= (density_2013.iloc[:,1][i] - rho_all.iloc[k][j])**2
    #                 cond = 1
    #             else:
    #                 j += 1
    #     return (alt_rho_all)
    # 
    # 
    # @jit
    # def NGN(alt_rho_all, spacing, output_depth_cm):
    #     rho_all_spacing=np.zeros(int(output_depth_cm/spacing))
    #     for i in range(0,output_depth_cm,spacing):
    #         rho_all_spacing[int(i/spacing)] = np.mean(alt_rho_all[i:i+spacing])
    #     return rho_all_spacing
    # =============================================================================
    
    
    # =============================================================================
    # 
    #     diff_dens_2013 = np.zeros(len(density_2013))
    #     alt_rho_all = np.zeros(len(density_2013))
    #     sqr_diffs_2013=0
    #     j=0
    #     for i in range(len(diff_dens_2013)):
    #         cond=0
    #         while cond==0:
    #             if j>number_of_layers_int:
    #                 cond=1
    #             elif i < Depth_dataframe.iloc[k][j]*100:
    #                 #diff_dens_2013[i] = density_2013.iloc[:,1][i] - rho_all.iloc[k][j]
    #                 alt_rho_all[i] = rho_all.iloc[k][j]
    #                 #sqr_diffs_2013= (density_2013.iloc[:,1][i] - rho_all.iloc[k][j])**2
    #                 cond = 1
    #             else:
    #                 j += 1
    # =============================================================================
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    time_scaling=100
    Number_of_time_points= int(np.floor(len(rho_all)/time_scaling))
    
    
    interval_0m_1m = np.zeros(Number_of_time_points)
    interval_1m_10m = np.zeros(Number_of_time_points)
    interval_10m_15m = np.zeros(Number_of_time_points)
    
    
    
    grad_10m_15m=np.zeros(Number_of_time_points)
    grad_5m_15m=np.zeros(Number_of_time_points)
    
    
    
    total_mean = np.zeros(Number_of_time_points)
    
    
    Delta_rhos_minus_total_mean = np.zeros(Number_of_time_points)
    
    
    
    Number_of_negative_gradients_sim_10cm=np.zeros(Number_of_time_points)
    Number_of_negative_gradients_sim_20cm=np.zeros(Number_of_time_points)
    Number_of_negative_gradients_sim_50cm=np.zeros(Number_of_time_points)
    Number_of_negative_gradients_sim_100cm=np.zeros(Number_of_time_points)
    Number_of_negative_gradients_sim_190cm=np.zeros(Number_of_time_points)
    
    
    
    grad_RMSD_10cm=np.zeros(Number_of_time_points)
    grad_RMSD_20cm=np.zeros(Number_of_time_points)
    grad_RMSD_50cm=np.zeros(Number_of_time_points)
    grad_RMSD_100cm=np.zeros(Number_of_time_points)
    grad_RMSD_190cm=np.zeros(Number_of_time_points)
    
    
    for K in range(Number_of_time_points):
        print(K/Number_of_time_points)
        k=K*time_scaling
        
    
    
                    
                    
                        
            
        
        #print(sqr_diffs_2012)
        
        
        
        
        
        
        #PLot to show, that the new density array actually fits the old one
        #plt.plot(alt_rho_all, density_2012[0], '*'), plt.plot(rho_all.iloc[target_hr], Depth_dataframe.iloc[target_hr], '.')
        
        
        
        
        
        #diff_dens_2012 =density_2012[1]- np.ones(len(density_2012))*300
        
    
        
    
    
    #plt.show()
    
    # =============================================================================
    # 
    # difference_array = np.zeros(len(density_2012))
    # for i in range(len(density_2012)):
    #     print(i, j)
    #     #Check if the given Depth from data mathces Depth from simulation
    #     #multiply by 100 to account for cm -> m conversion
    #     for j in range(number_of_layers_int):
    #         if i/100 < Depth_dataframe.iloc[target_hr][j] and i+1/100 > Depth_dataframe.iloc[target_hr][j]:
    #             difference_array[i] = density_2012.iloc[:,1][i] - rho_all.iloc[target_hr][j]
    # 
    #     
    # =============================================================================
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    # =============================================================================
    #     
    #         #This block of code generates an array, with differences between the simulation
    #     #And the experimental data. The density profiles from simulation and data does not 
    #     #have the same dimensions, and the simkulation have density in terms of layers
    #     #instead of depth. This Block should take take of this by, creating an-
    #     #array, by adding data points with same density as the layers it is in.
    #     diff_dens_2013 = np.zeros(len(density_2013))
    #     alt_rho_all = np.zeros(len(density_2013))
    #     sqr_diffs_2013=0
    #     j=0
    #     for i in range(len(diff_dens_2013)):
    #         cond=0
    #         while cond==0:
    #             if j>number_of_layers_int:
    #                 cond=1
    #             elif i < Depth_dataframe.iloc[k][j]*100:
    #                 #diff_dens_2013[i] = density_2013.iloc[:,1][i] - rho_all.iloc[k][j]
    #                 alt_rho_all[i] = rho_all.iloc[k][j]
    #                 #sqr_diffs_2013= (density_2013.iloc[:,1][i] - rho_all.iloc[k][j])**2
    #                 cond = 1
    #             else:
    #                 j += 1
    # =============================================================================
        
        alt_rho_all = even_array(rho_all_np[k], Depth_np[k], output_depth=len(density_2013))
        
        interval_0m_1m[K] = np.mean(alt_rho_all[:100])
        interval_1m_10m[K] = np.mean(alt_rho_all[100:1000])    
        interval_10m_15m[K] = np.mean(alt_rho_all[1000:1500])    
        
        grad_10m_15m[K] = np.mean(alt_rho_all[1000:1500])-np.mean(alt_rho_all[100:1000])
        grad_5m_15m[K] = np.mean(alt_rho_all[500:1500])-np.mean(alt_rho_all[100:500])
    
        #print(even_array(rho_all_np[k], Depth_np[k], output_depth=len(density_2013)) - alt_rho_all)
    
    
        total_mean[K] = np.mean(alt_rho_all[100:1500])  
        
        
        Delta_rhos_minus_total_mean[K] = grad_10m_15m[K] - total_mean[K]
        
        
        
        spacing=10
        Number_of_negative_gradients_sim_10cm[K]= sum(np.gradient(NGN(alt_rho_all, spacing, 1800))<0)
        
        spacing=20
        Number_of_negative_gradients_sim_20cm[K]= sum(np.gradient(NGN(alt_rho_all, spacing, 1800))<0)
        
        spacing=50
        Number_of_negative_gradients_sim_50cm[K]= sum(np.gradient(NGN(alt_rho_all, spacing, 1800))<0)
        
        spacing=100
        Number_of_negative_gradients_sim_100cm[K]= sum(np.gradient(NGN(alt_rho_all, spacing, 1800))<0)
        
        spacing=200
        Number_of_negative_gradients_sim_190cm[K]= sum(np.gradient(NGN(alt_rho_all, spacing, 1800))<0)
    # =============================================================================
    #     
    #     
    #     spacing=10
    #     rho_all_spacing_10cm=np.zeros(int(1900/spacing))
    #     meas_2013_10cm=np.zeros(int(1900/spacing))
    #     #Compute a meaned array for the density profiles:
    #     for i in range(0, 1900, spacing):
    #         rho_all_spacing_10cm[int(i/spacing)] = np.mean(alt_rho_all[i:i+spacing])
    #         meas_2013_10cm[int(i/spacing)] = np.mean(density_2013.iloc[:,1][i:i+spacing])
    #     #print(rho_all_spacing_10cm)
    # 
    # 
    # 
    # 
    # 
    # 
    # 
    #     grad_RMSD_10cm[K]= np.sum((np.gradient(rho_all_spacing_10cm)-np.gradient(meas_2013_10cm))**2)
    #     Number_of_negative_gradients_sim_10cm[K]= sum(np.gradient(rho_all_spacing_10cm)<0)
    # 
    # 
    #     
    #     spacing=20
    #     rho_all_spacing_20cm=np.zeros(int(1900/spacing))
    #     meas_2013_20cm=np.zeros(int(1900/spacing))
    #     #Compute a meaned array for the density profiles:
    #     for i in range(0, 1900, spacing):
    #         rho_all_spacing_20cm[int(i/spacing)] = np.mean(alt_rho_all[i:i+spacing])
    #         meas_2013_20cm[int(i/spacing)] = np.mean(density_2013.iloc[:,1][i:i+spacing])
    #     #print(rho_all_spacing_20cm)
    # 
    # 
    # 
    #     grad_RMSD_20cm[K]= np.sum((np.gradient(rho_all_spacing_20cm)-np.gradient(meas_2013_20cm))**2)
    #     Number_of_negative_gradients_sim_20cm[K]= sum(np.gradient(rho_all_spacing_20cm)<0) 
    #     
    #     
    #     
    #     
    #     
    #     
    #     spacing=50
    #     rho_all_spacing_50cm=np.zeros(int(1900/spacing))
    #     meas_2013_50cm=np.zeros(int(1900/spacing))
    #     #Compute a meaned array for the density profiles:
    #     for i in range(0, 1900, spacing):
    #         rho_all_spacing_50cm[int(i/spacing)] = np.mean(alt_rho_all[i:i+spacing])
    #         meas_2013_50cm[int(i/spacing)] = np.mean(density_2013.iloc[:,1][i:i+spacing])
    #     #print(rho_all_spacing_50cm)
    # 
    # 
    # 
    #     grad_RMSD_50cm[K]= np.sum((np.gradient(rho_all_spacing_50cm)-np.gradient(meas_2013_50cm))**2)
    #     Number_of_negative_gradients_sim_50cm[K]= sum(np.gradient(rho_all_spacing_50cm)<0)
    # 
    # 
    # 
    # 
    # 
    #     spacing=100
    #     rho_all_spacing_100cm=np.zeros(int(1900/spacing))
    #     meas_2013_100cm=np.zeros(int(1900/spacing))
    #     #Compute a meaned array for the density profiles:
    #     for i in range(0, 1900, spacing):
    #         rho_all_spacing_100cm[int(i/spacing)] = np.mean(alt_rho_all[i:i+spacing])
    #         meas_2013_100cm[int(i/spacing)] = np.mean(density_2013.iloc[:,1][i:i+spacing])
    # 
    # 
    # 
    #     grad_RMSD_100cm[K]= np.sum((np.gradient(rho_all_spacing_100cm)-np.gradient(meas_2013_100cm))**2)
    #     Number_of_negative_gradients_sim_100cm[K]= sum(np.gradient(rho_all_spacing_100cm)<0)
    #     
    #     
    #     
    #     
    #     spacing=190
    #     rho_all_spacing_190cm=np.zeros(int(1900/spacing))
    #     meas_2013_190cm=np.zeros(int(1900/spacing))
    #     #Compute a meaned array for the density profiles:
    #     for i in range(0, 1900, spacing):
    #         rho_all_spacing_190cm[int(i/spacing)] = np.mean(alt_rho_all[i:i+spacing])
    #         meas_2013_190cm[int(i/spacing)] = np.mean(density_2013.iloc[:,1][i:i+spacing])
    # 
    # 
    # 
    # 
    # 
    # 
    #     grad_RMSD_190cm[K]= np.sum((np.gradient(rho_all_spacing_190cm)-np.gradient(meas_2013_190cm))**2)
    #     Number_of_negative_gradients_sim_190cm[K]= sum(np.gradient(rho_all_spacing_190cm)<0)
    # 
    # =============================================================================
        #print(k, K/Number_of_time_points, Number_of_negative_gradients_sim_190cm[K])
    
    
    
    
    
    
    #The name of density profile to look nice on graphs
    initial_density_profile_for_plot_title=initial_density_profile
    #initial_density_profile_for_plot_title = '2012 site J composite'
    
   # initial_density_profile_for_plot_title = '2012 site J composite'
    
    
    
    #All the plots for the perturbations specific
# =============================================================================
#     fig = plt.figure(figsize=(8,6))
#     plt.plot([start_date+datetime.timedelta(hours=x*time_scaling) for x in range(Number_of_time_points)], Number_of_negative_gradients_sim_10cm, label='10cm')
#     plt.plot([start_date+datetime.timedelta(hours=x*time_scaling) for x in range(Number_of_time_points)], Number_of_negative_gradients_sim_20cm, label='20cm')
#     plt.plot([start_date+datetime.timedelta(hours=x*time_scaling) for x in range(Number_of_time_points)], Number_of_negative_gradients_sim_50cm, label='50cm')
#     plt.plot([start_date+datetime.timedelta(hours=x*time_scaling) for x in range(Number_of_time_points)], Number_of_negative_gradients_sim_100cm, label='100cm')
#     plt.plot([start_date+datetime.timedelta(hours=x*time_scaling) for x in range(Number_of_time_points)], Number_of_negative_gradients_sim_190cm, label='200cm')
#     #plt.legend(loc="lower center", bbox_to_anchor=(1.05, 0.3))
#     plt.axvline(x=start_date+datetime.timedelta(hours=(target_hr)), linestyle=':', color='grey', label='Date: '+str(datetime.date(1989,1,1)+datetime.timedelta(hours=target_hr)))
#     plt.axvline(x=start_date+datetime.timedelta(hours=(target_hr+8050)), linestyle='--', color='grey', label='Date: '+str(datetime.date(1989,1,1)+datetime.timedelta(hours=target_hr+8050)))
#     plt.legend()
#     plt.title(initial_density_profile_for_plot_title+ ' profile NGN (legend show sampling resolution)')
#     plt.ylabel(r'No. of negative gradients')
#     plt.xlabel(r'Years from '+str(start_date)+ ' to '+ str(datetime.timedelta(hours=Number_of_time_points*time_scaling)+start_date))
#     #plt.savefig('python_outputs/gradientsasfunctionoftime')
#     plt.savefig(path_to_outputs+'gradientsasfunctionoftime'+str(initial_density_profile)+'.png')
#     plt.close()
#     
#     
#     
#     
#     
#     
#     
#     plt.plot([start_date+datetime.timedelta(hours=x*time_scaling) for x in range(Number_of_time_points)], grad_10m_15m, label= r"$\Delta \rho = \rho_{(10m,15m]}-\rho_{[1m,10m]}$")
#     plt.axvline(x=start_date+datetime.timedelta(hours=(target_hr)), linestyle=':', color='grey', label='Date: '+str(datetime.date(1989,1,1)+datetime.timedelta(hours=target_hr)))
#     plt.axvline(x=start_date+datetime.timedelta(hours=(target_hr+8050)), linestyle='--', color='grey', label='Date: '+str(datetime.date(1989,1,1)+datetime.timedelta(hours=target_hr+8050)))
#     plt.legend(loc='upper left', prop={'size': 8})
#     plt.title(initial_density_profile_for_plot_title+' Difference in mean density'+ r' $\Delta \rho$')
#     plt.ylabel(r"Difference in mean density (${kg/m^3})$")
#     plt.xlabel(r'Years from '+str(start_date)+ ' to '+ str(datetime.timedelta(hours=Number_of_time_points*time_scaling)+start_date))
#     plt.savefig(path_to_outputs+'gradient_10m_15m'+str(initial_density_profile)+'.png')
#     plt.close()
#     
#     
#     plt.plot([start_date+datetime.timedelta(hours=x*time_scaling) for x in range(Number_of_time_points)], grad_5m_15m, label= r"$\Delta = \rho_{(5m,15m]}-\rho_{[1m,5m]}$")
#     plt.axvline(x=start_date+datetime.timedelta(hours=(target_hr)), linestyle=':', color='grey', label='Date: '+str(datetime.date(1989,1,1)+datetime.timedelta(hours=target_hr)))
#     plt.axvline(x=start_date+datetime.timedelta(hours=(target_hr+8050)), linestyle='--', color='grey', label='Date: '+str(datetime.date(1989,1,1)+datetime.timedelta(hours=target_hr+8050)))
#     plt.legend()
#     plt.title(initial_density_profile_for_plot_title+ ' profile of negative gradients')
#     plt.ylabel(r"Difference in mean density (${kg}/{m^3}$)")
#     plt.xlabel(r'Years from '+str(start_date)+ ' to '+ str(datetime.timedelta(hours=Number_of_time_points*time_scaling)+start_date))
#     plt.savefig(path_to_outputs+'grad_5m_15m'+str(initial_density_profile)+'.png')
#     plt.close()
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     plt.plot([start_date+datetime.timedelta(hours=x*time_scaling) for x in range(Number_of_time_points)], interval_0m_1m, label= r"$rho_{[0,1]}$")
#     plt.axvline(x=start_date+datetime.timedelta(hours=(target_hr)), linestyle=':', color='grey', label='Date: '+str(datetime.date(1989,1,1)+datetime.timedelta(hours=target_hr)))
#     plt.axvline(x=start_date+datetime.timedelta(hours=(target_hr+8050)), linestyle='--', color='grey', label='Date: '+str(datetime.date(1989,1,1)+datetime.timedelta(hours=target_hr+8050)))
#     plt.legend()
#     plt.title(initial_density_profile_for_plot_title+ ' mean for interval (0,1]')
#     plt.ylabel(r"mean density (${kg}/{m^3}$)")
#     plt.xlabel(r'Years from '+str(start_date)+ ' to '+ str(datetime.timedelta(hours=Number_of_time_points*time_scaling)+start_date))
#     plt.savefig(path_to_outputs+' mean for interval (0,1]'+str(initial_density_profile)+'.png')
#     plt.close()
#     
#     
#     plt.plot([start_date+datetime.timedelta(hours=x*time_scaling) for x in range(Number_of_time_points)], interval_1m_10m, label= r"$rho_{(1,10]}$")
#     plt.axvline(x=start_date+datetime.timedelta(hours=(target_hr)), linestyle=':', color='grey', label='Date: '+str(datetime.date(1989,1,1)+datetime.timedelta(hours=target_hr)))
#     plt.axvline(x=start_date+datetime.timedelta(hours=(target_hr+8050)), linestyle='--', color='grey', label='Date: '+str(datetime.date(1989,1,1)+datetime.timedelta(hours=target_hr+8050)))
#     plt.legend()
#     plt.title(initial_density_profile_for_plot_title+ ' mean for interval (1,10]')
#     plt.ylabel(r"mean density (${kg}/{m^3}$)")
#     plt.xlabel(r'Years from '+str(start_date)+ ' to '+ str(datetime.timedelta(hours=Number_of_time_points*time_scaling)+start_date))
#     plt.savefig(path_to_outputs+' mean for interval (1,10]'+str(initial_density_profile)+'.png')
#     plt.close()
#     
#     
#     plt.plot([start_date+datetime.timedelta(hours=x*time_scaling) for x in range(Number_of_time_points)], interval_10m_15m, label= r"$rho_{(10,15]}$")
#     plt.axvline(x=start_date+datetime.timedelta(hours=(target_hr)), linestyle=':', color='grey', label='Date: '+str(datetime.date(1989,1,1)+datetime.timedelta(hours=target_hr)))
#     plt.axvline(x=start_date+datetime.timedelta(hours=(target_hr+8050)), linestyle='--', color='grey', label='Date: '+str(datetime.date(1989,1,1)+datetime.timedelta(hours=target_hr+8050)))
#     plt.legend()
#     plt.title(initial_density_profile_for_plot_title+ ' mean for interval [10,15]')
#     plt.ylabel(r"mean density (${kg}/{m^3}$)")
#     plt.xlabel(r'Years from '+str(start_date)+ ' to '+ str(datetime.timedelta(hours=Number_of_time_points*time_scaling)+start_date))
#     plt.savefig(path_to_outputs+' mean for interval (10,15]'+str(initial_density_profile)+'.png')
#     plt.close()
#     
#     
#     
#     
#     
#     plt.plot([start_date+datetime.timedelta(hours=x*time_scaling) for x in range(Number_of_time_points)], interval_10m_15m, label= r"$rho_{[10,15]}$")
#     plt.plot([start_date+datetime.timedelta(hours=x*time_scaling) for x in range(Number_of_time_points)], interval_1m_10m, label= r"$rho_{[1,10]}$")
#     plt.plot([start_date+datetime.timedelta(hours=x*time_scaling) for x in range(Number_of_time_points)], interval_0m_1m, label= r"$rho_{[0,1]}$")
#     plt.axvline(x=start_date+datetime.timedelta(hours=(target_hr)), linestyle=':', color='grey', label='Date: '+str(datetime.date(1989,1,1)+datetime.timedelta(hours=target_hr)))
#     plt.axvline(x=start_date+datetime.timedelta(hours=(target_hr+8050)), linestyle='--', color='grey', label='Date: '+str(datetime.date(1989,1,1)+datetime.timedelta(hours=target_hr+8050)))
#     plt.ylim(250,900)
#     plt.legend(loc='lower left', prop={'size': 8})
#     plt.title(initial_density_profile_for_plot_title+ ' mean for different intervals')
#     plt.ylabel(r"mean density (${kg}/{m^3}$)")
#     plt.xlabel(r'Years from '+str(start_date)+ ' to '+ str(datetime.timedelta(hours=Number_of_time_points*time_scaling)+start_date))
#     plt.savefig(path_to_outputs+'intervals'+str(initial_density_profile)+'.png')
#     plt.close()
#     
# =============================================================================
    
    
    
    
    
    list_of_Delta_rhos[count]= np.copy(grad_10m_15m)
    list_of_Dates[count] = [start_date+datetime.timedelta(hours=x*time_scaling) for x in range(Number_of_time_points)]
    list_of_total_mean[count] = np.copy(total_mean)
    
    
    
    list_of_Delta_rhos_minus_total_mean = np.copy(Delta_rhos_minus_total_mean)
    
    
    
    count+=1
    




pertubation = 'homogenous'
 
    
    
    
#Plot Delta rho for more perturbs
plt.plot(list_of_Dates[0], list_of_Delta_rhos[0], label=pertubation+'10')    
plt.plot(list_of_Dates[1], list_of_Delta_rhos[1], label=pertubation+'20')
plt.plot(list_of_Dates[2], list_of_Delta_rhos[2], label=pertubation+'30')
plt.legend()
plt.title(' Difference in mean density'+ r' $\Delta \rho$')
plt.ylabel(r"Difference in mean density (${kg/m^3})$")
plt.xlabel(r'Years from '+str(start_date)+ ' to '+ str(datetime.timedelta(hours=Number_of_time_points*time_scaling)+start_date))
plt.axvline(x=start_date+datetime.timedelta(hours=(target_hr)), linestyle=':', color='grey', label='Date: '+str(datetime.date(1989,1,1)+datetime.timedelta(hours=target_hr)))
plt.axvline(x=start_date+datetime.timedelta(hours=(target_hr+8050)), linestyle='--', color='grey', label='Date: '+str(datetime.date(1989,1,1)+datetime.timedelta(hours=target_hr+8050)))
plt.savefig('python_outputs/'+pertubation+'DEltarho'+'.pdf')    
plt.close()
    
    
    
    
#Plot mean[:1500] rho for more perturbs
plt.plot(list_of_Dates[0], list_of_total_mean[0], label=pertubation+'10')    
plt.plot(list_of_Dates[1], list_of_total_mean[1], label=pertubation+'20')
plt.plot(list_of_Dates[2], list_of_total_mean[2], label=pertubation+'30')  
plt.ylabel(r"Mean density (${kg/m^3})$")
plt.xlabel(r'Years from '+str(start_date)+ ' to '+ str(datetime.timedelta(hours=Number_of_time_points*time_scaling)+start_date))
plt.axvline(x=start_date+datetime.timedelta(hours=(target_hr)), linestyle=':', color='grey', label='Date: '+str(datetime.date(1989,1,1)+datetime.timedelta(hours=target_hr)))
plt.axvline(x=start_date+datetime.timedelta(hours=(target_hr+8050)), linestyle='--', color='grey', label='Date: '+str(datetime.date(1989,1,1)+datetime.timedelta(hours=target_hr+8050)))
plt.legend()
plt.title('Mean density for the'+ r' (1:15] region')
plt.savefig('python_outputs/'+pertubation+'mean[1:15]'+'.pdf')    
plt.close()
    
    
    
    
    
# =============================================================================
# plt.plot(list_of_Dates[0], list_of_Delta_rhos_minus_total_mean[0], label=r'$Delta \rho'+': Thick10')    
# plt.plot(list_of_Dates[5], list_of_Delta_rhos_minus_total_mean[5], label=r'$Delta \rho'+': Thick120')
# plt.plot(list_of_Dates[10], list_of_Delta_rhos_minus_total_mean[10], label=r'$Delta \rho'+': Thick30')  
# plt.ylabel(r"Mean density (${kg/m^3})$")
# plt.xlabel(r'Years from '+str(start_date)+ ' to '+ str(datetime.timedelta(hours=Number_of_time_points*time_scaling)+start_date))
# plt.axvline(x=start_date+datetime.timedelta(hours=(target_hr)), linestyle=':', color='grey', label='Date: '+str(datetime.date(1989,1,1)+datetime.timedelta(hours=target_hr)))
# plt.axvline(x=start_date+datetime.timedelta(hours=(target_hr+8050)), linestyle='--', color='grey', label='Date: '+str(datetime.date(1989,1,1)+datetime.timedelta(hours=target_hr+8050)))
# plt.legend()
# plt.title(' Mean density '+ r'$ mean_{(1:15]} - \Delta \rho$')
# plt.savefig('python_outputs/'+'mean(:15] - minus Delta_rho'+'.pdf')    
# plt.close()   
# =============================================================================
    
    
    
    
    
    
    
  