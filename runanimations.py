#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 25 10:14:51 2020

@author: ufbu
"""
import pandas as pd
import netCDF4 as nc
import matplotlib.pyplot as plt
import firn_animations as fa

#create_animations('_IterationSite_J_perturbations1', 'python_outputs/Site_J_perturbations_Iteration1/1989-2019/100layers/')
matlab_output_folder = '_IterationSite_J_perturbations1'
outputfolder_matlabsim = matlab_output_folder

path_to_outputs = 'python_outputs/Site_J_perturbations_Iteration1/1989-2019/100layers/'


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
    
    
ani_rho = fa.density_animation(rho_all, Depth_dataframe)
ani_rho.save(path_to_outputs+'.mp4')

#ani_rho.event_source.stop()
#del ani_rho
plt.clf()


ani_T = fa.temperature_animation(T_ice_dataframe, Depth_dataframe)
ani_T.save(path_to_outputs+'.mp4')

plt.clf() 