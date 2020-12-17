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
import matplotlib.animation as animation
import firn_animations as fa
import pathlib
import datetime


number_of_layers_int = 100
start_date = datetime.date(1989,1,1)

#Define dictionary with all the information needed to specify which animation 
#was analysed
simulation_inf = {'matlab_output_folder': '_IterationSite_J_perturbations0', 'initial_density_profile': 'Site_J', 'number_of_layers': '100layers', 'time_period': '1989-2019'}




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



#Save the simulation information (dictionary) as text, for easy change of analysis later on
f = open(path_to_outputs+"simulation_inf.txt","w")
f.write( str(simulation_inf) )
f.close()


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


#Load initial_density profile from data
init_density = pd.read_csv('Input/Initial state/density/'+initial_density_profile+'.csv', sep=',', header='infer')


#Load 2012_density profile from data
density_2012 = pd.read_csv('Input/Initial state/density/RetMIP_density_KAN-U_2012(deleted_NaN).csv', sep=',', header=None)



#Load 2013_density profile from data
density_2013 = pd.read_csv('Input/Initial state/density/DensityProfile_KAN_U_2013_Machguth_et_al.csv', sep=';')


#%%



ani_rho = fa.density_animation(rho_all, Depth_dataframe)
ani_rho.save(path_to_outputs+'Animations/Densityprofile_animation'+name+'.mp4')

plt.clf() 


ani_T = fa.temperature_animation(T_ice_dataframe, Depth_dataframe)
ani_T.save(path_to_outputs+'Animations/Temperature_animation'+name+'.mp4')

plt.clf() 






#Plot initial density_profile from Data
plt.plot(init_density.iloc[:,1], init_density.iloc[:,0])
plt.xlabel(r'Density ($\frac{kg}{m^3}$)')
plt.ylabel(r'Depth ($m$)')
plt.ylim(0, 10)
plt.gca().invert_yaxis()
plt.savefig(path_to_outputs+'plots/Densityprofile1989_data'+name+'.pdf')
plt.close()

#%%
#Plot initial density_profile from Simulation
plt.plot(rho_all.iloc[0], Depth_dataframe.iloc[0])
plt.xlabel(r'Density ($\frac{kg}{m^3}$)')
plt.ylabel(r'Depth ($m$)')
plt.ylim(0, 20)

plt.gca().invert_yaxis()
plt.savefig(path_to_outputs+'plots/Densityprofile1989_simulation'+name+'.pdf')
plt.close()

#%%
#Plot initial density_profile from Simulation and Data sumoultanously
plt.plot(init_density.iloc[:,1], init_density.iloc[:,0], label= 'initial density')
plt.plot(rho_all.iloc[0], Depth_dataframe.iloc[0], label= 'initial density with model layering')
plt.legend()
plt.xlabel(r'Density ($\frac{kg}{m^3}$)')
plt.ylabel(r'Depth ($m$)')
plt.ylim(0, 20)
plt.gca().invert_yaxis()
plt.savefig(path_to_outputs+'plots/Densityprofile1989_simulation_and_data'+name+'.pdf')
plt.close()
#%%






if time_period=='1989-2019':
    #Plot 2012 density_profile from Data
    plt.plot(density_2012.iloc[:,1][:100], density_2012.iloc[:,0][:100], label='2012 Macguth et. al')
    plt.legend()
    plt.xlabel(r'Density ($\frac{kg}{m^3}$)')
    plt.ylim(0, 20)
    plt.gca().invert_yaxis()
    plt.ylabel(r'Depth ($m$)')
    


    plt.savefig(path_to_outputs+'plots/Densityprofile2012_data'+name+'.pdf')
    plt.close()


    #PLot density profile mid-May 2012
    plt.plot(rho_all.iloc[target_hr], Depth_dataframe.iloc[target_hr], label= 'simulation: ' +str(start_date+ target_hr_delta))
    plt.xlabel(r'Density ($\frac{kg}{m^3}$)')
    plt.legend()
    plt.ylim(0, 20)
    plt.gca().invert_yaxis()
    plt.ylabel(r'Depth ($m$)')
    


    plt.savefig(path_to_outputs+'plots/density_profile_may_2012_simulation'+name+'.pdf')
    plt.close()


    #Plot initial density_profile from Simulation and Data sumoultanously
    plt.plot(density_2012[1][:100], density_2012[0][:100], label='2012 Macguth et. al')
    plt.plot(rho_all.iloc[target_hr], Depth_dataframe.iloc[target_hr] , label= 'simulation: ' +str(start_date+ target_hr_delta))
    plt.legend()
    plt.xlabel(r'Density ($\frac{kg}{m^3}$)')
    plt.ylabel(r'Depth ($m$)')
    plt.ylim(0, 10)
    plt.text(400, 30, str(start_date+ target_hr_delta))
    plt.gca().invert_yaxis()
    plt.savefig(path_to_outputs+'plots/Densityprofile2012_simulation_and_data'+name+'_'+str(start_date+ target_hr_delta)+'.pdf')
    plt.close()

    #This block of code generates an array, with differences between the simulation
    #And the experimental data. The density profiles from simulation and data does not 
    #have the same dimensions, and the simkulation have density in terms of layers
    #instead of depth. This Block should take take of this by, creating an-
    #array, by adding data points with same density as the layers it is in.
    diff_dens_2012 = np.zeros(len(density_2012))
    alt_rho_all_2012 = np.zeros(len(density_2012))
    sqr_diffs_2012=0
    j=0
    for i in range(len(diff_dens_2012)):
        cond=0
        while cond==0:
            if i < Depth_dataframe.iloc[target_hr][j]*10 and i < Depth_dataframe.iloc[target_hr][j+1]*10:
                diff_dens_2012[i] = density_2012.iloc[:,1][i] - rho_all.iloc[target_hr][j]
                alt_rho_all_2012[i] = rho_all.iloc[target_hr][j]
                sqr_diffs_2012= (density_2012.iloc[:,1][i] - rho_all.iloc[target_hr][j])**2
                cond = 1
            else:
                j += 1
                
                
                    
        
    
    print(sqr_diffs_2012)
    
    
    
    
    
    
    #PLot to show, that the new density array actually fits the old one
    #plt.plot(alt_rho_all, density_2012[0], '*'), plt.plot(rho_all.iloc[target_hr], Depth_dataframe.iloc[target_hr], '.')
    
    
    
    
    
    #diff_dens_2012 =density_2012[1]- np.ones(len(density_2012))*300
    
    
        
    plt.plot(density_2012[1][:100], density_2012[0][:100], label='2012 Macguth et. al')
    plt.plot(rho_all.iloc[target_hr], Depth_dataframe.iloc[target_hr], label= 'simulation: ' +str(start_date+ target_hr_delta))
    plt.legend()
    plt.plot(diff_dens_2012, density_2012[0], label='difference of the two profiles')
    plt.xlabel(r'Density ($\frac{kg}{m^3}$)')
    plt.ylim(0, 20)
    plt.gca().invert_yaxis()
    plt.ylabel(r'Depth ($m$)')
    
    plt.savefig(path_to_outputs+'plots/Difference_in_Densityprofile2012_simulation_and_data'+name+'.pdf')
    plt.close()
    
    
    
    
    
    plt.plot(np.gradient(density_2012[1]), density_2012[0])
    plt.plot(np.gradient(rho_all.iloc[target_hr]), Depth_dataframe.iloc[target_hr])
    #plt.xlabel(r'Density ($\frac{kg}{m^3}$)')
    plt.ylim(0, 20)
    plt.gca().invert_yaxis()
    plt.ylabel(r'Depth ($m$)')
    
    plt.savefig(path_to_outputs+'plots/Gradients_of_Densityprofile2012_simulation_and_data'+name+'.pdf')
    plt.close()
    


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
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
        #This block of code generates an array, with differences between the simulation
    #And the experimental data. The density profiles from simulation and data does not 
    #have the same dimensions, and the simkulation have density in terms of layers
    #instead of depth. This Block should take take of this by, creating an-
    #array, by adding data points with same density as the layers it is in.
    diff_dens_2013 = np.zeros(len(density_2013))
    alt_rho_all_2013 = np.zeros(len(density_2013))
    sqr_diffs_2013=0
    j=0
    for i in range(len(diff_dens_2013)):
        cond=0
        while cond==0:
            if j>number_of_layers_int:
                cond=1
            elif i < Depth_dataframe.iloc[target_hr+8050][j]*100:
                diff_dens_2013[i] = density_2013.iloc[:,1][i] - rho_all.iloc[target_hr+8050][j]
                alt_rho_all_2013[i] = rho_all.iloc[target_hr+8050][j]
                sqr_diffs_2013= (density_2013.iloc[:,1][i] - rho_all.iloc[target_hr+8050][j])**2
                cond = 1
            else:
                j += 1
    
    
    
    #PLot to show, that the new density array actually fits the old one
    plt.plot(alt_rho_all_2013, density_2013.iloc[:,0], '*', label= 'Interpolation')
    plt.plot(rho_all.iloc[target_hr+8050], Depth_dataframe.iloc[target_hr+8050], '.', label = 'Simulation output')
    plt.legend()
    plt.title('Site J Simulation output and Interpolation')
    plt.ylim(0, 20)
    plt.gca().invert_yaxis()
    plt.xlabel(r'Density ($\frac{kg}{m^3}$)')
    plt.ylabel(r'Depth ($m$)')
    plt.savefig(path_to_outputs+'plots/Comparisonoftoarraywithsamevaluesbutamountofdatapoints.pdf')
         
    
   
    #Same plots but comparing with 2013 core.
            #Plot 2013 density_profile from Data
    plt.plot(density_2013.iloc[:,1], density_2013.iloc[:,0], label='2013 Macguth et. al')
    plt.xlabel(r'Density ($\frac{kg}{m^3}$)')
    plt.ylim(0, 20)
    plt.gca().invert_yaxis()
    plt.legend()
    plt.ylabel(r'Depth ($m$)')
    
    plt.savefig(path_to_outputs+'plots/Densityprofile2013_data'+name+'.pdf')
    plt.close()


    #PLot density profile mid-May 2013
    plt.plot(rho_all.iloc[target_hr+8050], Depth_dataframe.iloc[target_hr+8050], label='simulation '+str(start_date+ target_hr_delta+ datetime.timedelta(hours=8050)))
    plt.xlabel(r'Density ($\frac{kg}{m^3}$)')
    #plt.text(400, 30, str(start_date+ target_hr_delta+ datetime.timedelta(hours=8050)))
    plt.ylim(0, 20)
    plt.gca().invert_yaxis()
    plt.legend()        
    plt.ylabel(r'Depth ($m$)')
    plt.savefig(path_to_outputs+'plots/density_profile_may_2013_simulation'+name+'.pdf')
    plt.close()


    #Plot initial density_profile from Simulation and Data sumoultanously
    plt.plot(density_2013.iloc[:,1], density_2013.iloc[:,0], label='2013 Macguth et. al')
    plt.plot(rho_all.iloc[target_hr+8050], Depth_dataframe.iloc[target_hr+8050], label='simulation '+str(start_date+ target_hr_delta+ datetime.timedelta(hours=8050)))
    plt.xlabel(r'Density ($\frac{kg}{m^3}$)')
    plt.ylabel(r'Depth ($m$)')
    plt.ylim(0, 20)
    #plt.text(400, 30, str(start_date+ target_hr_delta))
    plt.gca().invert_yaxis()
    plt.legend()
    plt.savefig(path_to_outputs+'plots/Densityprofile2013_simulation_and_data'+name+'_'+str(start_date+ target_hr_delta+ datetime.timedelta(hours=8050))+'.pdf')
    plt.close()


    plt.plot(np.gradient(density_2013.iloc[:,1]), density_2013.iloc[:,0])
    plt.plot(np.gradient(alt_rho_all_2013), density_2013.iloc[:,0])
    plt.plot(np.gradient(rho_all.iloc[target_hr+8050]), Depth_dataframe.iloc[target_hr+8050])
    #plt.xlabel(r'Density ($\frac{kg}{m^3}$)')
    plt.ylim(0, 20)
    plt.gca().invert_yaxis()
    plt.ylabel(r'Depth ($m$)')
    plt.savefig(path_to_outputs+'plots/Gradients_of_Densityprofile2013_simulation_and_data'+name+'.pdf')
    plt.close()











    Number_of_negative_gradients_2013={}
    Number_of_negative_gradients_sim={}
    grad_RMSD={}
    
    
    spacing=10
    rho_all_2013_10cm=np.zeros(int(1900/spacing))
    meas_2013_10cm=np.zeros(int(1900/spacing))
    #Compute a meaned array for the density profiles:
    for i in range(0, 1900, spacing):
        rho_all_2013_10cm[int(i/spacing)] = np.mean(alt_rho_all_2013[i:i+spacing])
        meas_2013_10cm[int(i/spacing)] = np.mean(density_2013.iloc[:,1][i:i+spacing])
    print(rho_all_2013_10cm)
 #%%       
    #plt.plot(rho_all.iloc[target_hr+8050], Depth_dataframe.iloc[target_hr+8050], label='simulation '+str(start_date+ target_hr_delta+ datetime.timedelta(hours=8050)))
    #plt.plot(rho_all_2013_10cm, [x*(19.00/len(rho_all_2013_10cm)) for x in range(len(rho_all_2013_10cm))])
    #plt.plot(meas_2013_10cm, [x*(19.00/len(rho_all_2013_10cm)) for x in range(len(rho_all_2013_10cm))])
    
#%% 

    #plt.plot(np.gradient(density_2013.iloc[:,1], density_2013.iloc[:,0]))
    plt.plot(np.gradient(rho_all_2013_10cm, 19.00/len(rho_all_2013_10cm)), [x*(19.00/len(rho_all_2013_10cm)) for x in range(len(rho_all_2013_10cm))], label='simulation '+str(spacing)+'cm spacing', color='orange')
    plt.plot(np.gradient(meas_2013_10cm, 19.00/len(rho_all_2013_10cm)), [x*(19.00/len(rho_all_2013_10cm)) for x in range(len(rho_all_2013_10cm))], label='2013_Machguth_et_al'+str(spacing)+'cm spacing', color='blue')
    #plt.xlabel(r'Density ($\frac{kg}{m^3}$)')
    plt.ylim(0, 20)
    plt.gca().invert_yaxis()
    plt.xlabel(r'Density ($\frac{kg}{m^4}$)')
    plt.ylabel(r'Depth ($m$)')
    
    plt.legend()
    plt.savefig(path_to_outputs+'plots/Gradients_of_Densityprofile2013_simulation_and_data_'+str(spacing)+'cm'+name+'.pdf')
    plt.close()
    grad_RMSD.update({'spacing='+str(spacing)+'cm': np.sum((np.gradient(rho_all_2013_10cm)-np.gradient(meas_2013_10cm))**2)})
    Number_of_negative_gradients_2013.update({'spacing='+str(spacing)+'cm': sum(np.gradient(meas_2013_10cm)<0)})
    Number_of_negative_gradients_sim.update({'spacing='+str(spacing)+'cm': sum(np.gradient(rho_all_2013_10cm)<0)})


    
    spacing=20
    rho_all_2013_20cm=np.zeros(int(1900/spacing))
    meas_2013_20cm=np.zeros(int(1900/spacing))
    #Compute a meaned array for the density profiles:
    for i in range(0, 1900, spacing):
        rho_all_2013_20cm[int(i/spacing)] = np.mean(alt_rho_all_2013[i:i+spacing])
        meas_2013_20cm[int(i/spacing)] = np.mean(density_2013.iloc[:,1][i:i+spacing])
    print(rho_all_2013_20cm)


    #plt.plot(np.gradient(density_2013.iloc[:,1], density_2013.iloc[:,0]))
    plt.plot(np.gradient(rho_all_2013_20cm, 19.00/len(rho_all_2013_20cm)), [x*(19.00/len(rho_all_2013_20cm)) for x in range(len(rho_all_2013_20cm))], label='simulation '+str(spacing)+'cm spacing', color='orange')
    plt.plot(np.gradient(meas_2013_20cm, 19.00/len(rho_all_2013_20cm)), [x*(19.00/len(rho_all_2013_20cm)) for x in range(len(rho_all_2013_20cm))], label='2013_Machguth_et_al'+str(spacing)+'cm spacing', color='blue')
    #plt.xlabel(r'Density ($\frac{kg}{m^3}$)')
    plt.ylim(0, 20)
    plt.gca().invert_yaxis()
    plt.xlabel(r'Density ($\frac{kg}{m^4}$)')
    plt.ylabel(r'Depth ($m$)')
    
    plt.legend()
    
    plt.savefig(path_to_outputs+'plots/Gradients_of_Densityprofile2013_simulation_and_data_'+str(spacing)+'cm'+name+'.pdf')
    plt.close()
    grad_RMSD.update({'spacing='+str(spacing)+'cm': np.sum((np.gradient(rho_all_2013_20cm)-np.gradient(meas_2013_20cm))**2)})
    Number_of_negative_gradients_2013.update({'spacing='+str(spacing)+'cm': sum(np.gradient(meas_2013_20cm)<0)})
    Number_of_negative_gradients_sim.update({'spacing='+str(spacing)+'cm': sum(np.gradient(rho_all_2013_20cm)<0)})   
    
    
    
    
    
    
    spacing=50
    rho_all_2013_50cm=np.zeros(int(1900/spacing))
    meas_2013_50cm=np.zeros(int(1900/spacing))
    #Compute a meaned array for the density profiles:
    for i in range(0, 1900, spacing):
        rho_all_2013_50cm[int(i/spacing)] = np.mean(alt_rho_all_2013[i:i+spacing])
        meas_2013_50cm[int(i/spacing)] = np.mean(density_2013.iloc[:,1][i:i+spacing])
    print(rho_all_2013_50cm)


    #plt.plot(np.gradient(density_2013.iloc[:,1], density_2013.iloc[:,0]))
    plt.plot(np.gradient(rho_all_2013_50cm, 19.00/len(rho_all_2013_50cm)), [x*(19.00/len(rho_all_2013_50cm)) for x in range(len(rho_all_2013_50cm))], label='simulation '+str(spacing)+'cm spacing', color='orange')
    plt.plot(np.gradient(meas_2013_50cm, 19.00/len(rho_all_2013_50cm)), [x*(19.00/len(rho_all_2013_50cm)) for x in range(len(rho_all_2013_50cm))], label='2013_Machguth_et_al'+str(spacing)+'cm spacing', color='blue')
    #plt.xlabel(r'Density ($\frac{kg}{m^3}$)')
    plt.ylim(0, 20)
    plt.gca().invert_yaxis()
    plt.xlabel(r'Density ($\frac{kg}{m^4}$)')
    plt.ylabel(r'Depth ($m$)')
    
    plt.legend()
    
    plt.savefig(path_to_outputs+'plots/Gradients_of_Densityprofile2013_simulation_and_data_'+str(spacing)+'cm'+name+'.pdf')
    plt.close()
    grad_RMSD.update({'spacing='+str(spacing)+'cm': np.sum((np.gradient(rho_all_2013_50cm)-np.gradient(meas_2013_50cm))**2)})
    Number_of_negative_gradients_2013.update({'spacing='+str(spacing)+'cm': sum(np.gradient(meas_2013_50cm)<0)})
    Number_of_negative_gradients_sim.update({'spacing='+str(spacing)+'cm': sum(np.gradient(rho_all_2013_50cm)<0)})

    spacing=100
    rho_all_2013_100cm=np.zeros(int(1900/spacing))
    meas_2013_100cm=np.zeros(int(1900/spacing))
    #Compute a meaned array for the density profiles:
    for i in range(0, 1900, spacing):
        rho_all_2013_100cm[int(i/spacing)] = np.mean(alt_rho_all_2013[i:i+spacing])
        meas_2013_100cm[int(i/spacing)] = np.mean(density_2013.iloc[:,1][i:i+spacing])


    #plt.plot(np.gradient(density_2013.iloc[:,1], density_2013.iloc[:,0]))
    plt.plot(np.gradient(rho_all_2013_100cm, 19.00/len(rho_all_2013_100cm)), [x*(19.00/len(rho_all_2013_100cm)) for x in range(len(rho_all_2013_100cm))], label='simulation '+str(spacing)+'cm spacing', color='orange')
    plt.plot(np.gradient(meas_2013_100cm, 19.00/len(rho_all_2013_100cm)), [x*(19.00/len(rho_all_2013_100cm)) for x in range(len(rho_all_2013_100cm))], label='2013_Machguth_et_al'+str(spacing)+'cm spacing', color='blue')
    #plt.xlabel(r'Density ($\frac{kg}{m^3}$)')
    plt.ylim(0, 20)
    plt.gca().invert_yaxis()
    plt.xlabel(r'Density ($\frac{kg}{m^4}$)')
    plt.ylabel(r'Depth ($m$)')
    
    plt.legend()
    
    plt.savefig(path_to_outputs+'plots/Gradients_of_Densityprofile2013_simulation_and_data_'+str(spacing)+'cm'+name+'.pdf')
    plt.close()
    grad_RMSD.update({'spacing='+str(spacing)+'cm': np.sum((np.gradient(rho_all_2013_100cm)-np.gradient(meas_2013_100cm))**2)})
    Number_of_negative_gradients_2013.update({'spacing='+str(spacing)+'cm': sum(np.gradient(meas_2013_100cm)<0)})
    Number_of_negative_gradients_sim.update({'spacing='+str(spacing)+'cm': sum(np.gradient(rho_all_2013_100cm)<0)})
    
    
    
    
    spacing=190
    rho_all_2013_190cm=np.zeros(int(1900/spacing))
    meas_2013_190cm=np.zeros(int(1900/spacing))
    #Compute a meaned array for the density profiles:
    for i in range(0, 1900, spacing):
        rho_all_2013_190cm[int(i/spacing)] = np.mean(alt_rho_all_2013[i:i+spacing])
        meas_2013_190cm[int(i/spacing)] = np.mean(density_2013.iloc[:,1][i:i+spacing])


    #plt.plot(np.gradient(density_2013.iloc[:,1], density_2013.iloc[:,0]))
    plt.plot(np.gradient(rho_all_2013_190cm, 19.00/len(rho_all_2013_190cm)), [x*(19.00/len(rho_all_2013_190cm)) for x in range(len(rho_all_2013_190cm))], label='simulation '+str(spacing)+'cm spacing', color='orange')
    plt.plot(np.gradient(meas_2013_190cm, 19.00/len(rho_all_2013_190cm)), [x*(19.00/len(rho_all_2013_190cm)) for x in range(len(rho_all_2013_190cm))], label='2013_Machguth_et_al'+str(spacing)+'cm spacing', color='blue')
    #plt.xlabel(r'Density ($\frac{kg}{m^3}$)')
    plt.xlabel(r'Density ($\frac{kg}{m^4}$)')
    plt.ylabel(r'Depth ($m$)')
    plt.ylim(0, 20)
    plt.gca().invert_yaxis()

    plt.legend()
    
    plt.savefig(path_to_outputs+'plots/Gradients_of_Densityprofile2013_simulation_and_data_'+str(spacing)+'cm'+name+'.pdf')
    plt.close()
    grad_RMSD.update({'spacing='+str(spacing)+'cm': np.sum((np.gradient(rho_all_2013_190cm)-np.gradient(meas_2013_190cm))**2)})
    Number_of_negative_gradients_2013.update({'spacing='+str(spacing)+'cm': sum(np.gradient(meas_2013_190cm)<0)})
    Number_of_negative_gradients_sim.update({'spacing='+str(spacing)+'cm': sum(np.gradient(rho_all_2013_190cm)<0)})
print(grad_RMSD)
print(grad_RMSD)
print(Number_of_negative_gradients_2013)
print(Number_of_negative_gradients_sim)
f = open(path_to_outputs+"Number_of_negative_gradients_2013.txt","w")
f.write( str(Number_of_negative_gradients_2013) )
f.close()
f = open(path_to_outputs+"Number_of_negative_gradients_sim.txt","w")
f.write( str(Number_of_negative_gradients_sim) )
f.close()











