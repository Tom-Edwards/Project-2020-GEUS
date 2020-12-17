#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 21:25:26 2020

@author: ufbu
"""

import numpy as np
from numba import jit

number_of_layers_int=100

@jit
def even_array(rho_array, depth_array, output_depth=2000):
        #This block of code generates an array, with differences between the simulation
    #And the experimental data. The density profiles from simulation and data does not 
    #have the same dimensions, and the simkulation have density in terms of layers
    #instead of depth. This Block should take take of this by, creating an-
    #array, by adding data points with same density as the layers it is in.
    alt_rho_all = np.zeros(output_depth)
    j=0
    for i in range(output_depth):
        cond=0
        while cond==0:
            if j>number_of_layers_int:
                cond=1
            elif i < depth_array[j]*100:
                #diff_dens_2013[i] = density_2013.iloc[:,1][i] - rho_all.iloc[k][j]
                alt_rho_all[i] = rho_array[j]
                #sqr_diffs_2013= (density_2013.iloc[:,1][i] - rho_all.iloc[k][j])**2
                cond = 1
            else:
                j += 1
    return (alt_rho_all)


@jit
def NGN(alt_rho_all, spacing, output_depth_cm):
    rho_all_spacing=np.zeros(int(output_depth_cm/spacing))
    for i in range(0,output_depth_cm,spacing):
        rho_all_spacing[int(i/spacing)] = np.mean(alt_rho_all[i:i+spacing])
    return rho_all_spacing