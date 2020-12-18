# Project-2020-GEUS
GEUS SEB firn model density profile perturbers and data analysis. 

The main body of Matlab work is in the script Allinone_new.m

The script does the following:

- Imports the netCDF files from a simulation run in the SEB firn model
- Caclaultes the total density for each layer
- Extracts the denisty profile at the relevant year
- Imports Average density values for each region and resolution of the 2012 and 2013 DP's
- Discretises and interpolates the profile into 10cm intervals (also see ConvertDepth.m)
- Calculates the RMSD over different resolutions and regions between the final profile and the 2012/2013 profiles
- Calculates the RMSD over different resolutions and regions between the final profile and the intial profile
- Plots of perturbed profiles, initial profiles and simulated profiles
- Distribution of layers as a function of depths
- interpolation plots for a visual check on how distribution changes
- Perturbation checker

Other important scripts: 

RMSD_regions_2012 and 2013

- These scripts plot the RMSD values for each region and at different resolutions for each perturbation

RMSD_totals

- plots the RMSD for the whole of the top 10m of the density profile as opposed to each region

AVGD_calc

- Calculates the average density over different resolutions for each region in the profile

AVGDENS_linearly_distributed densities.m

- Density perturbing script

Ice_lens_generator for 2012 and 1989 profiles

- Generates and distributes ice lenses in a profile randomly

RandomDP generator

- Generates random Density profile by mixing up values in an existing one





Licensed under BY CC-SA
