# Project-2020-GEUS


GEUS SEB firn model density profile perturbers and data analysis. 
-----------------------------------------------------------------------------------------------------------------------------------

Project by Tom Edwards and Uffe Bundesen of roskilde university

edwards@ruc.dk
ufbu@ruc.dk

-----------------------------------------------------------------------------------------------------------------------------------
MATLAB scripts - Tom Edwards
-----------------------------------------------------------------------------------------------------------------------------------
The main body of Matlab work is in the script Allinone_new.m

The script is annotated but in brief does the following:

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

-----------------------------------------------------------------------------------------------------------------------------------

-----------------------------------------------------------------------------------------------------------------------------------
Python scripts - Uffe Bundesen
-----------------------------------------------------------------------------------------------------------------------------------

The two main scripts are GradientsasfunctionoftimeLOOP.py and Readoutput.py.

In both of these scripts, you need to specify the outputfolder from the SEB-firn model, the initial density profile. Also the number of layers and the time period the model will run for (weather data). 

In GradientsasfunctionoftimeLOOP.py you can ran through more outputs from the SEB firn model, in a for loop. THis was done for the homogenous and the thick pertubations.
The outputs are placed in a folder PythonOutputs with a subfolder structure /densityprofile/timeperiod/layers.


The GradientsasfunctionoftimeLOOP creates Delta rho plots
as well as NGN plots for a given perturbation. Lastly the script plots Delta rho and mean in the region (1,15} for the full time period for three different pertubations.



Gfunct is Called from GradientsasfunctionoftimeLOOP. This script does the interpolation of the density profile, and returns a evenly spaced density profile at a 10cm resolution.
do the interpolation.


The Readoutput script, does all the animations
and some plots for a given model
Readoutput. Among other plots are the 
inital profile, 
the initial together with the interpolation, 
the profile in 2012 compared with measured.
the profile in 2013 compared with measured.
And some plots related to the gradient of the density profile.
In order to run it, you need to specify the outputfolder
from the SEB-firn model, the initial density profile, the numbeer of layers
and the time period the model has been running. The script also creates a text file with the simulation information, so the script can easily be run again, if you want to change the plots. Readoutput places the plots it creates in a subfolder /densityprofile/timeperiod/layers/plots


The Code is presented as is. Since the scripts are not everywhere commented thouroughly feel free to ask, if any questions might arise.


Firn Animations is Called from Readoutput to create the animations. Two animations can be created. One for the density profile and the other for the temperature profile, over the time period of the simulation.




Licensed under BY CC-SA
