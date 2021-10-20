# Plotting routines and processed data for ACCESS-CM2 Specific Heat Calculations

This repository contains post-processed ACCESS-CM2 PI control CMIP6 climate model output and code used to produce Figs. 5, 6 and 9 in the published article:

McDougall, T., J., Barker, P.M.,Holmes, R.M., Pawlowicz, R., Griffies, S. and Durack, P. (2021): The interpretation of temperature and salinity variables in numerical ocean model output and the calculation of heat fluxes and heat content,Geoscientific Model Development,14, 1-21, [https://doi.org/10.5194/gmd-14-1-2021](https://doi.org/10.5194/gmd-14-1-2021)

The processed ACCESS-CM2 data is included as .mat files and is accompanied by Matlab processing routines (including code from the TEOS-10 Gibbs SeaWater Oceanographic Toolbox, [https://www.teos-10.org/software.htm#1](https://www.teos-10.org/software.htm#1)) to produce the figures.

ACCESS_SpecificHeat_PIcontrol_SWP.mat contains processed data from the ACCESS-CM2 PI control run. 

ACCESS_SpecificHeat.m contains the plotting code for Figures 5 & 6.

ACCESS_SpecificHeat_PIcontrol_CTptERROR.mat contains processed data for the plotting of Figure 9.

ACCESS_SpecificHeat_CTptERROR.m contains the plotting code for Figure 9.

While most of the raw data is freely available on the [CMIP6 ESGF website](https://esgf-node.llnl.gov/projects/cmip6/), the hfds surface heat flux variable for ACCESS-CM2 in that archive is not complete. Instead, a complete version of the raw ACCESS-CM2 output is used that is available on the National Computational Institute.
