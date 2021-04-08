# Plotting routines and processed data for ACCESS-CM2 Specific Heat Calculations

Contains data and code to plot specific heat properties from ACCESS-CM2 PI control output included as Figs. 5, 6 and 9 in the submitted article:

McDougall, T., J., Barker, P.M., Holmes, R.M., Pawlowicz, R., Griffies, S. and Durack, P.: The interpretation of temperature and salinity variables in numerical ocean model output, and the calculation of heat fluxes and heat content. Submitted to the Geoscientific Model Development.

ACCESS_SpecificHeat.m contains the plotting code for Figures 5 & 6.

ACCESS_SpecificHeat_PIcontrol_SWP.mat contains processed data from the ACCESS-CM2 PI control run. 

ACCESS_SpecificHeat_CTptERROR.m contains the plotting code for Figure 9.

ACCESS_SpecificHeat_PIcontrol_CTptERROR.mat contains processed data for the plotting of Figure 9.

While most of the raw data is freely available on the [CMIP6 ESGF website](https://esgf-node.llnl.gov/projects/cmip6/), the hfds surface heat flux variable for ACCESS-CM2 in that archive is not complete. Here we use a complete version from the raw ACCESS-CM2 output available on the National Computational Institute.

If you have any questions, please contact me at ryan.holmes@unsw.edu.au
