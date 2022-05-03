#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# executable.py
#
# Contains options for the user and executes the main script
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

import sys
import os
import main

data_dir = sys.argv[1]


#Input settings
out_label = "test"

#Calculate the stream function and transports
sf_zint_log = True
WG_bounds = (0,9999999,-9999999, 0) #bounds for the Weddell Gyre transport calculation in km (x_min, x_max, y_min, y_max) [tuple]
                 # = None for no bounds

#Calculate the residual overturning stream function
sf_xint_log = True
sf_xint_interp_log = True
nn_rhop = 101      # Maximum number of isopycnal levels

#Calculate eddy energies and isopycnal variability
eke_log = True   # Eddy kinetic energy
tvar_window = 'year'  # Time-averaging window for eddy energy calculation in days
                      # = None to use full time window

#Calculate full and partial zonal integrals of temperature, salinity, and density
tracer_xint_log = True
xmin_list = (3500, 500)
xmax_list = (3600, 600)
range_labels = ('range1', 'range2')

main.run_analysis(data_dir, out_label, sf_zint_log, WG_bounds, sf_xint_log, sf_xint_interp_log, nn_rhop, tvar_window, eke_log, tracer_xint_log, xmin_list, xmax_list, range_labels)


