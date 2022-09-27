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
out_label = "ucompr_test"

#Calculate the stream function and transports
sf_zint_log = False
WG_bounds = (0,9999999,-9999999, 0) #bounds for the Weddell Gyre transport calculation in km (x_min, x_max, y_min, y_max) [tuple]
                 # = None for no bounds

#Calculate the residual overturning stream function
sf_xint_log = False
nn_rhop = 11    # Maximum number of isopycnal levels
sf_xint_interp_log = False
nn_z = 31        # Number of z levels to interpolate to

#Calculate eddy energies and isopycnal variability
eke_log = False     # Eddy kinetic energy
tvar_window = None  # Time-averaging window for eddy energy calculation in days
                      # = None to use full time window

#Calculate full and partial zonal means of temperature, salinity, and density
tracer_xint_log = False

#Calculate full and partial zonal means of x and y velocities
vel_xint_log = False

#Calculate full and partial zonal means of eddy kinetic energy
eke_xint_log = False

#Specify ranges for tracer_xint and vel_xint
xmin_list = (3500, 500, 2000, -2000)
xmax_list = (4000, 1000, 2500, -1000)
range_labels = ('noridge', 'ridgewest', 'ridgeeast', 'channel')

#ACC decomp
ACC_decomp_log = False

#Separate the time-averaged flow into compressible and incompressible parts using an elliptical solver
#Use this as a correction for all calculations of stream function
WG_decomp_log = False

main.run_analysis(data_dir, out_label, sf_zint_log, WG_bounds, sf_xint_log, sf_xint_interp_log, nn_rhop, nn_z, tvar_window, eke_log, tracer_xint_log, vel_xint_log, eke_xint_log,
                  xmin_list, xmax_list, range_labels, ACC_decomp_log, WG_decomp_log)
