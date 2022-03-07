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

sf_xint_log = True

#Calculate the relative vorticity field and vorticity of the depth-integrated velocity
rel_vort_log = False


#Calculate eddy energies and isopycnal variability
z_rhop_log = True # Depth of isopycnals
nn_rhop = 31      # Maximum number of isopycnal levels
eke_log = True  # Eddy kinetic energy
epe_log = True    # Eddy available potential energy (requires a calculation of z_rhop)
tvar_window = None  # Time-averaging window for eddy energy calculation in days
                  # = None to use full time window

main.run_analysis(data_dir, out_label, sf_zint_log, WG_bounds, sf_xint_log, rel_vort_log, z_rhop_log, nn_rhop, tvar_window, eke_log, epe_log)


