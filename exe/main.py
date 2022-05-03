#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# main.py
#
# Main body of the analysis. Calls functions from lib to
# perform specific analysis
#
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

def run_analysis(data_dir, out_label, sf_zint_log, WG_bounds, sf_xint_log, sf_xint_interp_log,
                nn_rhop, tvar_window, eke_log, tracer_xint_log, 
                xmin_list, xmax_list, range_labels):
    """
    Run da code. Waaah
    """
    import iris
    import numpy as np
    import sys
    import os

    sys.path.append(os.path.abspath("../lib"))

    import cubeprep
    import eddy_energy
    import streamfunction
    import tracers
    
    out_dir = data_dir + '/OUTPUT.' + out_label + '/'
    if not os.path.exists(out_dir): os.mkdir(out_dir)

    data_list = iris.load(data_dir + "/*grid*.nc")
    mask_list = iris.load(data_dir + "/mesh_mask.nc")

    print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
    print("run_analysis")
    print("Outputs will be found at: ", out_dir)
    print("Carrying out the following operations")
    print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")

    if sf_zint_log == True: 
        print("")
        print("Stream function and transport calculation")
        print("Bounds for Weddell Gyre transport (x_min, x_max, y_min, y_max):", WG_bounds)

    if sf_xint_log == True:
        print("")
        print("Residual overturning stream function and transport calculation")

    if sf_xint_interp_log == True:
        print("")
        print("Interpolate overturning stream function onto depth space")

        if sf_xint_log == False: print("Will load res_ov.nc file")

    if eke_log == True: 
        print("")
        print("Eddy kinetic energy")
        print("Rolling window length = ", tvar_window, " days")

    if tracer_xint_log == True: 
        print("")
        print("Zonal integrals of tracers")
        for i in range(len(xmin_list)):
            print(" Range: ", range_labels[i] )
            print(" xmin = ", str(xmin_list[i]))
            print(" xmax = ", str(xmax_list[i]))

    #Load the dictionary of variable names for this data
    varname_dict = cubeprep.var_names()

    #Stream function calculations
    if sf_zint_log == True:
        print("")
        print("Calculating stream function of the depth-integrated flow >>>")
        sf_zint_cube, WG_transport_cube, ACC_transport_cube = streamfunction.DepIntSf(data_list, mask_list, varname_dict, WG_bounds=WG_bounds)

        iris.save(sf_zint_cube, out_dir + '/sf_zint.nc')
        iris.save(WG_transport_cube, out_dir + '/WG_transport.nc')
        iris.save(ACC_transport_cube, out_dir + '/ACC_transport.nc')

    if sf_xint_log == True:
        print("")
        print("Calculating stream function of the x-integrated flow (Overturning) >>>")
        res_ov_cube = streamfunction.ResidualOverturning(data_list, mask_list, nn_rhop, varname_dict)
        iris.save(res_ov_cube, out_dir + '/res_ov.nc')

    if sf_xint_interp_log == True:
        print("")
        print("Interpolating the residual overturning stream function onto depth space")
        if sf_xint_log != True:
            try:
                res_ov_cube = iris.load( out_dir + '/res_ov.nc')[0]   
            except:
                print("sf_xint_log == False and cannot load res_ov_cube")
                print("Skipping interpolation of residual overturning stream function !")

        try:
            res_ov_depth_cube = streamfunction.ResOv2depth(res_ov_cube, data_list, mask_list, varname_dict)
            iris.save(res_ov_depth_cube, out_dir + '/res_ov_depth.nc')
        except: 
            print("Failed to convert to depth coordinates")


    if eke_log == True:
        print("")
        print("Calculating the Eddy Kinetic Energy >>>")
        eke_cube, eke_zint_cube = eddy_energy.eddy_kinetic_energy(data_list, mask_list, tvar_window, varname_dict)

        iris.save(eke_cube, out_dir + "/eke." + str(tvar_window) +".nc")
        iris.save(eke_zint_cube, out_dir + "/eke_zint." + str(tvar_window) +".nc")
        print("Complete -------------------------")

    #Calculate the zonal mean of tracers over specified ranges
    if tracer_xint_log == True:
        print("")
        print("Calculating the zonal mean of tracers over specified ranges")

        temp_xint_cube_list, sal_xint_cube_list, rhop_xint_cube_list = tracers.tracer_xint(data_list, mask_list, varname_dict, xmin_list, xmax_list, range_labels)

        print(temp_xint_cube_list)

        iris.save(temp_xint_cube_list, out_dir + "/temp_xint.nc")
        iris.save(sal_xint_cube_list, out_dir + "/sal_xint.nc")
        iris.save(rhop_xint_cube_list, out_dir + "/rhop_xint.nc")

    return