#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# main.py
#
# Main body of the analysis. Calls functions from lib to
# perform specific analysis
#
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

def run_analysis(data_dir, out_label, sf_zint_log, WG_bounds, sf_xint_log, sf_xint_interp_log,
                nn_rhop, nn_z, sponge_sample_dict, tvar_window, eke_log, tracer_xint_log,
                vel_xint_log, eke_xint_log, xmin_list, xmax_list, range_labels, ACC_decomp_log,
                WG_decomp_log):
    """
    Runs the program from settings in executable.py
    """
    import dask
    import dask.array as da
    import xarray as xr
    import numpy as np
    import sys
    import os
    sys.path.append(os.path.abspath("../lib"))
    import cubeprep
    import eddy_energy
    import streamfunction
    import tracers
    
    #Define output directory and create it if needed
    out_dir = data_dir + '/OUTPUT.' + out_label + '/'
    if not os.path.exists(out_dir): os.mkdir(out_dir)

    #Open the model output file and the mesh_mask file
    data_list = xr.open_mfdataset(data_dir + "/CANAL_grid_*.nc", chunks={"time_counter":1, "x":750, "y":750}, compat='override')
    mask_list = xr.open_dataset(data_dir + "/mesh_mask.nc", chunks={"x":750, "y":750})

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

    if vel_xint_log == True: 
        print("")
        print("Zonal integrals of x and y velocities")
        for i in range(len(xmin_list)):
            print(" Range: ", range_labels[i] )
            print(" xmin = ", str(xmin_list[i]))
            print(" xmax = ", str(xmax_list[i]))

    if ACC_decomp_log == True:
        print("")
        print("Decomposition of ACC transport")

    if WG_decomp_log == True:
        print("")
        print("Decomposition of Weddell Gyre stream function")


    #Load the dictionary of variable names for this data
    varname_dict = cubeprep.var_names()

    #Stream function calculations
    if sf_zint_log == True:
        print("")
        print("Calculating stream function of the depth-integrated flow >>>")
        sf_zint_cube, WG_transport_cube, ACC_transport_cube = streamfunction.DepIntSf(data_list, mask_list, varname_dict, WG_bounds=WG_bounds)

        sf_zint_cube.to_netcdf( out_dir + '/sf_zint.nc')
        WG_transport_cube.to_netcdf( out_dir + '/WG_transport.nc')
        ACC_transport_cube.to_netcdf( out_dir + '/ACC_transport.nc')

    if sf_xint_log == True:
        print("")
        print("Calculating stream function of the x-integrated flow (Overturning) >>>")
        res_ov_cube, rhop_depth_cube = streamfunction.ResidualOverturning(data_list, mask_list, nn_rhop, sponge_sample_dict, varname_dict)
        xr.merge([res_ov_cube, rhop_depth_cube]).to_netcdf(out_dir + '/res_ov.nc')

    if sf_xint_interp_log == True:
        print("")
        print("Interpolating the residual overturning stream function onto depth space")
        if sf_xint_log != True:
            try: 
                res_ov_cube = xr.open_dataset(out_dir + "/res_ov.nc")['res_ov']
                rhop_depth_cube = xr.open_dataset(out_dir + "/res_ov.nc")['rhop_depth']
            except:
                print("sf_xint_log == False and cannot load res_ov_cube")
                print("Skipping interpolation of residual overturning stream function !")

        # try:
        res_ov_depth_cube = streamfunction.ResOv2depth(res_ov_cube, rhop_depth_cube, data_list, mask_list, varname_dict, nn_z=nn_z)
        res_ov_depth_cube.to_netcdf(out_dir + '/res_ov_depth.nc')

        # except: 
        # print("Failed to convert to depth coordinates")


    if eke_log == True:
        print("")
        print("Calculating the Eddy Kinetic Energy >>>")
        eke_cube, eke_zint_cube, eke_zmean_cube = eddy_energy.eddy_kinetic_energy(data_list, mask_list, varname_dict)

        eke_cube.to_netcdf(out_dir + "/eke.nc")
        eke_zint_cube.to_netcdf(out_dir + "/eke_zint.nc")
        eke_zmean_cube.to_netcdf(out_dir + "/eke_zmean.nc")
        print("Complete -------------------------")

    #Calculate the zonal mean of tracers over specified ranges
    if tracer_xint_log == True:
        print("")
        print("Calculating the zonal mean of tracers over specified ranges")

        temp_xint_cube_list, sal_xint_cube_list, rhop_xint_cube_list = tracers.tracer_xint(data_list, mask_list, varname_dict, xmin_list, xmax_list, range_labels)

        xr.merge(temp_xint_cube_list).to_netcdf(out_dir + "/temp_xint.nc")
        xr.merge(sal_xint_cube_list).to_netcdf(out_dir + "/sal_xint.nc")
        xr.merge(rhop_xint_cube_list).to_netcdf(out_dir + "/rhop_xint.nc")
        
    if vel_xint_log == True:
        print("")
        print("Calculating the zonal mean of x and y velocities over specified ranges")

        u_xint_cube_list = tracers.uvar_xint(data_list, mask_list, varname_dict, xmin_list, xmax_list, range_labels)
        xr.merge(u_xint_cube_list).to_netcdf(out_dir + "/u_xint.nc")

        v_xint_cube_list = tracers.vvar_xint(data_list, mask_list, varname_dict, xmin_list, xmax_list, range_labels)
        xr.merge(v_xint_cube_list).to_netcdf(out_dir + "/v_xint.nc")

    if eke_xint_log == True:
        print("")
        print("Calculating the zonal mean of eddy kinetic energy over specified ranges")
        if eke_log == False:
            try:
                eke_cube = xr.open_dataset(out_dir + "/eke.nc")['eke']
            except: 
                "Unable to load eddy kinetic energy data"

        eke_xint_cube_list = tracers.eke_xint(eke_cube, mask_list, varname_dict, xmin_list, xmax_list, range_labels)
        xr.merge(eke_xint_cube_list).to_netcdf(out_dir + "/eke_xint.nc")

    if ACC_decomp_log == True:
        print("")
        print("Decomposing the ACC transport into barotropic and baroclinic parts")

        acc_decomp_cube_list = streamfunction.acc_decomp(data_list, mask_list, varname_dict)
        xr.merge(acc_decomp_cube_list).to_netcdf(out_dir + "/acc_decomp.nc")

    if WG_decomp_log == True:
        print("")
        print("Decomposing the WG transport into compressible and incompressible parts")
        print("... also into barotropic and baroclinic parts")

        if ACC_decomp_log != True:
            try:
                acc_decomp_cube_list = xr.open_dataset(out_dir + "/acc_decomp.nc")
            except:
                print("ACC_decomp_log == False and cannot load acc_decomp.nc")
                print("Skipping decomposition of WG")
        else:
            acc_decomp_cube_list = xr.merge(acc_decomp_cube_list)

        if sf_zint_log != True:
            try:
              sf_zint_cube = xr.open_dataset(out_dir + "/sf_zint.nc")
            except:
              print("sf_zint_cube == False and cannot load sf_zint.nc")
              print("Skipping decomposition of WG")


        phi_cube_list = streamfunction.WG_decomp( data_list, acc_decomp_cube_list,
                                                    mask_list, varname_dict )

        phi_cube_list = xr.merge(phi_cube_list)
        phi_cube_list.to_netcdf(out_dir + '/phi.nc')

        u_corr_list = streamfunction.u_corr(  data_list, acc_decomp_cube_list, phi_cube_list, sf_zint_cube, WG_bounds,
                                              mask_list, varname_dict )

        u_corr_list = xr.merge(u_corr_list)
        u_corr_list.to_netcdf(out_dir + '/u_corr.nc')

    return
