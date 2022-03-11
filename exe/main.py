#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# main.py
#
# Main body of the analysis. Calls functions from lib to
# perform specific analysis
#
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

def run_analysis(data_dir, out_label, sf_zint_log, WG_bounds, sf_xint_log, rel_vort_log,
                z_rhop_log, nn_rhop, tvar_window, eke_log, epe_log):
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
    
    out_dir = data_dir + '/OUTPUT.' + out_label + '/'
    if not os.path.exists(out_dir): os.mkdir(out_dir)

    data_list = iris.load(data_dir + "/*grid*.nc")
    mask_list = iris.load(data_dir + "/mesh_mask.nc")

    print(mask_list)

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
        print("Overturning stream function and transport calculation")


    # if rel_vort_log == True: 
    #     print("")
    #     print("Relative vorticity")

    if z_rhop_log == True: 
        print("")
        print("Interpolated depth of the isopycnals")
        print("Maximum number of isopycnal levels =", nn_rhop)

    if eke_log == True: 
        print("")
        print("Eddy kinetic energy")
        print("Rolling window length = ", tvar_window, " days")

    if epe_log == True: 
        print("")
        print("Eddy potential energy")
        print("Rolling window length = ", tvar_window, " days")

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
        # sf_xint_cube = streamfunction.ZonIntSF(data_list, mask_list, varname_dict)

        # iris.save(sf_xint_cube, out_dir + '/sf_xint.nc')
        res_ov_cube = streamfunction.ResidualOverturning(data_list, mask_list, nn_rhop, varname_dict)
        print(res_ov_cube)
        # print(v_rhop_cube)
        iris.save(res_ov_cube, out_dir + '/res_ov.nc')
        # iris.save(v_rhop_cube, out_dir + '/v_rhop.nc')




    # Eddy energy calculations
    if z_rhop_log == True:
        print("")
        print("Calculating depth of isopycnals>>>")
        z_rhop_cube = eddy_energy.z_rhop(data_list, mask_list, nn_rhop, varname_dict)  

        iris.save(z_rhop_cube, out_dir + "/z_rhop2.nc")
        print("Complete -------------------------")


    if eke_log == True:
        print("")
        print("Calculating the Eddy Kinetic Energy >>>")
        eke_cube, eke_zint_cube = eddy_energy.eddy_kinetic_energy(data_list, mask_list, tvar_window, varname_dict)

        iris.save(eke_cube, out_dir + "/eke." + str(tvar_window) +".nc")
        iris.save(eke_zint_cube, out_dir + "/eke_zint." + str(tvar_window) +".nc")
        print("Complete -------------------------")

    if epe_log == True:
        print("")
        print("Calculating the Eddy Kinetic Energy >>>")

        if z_rhop_log == False: z_rhop_cube = iris.load(out_dir + "/z_rhop.nc")[0]

        epe_cube, epe_rhoint_cube = eddy_energy.eddy_potential_energy(z_rhop_cube, data_list, mask_list, tvar_window, varname_dict)

        iris.save(epe_cube, out_dir + "/epe." + str(tvar_window) +".nc")
        iris.save(epe_rhoint_cube, out_dir + "/epe_rhoint." + str(tvar_window) +".nc")
        print("Complete -------------------------")


    return