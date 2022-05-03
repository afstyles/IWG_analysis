"""
          _     _                                                        
         | |   | |                                                       
  ___  __| | __| |_   _     ___ _ __   ___ _ __ __ _ _   _   _ __  _   _ 
 / _ \/ _` |/ _` | | | |   / _ \ '_ \ / _ \ '__/ _` | | | | | '_ \| | | |
|  __/ (_| | (_| | |_| |  |  __/ | | |  __/ | | (_| | |_| |_| |_) | |_| |
 \___|\__,_|\__,_|\__, |   \___|_| |_|\___|_|  \__, |\__, (_) .__/ \__, |
                   __/ |_____                   __/ | __/ | | |     __/ |
                  |___/______|                 |___/ |___/  |_|    |___/ 
eddy_energy.py

Calculations of eddy energies

Contains methods:
    eddy_kinetic_energy - Calculates the eddy kinetic energy over a specified time window
    im1 --> Rolls the i index of an array by one place
    jm1 --> Rolls the j index of an array by one place
    km1 --> Rolls the k index of an array by one place
    z_rhop --> Estimates the depth of a number of isopycnal levels (unused)
    eddy_potential_energy --> Calculates the eddy potential energy using density levels from z_rhop (unused and needs updating)
"""

from cubeprep import CubeListExtract as CLE
import iris
import iris.coord_categorisation as cat
from iris.cube import Cube
import sys
import numpy as np
from iris.coords import DimCoord, AuxCoord
from scipy.interpolate import interp1d
from fortran_lib import fortran_lib


def eddy_kinetic_energy(data_list, mask_list, tvar_window, var_dict):
    """
    eddy_kinetic_energy(data_list, mask_list, tvar_window, var_dict)

    Calculates the eddy kinetic energy per unit mass over a specified time window.

    INPUT variables
    data_list   - CubeList object of data from NEMO run
    mask_list   - CubeList object of data from mesh_mask
    tvar_window - Describes the time-window for time-averaging
                  == None        --> Use whole time range 
                  == 'month'     --> Aggregate by month
                  == 'season'    --> Aggregate by season
                  == 'year'      --> Aggregate by year
                  == n (integer) --> Rolling window of n time steps

    var_dict    - Dictionary of variable names for needed data


    OUTPUT variables
    EKE_cube - IRIS cube of the eddy kinetic energy per unit mass (t,z,y,x) [m2/s2]
             --> The cube will always have a time axis even if tvar_window = None
    EKE_zint_cube - IRIS cube of the depth-integrated eddy kinetic energy per unit mass (t,y,x) [m3/s2]
             --> The cube will always have a time axis even if tvar_window = None  
    """
    
    from cubeprep import CubeListExtract as CLE
    import iris
    import iris.coord_categorisation as cat
    from iris.cube import Cube
    import sys
    import numpy as np

    u_cube = CLE(data_list, var_dict['u'])
    v_cube = CLE(data_list, var_dict['v'])
    w_cube = CLE(data_list, var_dict['w'])
    rhop_cube = CLE(data_list, var_dict['rho'])

    umask = iris.util.squeeze(CLE(mask_list, var_dict['umask'])).data
    vmask = iris.util.squeeze(CLE(mask_list, var_dict['vmask'])).data
    tmask = iris.util.squeeze(CLE(mask_list, var_dict['tmask'])).data

    e1u = iris.util.squeeze(CLE(mask_list, var_dict['e1u'])).data
    e2u = iris.util.squeeze(CLE(mask_list, var_dict['e2u'])).data
    e3u = iris.util.squeeze(CLE(mask_list, var_dict['e3u'])).data

    e1v = iris.util.squeeze(CLE(mask_list, var_dict['e1u'])).data
    e2v = iris.util.squeeze(CLE(mask_list, var_dict['e2v'])).data
    e3v = iris.util.squeeze(CLE(mask_list, var_dict['e3v'])).data

    e1t = iris.util.squeeze(CLE(mask_list, var_dict['e1t'])).data
    e2t = iris.util.squeeze(CLE(mask_list, var_dict['e2t'])).data
    e3t = iris.util.squeeze(CLE(mask_list, var_dict['e3t'])).data

    e3w = iris.util.squeeze(CLE(mask_list, var_dict['e3w'])).data

    aux_time = u_cube.aux_coords[0]
    aux_time.rename('aux_time')

    aux_time = v_cube.aux_coords[0]
    aux_time.rename('aux_time')

    aux_time = w_cube.aux_coords[0]
    aux_time.rename('aux_time')



    #Calculate the variances over a specified time window
    if tvar_window == None:
        #Calculate variance over all time
        uvar_cube = u_cube.collapsed('time', iris.analysis.VARIANCE)
        vvar_cube = v_cube.collapsed('time', iris.analysis.VARIANCE)
        wvar_cube = w_cube.collapsed('time', iris.analysis.VARIANCE)

        uvar_cube = iris.util.as_compatible_shape(uvar_cube, u_cube)
        vvar_cube = iris.util.as_compatible_shape(vvar_cube, v_cube)
        wvar_cube = iris.util.as_compatible_shape(wvar_cube, w_cube)


    elif isinstance(tvar_window,str): #Aggregate by month, season, or year
        if tvar_window.lower() == 'month':
            #
            cat.add_month(u_cube, "time", name="month")
            cat.add_month(v_cube, "time", name="month")
            cat.add_month(w_cube, "time", name="month")
            #
            uvar_cube = u_cube.aggregated_by(["month"], iris.analysis.VARIANCE)
            vvar_cube = v_cube.aggregated_by(["month"], iris.analysis.VARIANCE)
            wvar_cube = w_cube.aggregated_by(["month"], iris.analysis.VARIANCE)
            #
        elif tvar_window.lower() == 'season':
            #
            cat.add_season(u_cube, "time", name="season")
            cat.add_season(v_cube, "time", name="season")
            cat.add_season(w_cube, "time", name="season")
            #
            uvar_cube = u_cube.aggregated_by(["season"], iris.analysis.VARIANCE)
            vvar_cube = v_cube.aggregated_by(["season"], iris.analysis.VARIANCE)
            wvar_cube = w_cube.aggregated_by(["season"], iris.analysis.VARIANCE)

        elif tvar_window.lower() == 'year':
            #
            cat.add_year(u_cube, "time", name="year")
            cat.add_year(v_cube, "time", name="year")
            cat.add_year(w_cube, "time", name="year")
            #
            uvar_cube = u_cube.aggregated_by(["year"], iris.analysis.VARIANCE)
            vvar_cube = v_cube.aggregated_by(["year"], iris.analysis.VARIANCE)
            wvar_cube = w_cube.aggregated_by(["year"], iris.analysis.VARIANCE)


    elif isinstance(tvar_window, int): #Rolling window of n time steps. n=tvar_window
        #
        uvar_cube = u_cube.rolling_window('time', iris.analysis.VARIANCE, tvar_window)
        vvar_cube = v_cube.rolling_window('time', iris.analysis.VARIANCE, tvar_window)
        wvar_cube = w_cube.rolling_window('time', iris.analysis.VARIANCE, tvar_window)

    #The variances of u,v, and w must be centred on the same T point
    eke = 0.5*tmask*( ( im1(uvar_cube.data*e1u*e2u*e3u*umask)  + (uvar_cube.data*e1u*e2u*e3u*umask) )
          + ( jm1(vvar_cube.data*e1v*e2v*e3v*vmask)  + (vvar_cube.data*e1v*e2v*e3v*vmask) )
          + ( km1(wvar_cube.data*e1t*e2t*e3w*tmask)  + (wvar_cube.data*e1t*e2t*e3w*tmask) ) ) / (2*e1t*e2t*e3t)

    eke  = np.ma.masked_array(eke, mask=np.broadcast_to(~np.ma.make_mask(tmask),eke.shape))

    #Save output array as an IRIS cube
    time_coord = uvar_cube.coord("time")
    eke_cube = Cube(eke, dim_coords_and_dims=[(time_coord,0)])
    eke_cube.var_name = "eke_" + str(tvar_window)
    eke_cube.long_name = "eddy kinetic energy"
    eke_cube.units = "m2 s-2"
    eke_cube.attributes = {'tvar_window':str(tvar_window)}

    eke_zint_cube = Cube(np.sum(eke*e3t,axis=-3), dim_coords_and_dims=[(time_coord,0)])
    eke_zint_cube.var_name = "eke_zint_" + str(tvar_window)
    eke_zint_cube.long_name = "eddy kinetic energy (zint)"
    eke_zint_cube.units = "m3 s-2"
    eke_zint_cube.attributes = {'tvar_window':str(tvar_window)}

    return eke_cube, eke_zint_cube

def im1(M): 
    import numpy as np
    output = np.roll(M, 1, axis=-1)
    return output

def jm1(M): 
    import numpy as np
    output = np.roll(M, 1, axis=-2)
    output[...,0,:] = 0.
    return output

def km1(M): 
    import numpy as np
    output = np.roll(M, 1, axis=-3)
    output[...,0,:,:] = 0.
    return output

# VVVV Currently unused functions listed below VVVVV

def z_rhop( data_list, mask_list, nn_rhop, var_dict):
    """
    Use linear interpolation to estimate the depth of isopycnals.
    """
    from cubeprep import CubeListExtract as CLE
    import iris
    from iris.coords import DimCoord, AuxCoord
    from iris.cube import Cube
    import numpy as np
    from scipy.interpolate import interp1d
    from fortran_lib import fortran_lib
    import sys

    rhop_cube = CLE(data_list, var_dict['rho'] )
    deptht_cube = iris.util.squeeze(CLE(mask_list, var_dict['deptht']))
    e3t_cube = iris.util.squeeze(CLE(mask_list, var_dict['e3t']))

    #tmask - 1 = unmasked, 0 = masked (Fortran convention)
    tmask = np.ma.make_mask(iris.util.squeeze(CLE(mask_list, var_dict['tmask'])).data)
    # tmask_zint = np.sum(tmask, axis=-3)


    rhop_min = np.min(rhop_cube.data[:,tmask])
    rhop_max = np.max(rhop_cube.data[:,tmask])
    rhop_coord = np.linspace(rhop_min, rhop_max, num=nn_rhop)
    drhop = rhop_coord[1] - rhop_coord[0]

    #Also calculate the bounds for the isopycnals
    rhop_bounds = rhop_coord - drhop/2
    rhop_bounds = np.append(rhop_bounds, rhop_coord[-1] + drhop/2)


    z_rhop, z_rhop_mask = fortran_lib.z2rhop( np.broadcast_to(deptht_cube.data.data, rhop_cube.shape), rhop_cube.data.data,
                                             tmask, rhop_coord, "extrapolate")

    bounds_rhop, bounds_rhop_mask = fortran_lib.z2rhop( np.broadcast_to(deptht_cube.data.data, rhop_cube.shape), rhop_cube.data.data,
                                             tmask, rhop_bounds, "extrapolate")


    z_rhop = np.ma.masked_array(z_rhop, mask=z_rhop_mask)
    bounds_rhop = np.ma.masked_array(bounds_rhop, mask=bounds_rhop_mask)

    z_rhop_mask_zint = np.sum(~z_rhop.mask, axis=-3)
    single_point_column_mask = np.swapaxes(np.broadcast_to((z_rhop_mask_zint == 1), (z_rhop.shape[1], z_rhop.shape[0], z_rhop.shape[2], z_rhop.shape[3])),0,1)
    single_point_mask = single_point_column_mask * ~z_rhop.mask

    max_depth = np.broadcast_to(np.sum(e3t_cube.data * tmask, axis=-3), z_rhop.shape)

    mask1 = np.ma.make_mask( ( z_rhop > max_depth ) * ~z_rhop.mask )
    mask2 = np.ma.make_mask(  ( z_rhop < 0         ) * ~z_rhop.mask )

    z_rhop[mask1] = max_depth[mask1]
    z_rhop[mask2] = 0.


    max_depth_bounds = np.broadcast_to(np.sum(e3t_cube.data * tmask, axis=-3), bounds_rhop.shape)

    mask1 = np.ma.make_mask( bounds_rhop >= max_depth_bounds )
    mask2 = np.ma.make_mask( bounds_rhop <= 0 )     


    bounds_rhop[mask1] = max_depth_bounds[mask1]
    bounds_rhop[mask2] = 0.

    e3t_rhop = np.diff(bounds_rhop, axis=-3)

    e3t_rhop[single_point_column_mask] = 0.
    e3t_rhop[single_point_mask] = max_depth[single_point_mask]

    incrop_mask = e3t_rhop <= 0
    
    z_rhop = np.ma.masked_array(z_rhop, mask=z_rhop.mask + incrop_mask)


    #Save output array as an IRIS cube
    rhop_coord = DimCoord(rhop_coord, long_name = 'Density_level', units = rhop_cube.units)
    e3t_rhop_coord = AuxCoord(e3t_rhop, long_name='IsopycnalThickness', units=deptht_cube.units)

    try: 
        time_coord = rhop_cube.coord("time")
    except:
        aux_time = rhop_cube.aux_coords[0]
        aux_time.rename("aux_time")
        time_coord = rhop_cube.coord("time")


    z_rhop_cube = Cube(z_rhop, dim_coords_and_dims=[(time_coord,0),(rhop_coord,1)])
    z_rhop_cube.var_name = deptht_cube.var_name + "_rhop2"

        
    z_rhop_cube.add_aux_coord(e3t_rhop_coord, [0,1,2,3])

    try:
        lat = rhop_cube.coord('latitude')
        lon = rhop_cube.coord('longitude')

        z_rhop_cube.add_aux_coord(lat, [2,3])
        z_rhop_cube.add_aux_coord(lon, [2,3])

    except: print("unable to add lat/lon data")

    try: 
        z_rhop_cube.standard_name = deptht_cube.standard_name
    except TypeError:
        print("Z_TO_RHOP: Original cube does not have a standard name: [", deptht_cube.var_name, "]")

    try:
        z_rhop_cube.long_name = deptht_cube.long_name + " (dense_coord)"
    except TypeError:
        print("Z_TO_RHOP: Original cube does not have a long name: [", deptht_cube.var_name,"]")

    try:
        z_rhop_cube.units = deptht_cube.units
    except TypeError:
        print("Z_TO_RHOP: Original cube does not have units [", deptht_cube.var_name,"]")

    return z_rhop_cube


def eddy_potential_energy(z_rhop_cube, data_list, mask_list, tvar_window, var_dict):
    """
    Calculate the eddy potential energy
    """
    from cubeprep import CubeListExtract as CLE
    import iris
    import iris.coord_categorisation as cat
    from iris.cube import Cube
    import numpy as np
    
    #Remove the isopycnal thickness auxiliary coordinate as it is not needed for this calculation.
    z_rhop_cube.remove_coord("IsopycnalThickness")


    #Calculate the variances over a specified time window
    if tvar_window == None:
        #Calculate variance over all time

        z_rhop_var_cube = z_rhop_cube.collapsed('time', iris.analysis.VARIANCE)
        z_rhop_var_cube = iris.util.as_compatible_shape(z_rhop_var_cube, z_rhop_cube)

    elif isinstance(tvar_window,str): #Aggregate by month, season, or year
        if tvar_window.lower() == 'month':
            #
            cat.add_month(z_rhop_cube, "time", name="month")
            #
            z_rhop_var_cube = z_rhop_cube.aggregated_by(["month"], iris.analysis.VARIANCE)
            #
        elif tvar_window.lower() == 'season':
            #
            cat.add_season(z_rhop_cube, "time", name="season")
            #
            z_rhop_var_cube = z_rhop_cube.aggregated_by(["season"], iris.analysis.VARIANCE)

        elif tvar_window.lower() == 'year':
            #
            cat.add_year(z_rhop_cube, "time", name="season")
            #
            z_rhop_var_cube = z_rhop_cube.aggregated_by(["year"], iris.analysis.VARIANCE)



    elif isinstance(tvar_window, int): #Rolling window of n time steps. n=tvar_window
        #
        z_rhop_var_cube = z_rhop_cube.rolling_window('time', iris.analysis.VARIANCE, tvar_window)

    epe = 0.5*9.81*z_rhop_var_cube.data

    drhop = np.mean(np.diff(z_rhop_cube.coord("Density_level").points))
    epe_rhoint = np.sum(epe * drhop, axis=-3)

    #Save output array as an IRIS cube
    time_coord = z_rhop_var_cube.coord("time")
    rhop_coord = z_rhop_var_cube.coord("Density_level")
    epe_cube = Cube(epe, dim_coords_and_dims=[(time_coord,0),(rhop_coord,1)])
    epe_cube.var_name = "epe_" + str(tvar_window)
    epe_cube.long_name = "eddy potential energy"
    epe_cube.units = "J m kg-1"
    epe_cube.attributes = {'tvar_window':str(tvar_window)} 

    epe_rhoint_cube = Cube(epe_rhoint, dim_coords_and_dims=[(time_coord,0)])
    epe_rhoint_cube.var_name = "epe_rhoint_" + str(tvar_window)
    epe_rhoint_cube.long_name = "eddy potential energy (rhoint)"
    epe_rhoint_cube.units = "J m-2"
    epe_rhoint_cube.attributes = {'tvar_window':str(tvar_window)} 


    return epe_cube, epe_rhoint_cube
