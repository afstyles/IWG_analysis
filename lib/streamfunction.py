"""
 _____ _                             __                  _   _                           
/  ___| |                           / _|                | | (_)                          
\ `--.| |_ _ __ ___  __ _ _ __ ___ | |_ _   _ _ __   ___| |_ _  ___  _ __    _ __  _   _ 
 `--. \ __| '__/ _ \/ _` | '_ ` _ \|  _| | | | '_ \ / __| __| |/ _ \| '_ \  | '_ \| | | |
/\__/ / |_| | |  __/ (_| | | | | | | | | |_| | | | | (__| |_| | (_) | | | |_| |_) | |_| |
\____/ \__|_|  \___|\__,_|_| |_| |_|_|  \__,_|_| |_|\___|\__|_|\___/|_| |_(_) .__/ \__, |
                                                                            | |     __/ |
                                                                            |_|    |___/ 
                                                                            
Streamfunction.py

Calculations of horizontal and overturning stream functions.

Contains methods:
   DepIntSf --> Calculates the horizontal stream function of the depth-integrated flow
   ResidualOverturning --> Calculates the residual overturning stream function
   ResOv2depth --> Projects the residual overturning circulation to depth space through interpolation
   jp1 --> Rolls the j index of an array by one place
   ZonIntSf --> Calculates the eulerian overturning stream function (unused)

"""

from cubeprep import CubeListExtract as CLE
import iris
from iris.cube import Cube
from iris.coords import AuxCoord, DimCoord
import numpy as np
from fortran_lib import fortran_lib
from scipy.interpolate import griddata


def DepIntSf(data_list, mask_list, var_dict, WG_bounds=None ):
    """
    DepIntSf(data_list, mask_list, var_dict, WG_bounds=None )

    Calculates the horizontal stream function for the depth-integrated flow 
    and calculates the gyre and ACC transport.

    INPUT variables
    data_list - CubeList object of data from NEMO run
    mask_list - CubeList object of data from mesh_mask
    var_dict - Dictionary of variable names for needed data
    WG_bounds - Tuple object (xmin, xmax, ymin, ymax) [km] describing the rectangular area which the gyre transport is calculated within
                Set to None if you want an unbounded calculation of the gyre transport 

    OUTPUT variables
    sf_zint_cube - IRIS cube of the depth-integrated stream function (t,y,x) [Sv]
    WG_transport_cube - IRIS cube of the gyre transport (t) [Sv]
    ACC_transport_cube - IRIS cube of the ACC transport (t) [Sv]
    """

    u_cube = CLE(data_list, var_dict['u']) 

    umask = iris.util.squeeze(CLE(mask_list, var_dict['umask'])).data
    e1u = iris.util.squeeze(CLE(mask_list, var_dict['e1u'])).data
    e2u = iris.util.squeeze(CLE(mask_list, var_dict['e2u'])).data
    e3u = iris.util.squeeze(CLE(mask_list, var_dict['e3u'])).data


    u = np.ma.masked_array(u_cube.data, mask=np.broadcast_to(~np.ma.make_mask(umask),u_cube.shape))
    u_zint = np.sum(u*e3u, axis=-3)

    integrand = -u_zint*e2u

    sf_zint = np.cumsum(integrand, axis=-2)

    lat = u_cube.coord("latitude")
    lon = u_cube.coord("longitude")

    #Calculate the Weddell Gyre and ACC transport for each time interval
    if isinstance(WG_bounds,tuple):
        xmin, xmax, ymin, ymax = WG_bounds
        WG_mask = ((lon.points >= xmin)*(lon.points <= xmax)*(lat.points >=ymin)*(lat.points <=ymax))
        WG_transport = np.max(sf_zint * WG_mask, axis=(-1,-2))

    else:
        WG_mask = 1.
        WG_transport = np.max(sf_zint, axis=(-1,-2))
    

    WG_transport = np.max(sf_zint * WG_mask, axis=(-1,-2))
    ACC_transport = np.median(np.sum(u_zint*e2u, axis=-2),axis=-1)

    try: 
        time_coord = u_cube.coord("time")
    except:
        aux_time = u_cube.aux_coords[0]
        aux_time.rename("aux_time")
        time_coord = u_cube.coord("time")

    lat = u_cube.coord("latitude")
    lon = u_cube.coord("longitude")

    sf_zint_cube = Cube(sf_zint/(10**6), dim_coords_and_dims=[(time_coord,0)])
    sf_zint_cube.standard_name = 'ocean_barotropic_streamfunction'
    sf_zint_cube.long_name = 'Stream function of depth-integrated flow'
    sf_zint_cube.var_name = 'sf_zint'
    sf_zint_cube.units = 'Sv'
    sf_zint_cube.description = 'ocean streamfunction variables'
    sf_zint_cube.add_aux_coord(lat, [1,2])
    sf_zint_cube.add_aux_coord(lon, [1,2])

    WG_transport_cube = Cube(WG_transport/1e6, dim_coords_and_dims=[(time_coord,0)])
    WG_transport_cube.long_name = 'Gyre transport'
    WG_transport_cube.var_name = 'WG_transport'
    WG_transport_cube.units = 'Sv'
    WG_transport_cube.attributes = {'WG_bounds': str(WG_bounds)}

    ACC_transport_cube = Cube(ACC_transport/1e6, dim_coords_and_dims=[(time_coord,0)])
    ACC_transport_cube.long_name = 'ACC transport'
    ACC_transport_cube.var_name = 'ACC_transport'
    ACC_transport_cube.units = 'Sv'

    return sf_zint_cube, WG_transport_cube, ACC_transport_cube 

def ResidualOverturning(data_list, mask_list, nn_rhop, var_dict):
    """
    Calculates the residual overturning stream function

    INPUT variables
    data_list - CubeList object of data from NEMO run
    mask_list - CubeList object of data from mesh_mask
    nn_rhop - Integer describing number of density levels to use 
              when calculating the residual overturning
    var_dict - Dictionary of variable names for needed data

    OUTPUT variables
    res_ov_cube - IRIS cube of the residual overturning stream function integrand (t,rho,y,x) [Sv] 
                  --> Includes information of the density coordinate used
                  --> Zonally integrate to calculate the residual overturning
    """

    v_cube = CLE(data_list, var_dict['v'])
    rhop_cube = CLE(data_list, var_dict['rho'] )

    tmask = iris.util.squeeze(CLE(mask_list, var_dict['tmask'])).data
    vmask = iris.util.squeeze(CLE(mask_list, var_dict['vmask'])).data

    e1v = iris.util.squeeze(CLE(mask_list, var_dict['e1v'])).data
    e2v = iris.util.squeeze(CLE(mask_list, var_dict['e2v'])).data
    e3v = iris.util.squeeze(CLE(mask_list, var_dict['e3v'])).data

    e1t = iris.util.squeeze(CLE(mask_list, var_dict['e1v'])).data
    e2t = iris.util.squeeze(CLE(mask_list, var_dict['e2v'])).data
    e3t = iris.util.squeeze(CLE(mask_list, var_dict['e3v'])).data

    # #We need to calculate the density and depth centred on the V points
    rhop_v = (tmask*rhop_cube.data*e1t*e2t + jp1(tmask*rhop_cube.data*e1t*e2t))/((tmask+jp1(tmask))*e1v*e2v)  
    rhop_v = np.ma.masked_array( rhop_v, mask = np.broadcast_to(~np.ma.make_mask(vmask), rhop_v.shape))

    rhop_min = np.min(rhop_v)
    rhop_max = np.max(rhop_v)
    rhop_coord = np.linspace(rhop_min-0.2, rhop_max+0.2, num=nn_rhop)
    drhop = rhop_coord[1] - rhop_coord[0]

    res_ov, rhop_depth, output_mask = fortran_lib.residual_overturning(v_cube.data, rhop_v.data, rhop_coord, e3v, e1v, ~np.ma.make_mask(vmask))

    res_ov = np.ma.masked_array(res_ov, mask=output_mask)
    rhop_depth = np.ma.masked_array(rhop_depth, mask=output_mask)

    rhop_coord = DimCoord(rhop_coord, long_name = 'Density_level', units = rhop_cube.units)
    rhop_depth = AuxCoord(rhop_depth, long_name='IsopycnalDepth', units='m')
    lat = v_cube.coord("latitude")
    lon = v_cube.coord("longitude")

    try: 
        time_coord = v_cube.coord("time")
    except:
        aux_time = v_cube.aux_coords[0]
        aux_time.rename("aux_time")
        time_coord = v_cube.coord("time")

    res_ov_cube = Cube(res_ov/1e6, dim_coords_and_dims=[(time_coord,0),(rhop_coord,1)])
    res_ov_cube.long_name = 'Residual overturning function'
    res_ov_cube.var_name = 'res_ov'
    res_ov_cube.units = 'Sv'
    res_ov_cube.add_aux_coord(rhop_depth, [0,1,2,3])
    res_ov_cube.add_aux_coord(lon, [2,3])
    res_ov_cube.add_aux_coord(lat, [2,3])

    return res_ov_cube

def ResOv2depth(res_ov_cube, data_list, mask_list, var_dict, nn_z=201):
    """
    Function to project the residual overturning function on to the zonal mean of the isopycnal depth

    INPUT variables
    res_ov_cube - IRIS cube of the residual overturning stream function integrand (t,rho,y,x) [Sv] 
                  --> Includes information of the density coordinate used
                  --> Zonally integrate to calculate the residual overturning
    data_list - CubeList object of data from NEMO run
    mask_list - CubeList object of data from mesh_mask
    var_dict  - Dictionary of variable names for needed data
    nn_z      - Integer describing number of regular depth coordinates to interpolate to


    OUTPUT variables
    res_ov_depth_cube - IRIS cube of the residual overturning stream function interpolated on depth space (nn_z,y) [Sv] 
                  --> Includes interpolated depth coordinate and zonal mean of the y coordinate
    """
    IsopDepth = np.ma.masked_array(res_ov_cube.coord("IsopycnalDepth").points, mask=res_ov_cube.data.mask)
    IsopDepth = np.ma.masked_greater(IsopDepth, 1e10)
    IsopDepthXmean = np.mean(np.mean(IsopDepth,axis=-1), axis=0)

    tm_res_ov = np.mean(np.sum(res_ov_cube.data, axis=-1), axis=0)
    tm_res_ov = np.ma.masked_greater(tm_res_ov, 1e10)
    tm_res_ov = np.ma.masked_less(tm_res_ov, -1e10)

    lat_coord = np.mean( res_ov_cube.coord("latitude").points, axis=-1)
    lat = np.broadcast_to( lat_coord, tm_res_ov.shape)

    depth_points = IsopDepthXmean[~tm_res_ov.mask]
    lat_points = lat[~tm_res_ov.mask]
    res_ov_points = tm_res_ov[~tm_res_ov.mask]

    depth_coord = np.linspace(0, np.max(depth_points), num=nn_z)
    
    lat_grid, depth_grid = np.meshgrid(lat_coord, depth_coord)

    out_grid = griddata( (lat_points, depth_points), res_ov_points, (lat_grid,depth_grid), method='linear')
    out_grid = np.ma.masked_invalid(out_grid) #Values that are extrapolated from the hull of the data is masked


    depth_coord = DimCoord(depth_coord, long_name = 'InterpDepth', units = res_ov_cube.coord("IsopycnalDepth").units)
    lat_coord = DimCoord(lat_coord, long_name='y_xmean', units=res_ov_cube.coord("latitude").units)

    res_ov_depth_cube = Cube(out_grid, dim_coords_and_dims=[(depth_coord,0), (lat_coord,1)])
    res_ov_depth_cube.long_name = 'Residual overturning function (depth interpolation)'
    res_ov_depth_cube.var_name = 'res_ov_depth'
    res_ov_depth_cube.units = res_ov_cube.units

    return res_ov_depth_cube



def jp1(M): 
    import numpy as np
    output = np.roll(M, -1, axis=-2)
    output[...,-1,:] = 0.

    return output

#VVVVV Currently unused functions below VVVVV


def ZonIntSF(data_list, mask_list, var_dict):
    """
    ZonIntSF(data_list, mask_list, var_dict)

    Calculates the eulerian overturning stream function

    INPUT variables
    data_list - CubeList object of data from NEMO run
    mask_list - CubeList object of data from mesh_mask
    var_dict - Dictionary of variable names for needed data

    OUTPUT variables
    sf_xint_cube - IRIS cube of the overturning of the zonally integrated flow (t,z,y) [Sv]
    """


    v_cube = CLE(data_list, var_dict['v']) 

    vmask = iris.util.squeeze(CLE(mask_list, var_dict['vmask'])).data
    vmask2d = iris.util.squeeze(CLE(mask_list, var_dict['vmask2d'])).data

    e1v = iris.util.squeeze(CLE(mask_list, var_dict['e1v'])).data
    e2v = iris.util.squeeze(CLE(mask_list, var_dict['e2v'])).data
    e3v = iris.util.squeeze(CLE(mask_list, var_dict['e3v'])).data

    v = np.ma.masked_array(v_cube.data, mask=np.broadcast_to(~np.ma.make_mask(vmask),v_cube.shape))
    vflux_xint = np.sum(v*e3v*e1v, axis=-1)

    sf_xint = np.cumsum(np.flip(vflux_xint,axis=-2), axis=-2) #We flip the depth axis so we can integrate from the bottom 

    sf_xint = np.flip(sf_xint, axis=-2) #Restore the original order of the depth axis by flipping again

    try: 
        time_coord = v_cube.coord("time")
    except:
        aux_time = v_cube.aux_coords[0]
        aux_time.rename("aux_time")
        time_coord = v_cube.coord("time")

    sf_xint_cube = Cube(sf_xint/(10**6), dim_coords_and_dims=[(time_coord,0)])
    sf_xint_cube.long_name = 'Stream function of zonally-integrated flow'
    sf_xint_cube.var_name = 'sf_xint'
    sf_xint_cube.units = 'Sv'

    return sf_xint_cube


