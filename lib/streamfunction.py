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
# import iris
# from iris.cube import Cube
# from iris.coords import AuxCoord, DimCoord
import dask
import dask.array as da
import xarray as xr
import numpy as np
from fortran_lib import fortran_lib
import warnings
from scipy.interpolate import griddata, interp1d
from SOR import sor

def DepIntSf(data_list, mask_list, var_dict, WG_bounds=None ):
    """
    DepIntSf(data_list, mask_list, var_dict, WG_bounds=None )

    Calculates the horizontal stream function for the depth-integrated flow 
    and calculates the gyre and ACC transport.

    INPUT variables
    data_list - Dataset from NEMO run
    mask_list - Dataset from mesh_mask
    var_dict - Dictionary of variable names for needed data
    WG_bounds - Tuple object (xmin, xmax, ymin, ymax) [km] describing the rectangular area which the gyre transport is calculated within
                Set to None if you want an unbounded calculation of the gyre transport 

    OUTPUT variables
    sf_zint_cube - Xarray of the depth-integrated stream function (t,y,x) [Sv]
    WG_transport_cube - Xarray of the gyre transport (t) [Sv]
    ACC_transport_cube - Xarray of the ACC transport (t) [Sv]
    """

    u_cube = data_list[var_dict['u']]

    umask = mask_list[var_dict['umask']].data
    e1u = mask_list[var_dict['e1u']].data
    e2u = mask_list[var_dict['e2u']].data
    e3u = mask_list[var_dict['e3u']].data
    lat = mask_list[var_dict['y']].data
    lon = mask_list[var_dict['x']].data

    u = u_cube.data * umask
    u_zint = (u*e3u).sum(axis=-3)

    integrand = -u_zint * e2u
    sf_zint = integrand.cumsum(axis=-2)/1e6

    sf_zint_mask = da.ma.getmaskarray(da.ma.masked_less_equal((umask.sum(axis=-3)),0.1))
    sf_zint = da.ma.masked_array( sf_zint, mask=da.broadcast_to(sf_zint_mask, sf_zint.shape))

    #Calculate the Weddell Gyre and ACC transport for each time interval
    if isinstance(WG_bounds,tuple):
        xmin, xmax, ymin, ymax = WG_bounds

        WG_mask = da.ma.getmaskarray(da.ma.masked_less(lat,ymax)
                                  *da.ma.masked_greater(lat,ymin)
                                  *da.ma.masked_less(lon,xmax)
                                  *da.ma.masked_greater(lat,xmin))

        WG_transport = (sf_zint * WG_mask).max(axis=(-1,-2))


    else:
        WG_mask = 1.
        WG_transport = sf_zint.max(axis=(-1,-2))
    
    ACC_transport = -da.median(integrand.sum(axis=-2),axis=-1)/1e6

    time_coord = u_cube.coords["time_counter"]

    sf_zint_cube = xr.DataArray(sf_zint.compute(), 
                                dims=["time_counter", "y", "x"], 
                                coords={'time_counter':time_coord.values },
                                name="sf_zint", 
                                attrs={'long_name':'Stream function of depth-integrated flow',
                                       'standard_name' : 'ocean_barotropic_streamfunction',
                                       'units' : 'Sv'})    
    

    WG_transport_cube = xr.DataArray(WG_transport.compute(),
                                    dims=["time_counter"],
                                    coords={'time_counter':time_coord.values},
                                    name="WG_transport",
                                    attrs={'units': 'Sv', 'WG_bounds': str(WG_bounds)})

    ACC_transport_cube = xr.DataArray(ACC_transport.compute(),
                                    dims=["time_counter"],
                                    coords={'time_counter':time_coord.values},
                                    name="ACC_transport",
                                    attrs={'units': 'Sv'})

    return sf_zint_cube, WG_transport_cube, ACC_transport_cube 

def ResidualOverturning(data_list, mask_list, nn_rhop, sponge_sample_dict, var_dict):
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

    if sponge_sample_dict["eiv_log"] == True:
        v_cube = data_list[var_dict['v']] + data_list[var_dict['v_eiv']] 
    else:
        v_cube = data_list[var_dict['v']] 

    rhop_cube = data_list[var_dict['rho']]
    e1v = da.squeeze(mask_list[var_dict['e1v']].data)
    e2v = da.squeeze(mask_list[var_dict['e2v']].data)
    e3v = da.squeeze(mask_list[var_dict['e3v']].data)
    e1t = da.squeeze(mask_list[var_dict['e1t']].data)
    e2t = da.squeeze(mask_list[var_dict['e2t']].data)
    e3t = da.squeeze(mask_list[var_dict['e3t']].data)
    deptht = da.squeeze(mask_list[var_dict['deptht']].data)
    tmask = da.squeeze(mask_list[var_dict['tmask']].data).astype(bool)
    vmask = da.squeeze(mask_list[var_dict['vmask']].data).astype(bool)
    y_t = da.squeeze(mask_list[var_dict['y']].data)

    #Calculate density centred on v points
    rhop_v = (tmask * rhop_cube.data * e1t * e2t) + da.roll(tmask * rhop_cube.data * e1t * e2t, -1, axis=-2)
    rhop_v = rhop_v / ((tmask.astype(int) + da.roll(tmask,-1,axis=-2).astype(int))*e1v*e2v)
    roll_mask = np.ones(vmask.shape, dtype=bool)
    roll_mask[...,-1,:] = False
    roll_mask = da.from_array(roll_mask, chunks='auto')
    rhop_v = rhop_v * roll_mask * vmask
    rhop_v = da.ma.masked_array(rhop_v, mask= da.broadcast_to(~(vmask*roll_mask),rhop_v.shape) )
    rhop_v = da.ma.masked_invalid(rhop_v)

    if sponge_sample_dict["sponge_sample_log"] == True:
        z = np.flip(np.linspace(-sponge_sample_dict["depthmax"], -sponge_sample_dict["depthmin"], num=nn_rhop))
        rhop_coord = rhop_profile(z, T_top=sponge_sample_dict["T_top"], delta_z=sponge_sample_dict["delta_z"],
                                  a0=sponge_sample_dict["a0"], T0=sponge_sample_dict["T0"], S0=sponge_sample_dict["S0"],
                                  rau0=sponge_sample_dict["rau0"])

    else:
        rhop_min = da.min(rhop_v).compute()
        rhop_max = da.max(rhop_v).compute()
        rhop_coord = np.linspace(rhop_min-0.2, rhop_max+0.2, num=nn_rhop)
        

    # drhop = rhop_coord[1] - rhop_coord[0]

    nt = rhop_cube.shape[0]
    res_ov_list = []
    rhop_depth_list = []

    @da.as_gufunc(signature="(z),(z),(p),(z),(),(z)->(p),(p),(p)", axes=[(-3),(-3),(0),(-3),(),(-3),(-3),(-3),(-3)], allow_rechunk=True)
    def gufoo(v,rho,rhoc, e3v, e1v, vmask):

        # gufunc switches axis order for the operation. When entering into the f2py module, restore original axis order
        v = np.swapaxes(v, -3, -1)
        rho = np.swapaxes(rho, -3, -1)
        e3v = np.swapaxes(e3v, -3, -1)
        vmask = np.swapaxes(vmask, -3, -1)
        e1v = np.swapaxes(e1v, -1, -2)

        res_ov, rhop_depth, output_mask = fortran_lib.residual_overturning(v, rho, rhoc, e3v, e1v, ~vmask)

        #Output back to the gufunc's desired axis order (density at the end)
        res_ov = np.swapaxes(res_ov, -1,-3)
        rhop_depth = np.swapaxes(rhop_depth, -1,-3)
        output_mask = np.swapaxes(output_mask, -1,-3)

        return res_ov, rhop_depth, output_mask.astype(bool)
    


    res_ov, rhop_depth, output_mask = gufoo(v_cube.data, rhop_v, rhop_coord, e3v, e1v, vmask)
    res_ov = da.ma.masked_array( res_ov, mask=output_mask  )
    rhop_depth = da.ma.masked_array( rhop_depth, mask=output_mask)

    time_coord = v_cube.coords["time_counter"]

    res_ov_cube = xr.DataArray((res_ov/1e6).compute(),
                               dims=["time_counter", "Density_level", "y", "x"],
                               coords={"time_counter":time_coord.values, "Density_level":rhop_coord},
                               name='res_ov',
                               attrs={'units':'Sv' })

    rhop_depth_cube = xr.DataArray(rhop_depth.compute(),
                                   dims=["time_counter", "Density_level", "y", "x"],
                                   coords={"time_counter":time_coord.values, "Density_level":rhop_coord},
                                   name='rhop_depth',
                                   attrs={'units':'m'})

    return res_ov_cube, rhop_depth_cube

def ResOv2depth(res_ov_cube, rhop_depth_cube, data_list, mask_list, var_dict, nn_z=201):
    """
    Function to project the residual overturning function on to the zonal mean of the isopycnal depth

    INPUT variables
    res_ov_cube - Xarray of the residual overturning stream function integrand (t,rho,y,x) [Sv] 
                  --> Includes information of the density coordinate used
                  --> Zonally integrate to calculate the residual overturning
    data_list - Dataset from NEMO run
    mask_list - Dataset from mesh_mask
    var_dict  - Dictionary of variable names for needed data
    nn_z      - Integer describing number of regular depth coordinates to interpolate to


    OUTPUT variables
    res_ov_depth_cube - Xarray of the residual overturning stream function interpolated on depth space (nn_z,y) [Sv] 
                  --> Includes interpolated depth coordinate and zonal mean of the y coordinate
    """
    rhop_depth = da.ma.masked_greater(da.ma.masked_invalid(rhop_depth_cube.data), 4.5e3)
    rhop_depth_xmean = da.mean(da.mean(rhop_depth, axis=-1),axis=0)
    
    e3t = da.squeeze(mask_list[var_dict['e3t']].data)
    tmask = da.squeeze(mask_list[var_dict['tmask']].data).astype(bool)


    depth_array = (e3t * tmask).sum(axis=-3)
    depth_array = da.ma.masked_less_equal(depth_array, 0 )
    depth_xmax = da.max(depth_array, axis=-1)
    depth_coord = da.linspace(0, da.max(depth_xmax).compute(), num=nn_z)


    tm_res_ov = da.mean( da.sum(da.ma.masked_invalid(res_ov_cube.data),axis=-1), axis=0 )

    tm_res_ov = da.ma.masked_greater(tm_res_ov, 1e10)
    tm_res_ov = da.ma.masked_less(tm_res_ov, -1e10)

    lat_1d = da.mean( da.squeeze(mask_list[var_dict['y']].data), axis=-1 )
    # lat = da.broadcast_to(lat_1d, tm_res_ov.shape)

    # depth_points = rhop_depth_xmean[~da.ma.getmaskarray(tm_res_ov)]
    # lat_points = lat[~da.ma.getmaskarray(tm_res_ov)]
    # res_ov_points = tm_res_ov[~da.ma.getmaskarray(tm_res_ov)]

    # depth_coord = da.linspace(0, da.max(depth_points).compute(), num=nn_z)
    # lat_grid, depth_grid = da.meshgrid(lat_1d, depth_coord)

    @da.as_gufunc(signature="(p),(p),(z)->(z)", axes=[(-2),(-2),(0),(-2)], allow_rechunk=True)
    def gufoo(tm_res_ov, rhop_depth_xmean, depth_coord):
        ny = tm_res_ov.shape[0]
        nz = len(depth_coord)
        out_grid = np.full([ny,nz], np.nan)


        for j in range(ny):
            res_ov_col = tm_res_ov[j,:]
            rhop_depth_xmean_col = rhop_depth_xmean[j,:]

            tmp_mask = np.ma.make_mask([False] + [ (rhop_depth_xmean_col[k] < np.max(rhop_depth_xmean_col[:k])) for k in range(1,len(rhop_depth_xmean_col))])
            tmp_mask = tmp_mask + np.ma.getmaskarray(res_ov_col)

            res_ov_points = res_ov_col[~tmp_mask] 
            res_ov_points = np.append(res_ov_points,[0.])

           


            rhop_depth_xmean_points = rhop_depth_xmean_col[~tmp_mask] 
            rhop_depth_xmean_points = np.append(rhop_depth_xmean_points, [0.])



            if not np.ma.is_masked(depth_xmax[j].compute()): 
                rhop_depth_xmean_points = np.append(rhop_depth_xmean_points, [depth_xmax[j].compute()])
                res_ov_points = np.append(res_ov_points,[0.])

            npoints = np.min([len(rhop_depth_xmean_points), len(res_ov_points)])

            if npoints <= 1: continue

            interp_func = interp1d(rhop_depth_xmean_points, res_ov_points, kind='linear', fill_value=np.nan, bounds_error = False)
            out_col = interp_func(depth_coord)
            out_grid[j,:] = out_col

        out_grid = np.ma.masked_invalid(out_grid)

        return out_grid

    out_grid = gufoo(tm_res_ov, rhop_depth_xmean, depth_coord)

    # print(out_grid.shape)

    # out_grid = griddata( (lat_points.compute(), depth_points.compute()), res_ov_points.compute(),
    #                      (lat_grid.compute(),depth_grid.compute()), method='linear')
    # out_grid = np.ma.masked_invalid(out_grid) #Values that are extrapolated from the hull of the data is masked

    res_ov_depth_cube = xr.DataArray(out_grid,
                                     dims=["InterpDepth", "y_xmean"],
                                     coords={"InterpDepth":depth_coord, "y_xmean":lat_1d},
                                     name='res_ov_depth',
                                     attrs={'units':res_ov_cube.attrs['units']})

    return res_ov_depth_cube


def acc_decomp(data_list, mask_list, var_dict, g=9.8, rhop_0=1027., rhop_ref=1027., lat0=-65, R=6371e3, T=86400):
    """
    Function to decompose the acc into barotropic and baroclinic parts via thermal wind.

    INPUT variables
    data_list - Dataset from NEMO run
    mask_list - Dataset from mesh_mask
    var_dict  - Dictionary of variable names for needed data


    OUTPUT variables
    Dataset including the following variables:

    Variable name | Description
    ______________|___________________________________________________________________________________
    utop_zint     : x-velocity at free surface times free surface height
    ubot_zint     : x-velocity at the sea floor multiplied by sea floor depth
    uthw_zint     : Depth-integrated x-velocity calculated from thermal wind (relative to sea floor)
    vtop_zint     : y-velocity at free surface times free surface height
    vbot_zint     : y-velocity at the sea floor multiplied by sea floor depth
    vthw_zint     : Depth-integrated y-velocity calculated from thermal wind (relative to sea floor)
    acc_top       : ACC transport calculated from utop_zint
    acc_bot       : ACC transport calculated from ubot_zint
    acc_thw       : ACC transport calculated from uthw_zint
    rhoS          :  Density on the southern boundary of the model
    rhoN          :     Density on the northern boundary of the model
    ff_tN_cube    : Coriolis parameter on the northern boundary of the model
    ff_tS_cube    : Coriolis paramater on the southern boundary of the model
    intbound_cube : Density on the northern boundary minus density on the southern boundary


    
    """
    u_cube = data_list[var_dict['u']][...,1:-1]
    v_cube = data_list[var_dict['v']][...,1:-1]
    rhop_cube = data_list[var_dict['rho']][...,1:-1]
    ssh_cube = data_list[var_dict['ssh']][...,1:-1]

    tmask = da.squeeze(mask_list[var_dict['tmask']].data).astype(bool)[...,1:-1]
    tmask2d = da.squeeze(mask_list[var_dict['tmask2d']].data).astype(bool)[...,1:-1]
    umask = da.squeeze(mask_list[var_dict['umask']].data).astype(bool)[...,1:-1]
    umask2d = da.squeeze(mask_list[var_dict['umask2d']].data).astype(bool)[...,1:-1]
    vmask = da.squeeze(mask_list[var_dict['vmask']].data).astype(bool)[...,1:-1]
    vmask2d = da.squeeze(mask_list[var_dict['vmask2d']].data).astype(bool)[...,1:-1]

    e1u = da.squeeze(mask_list[var_dict['e1u']].data)[...,1:-1]
    e2u = da.squeeze(mask_list[var_dict['e2u']].data)[...,1:-1]
    e3u = da.squeeze(mask_list[var_dict['e3u']].data)[...,1:-1]
    e1v = da.squeeze(mask_list[var_dict['e1v']].data)[...,1:-1]
    e2v = da.squeeze(mask_list[var_dict['e2v']].data)[...,1:-1]
    e3v = da.squeeze(mask_list[var_dict['e3v']].data)[...,1:-1]
    e1t = da.squeeze(mask_list[var_dict['e1t']].data)[...,1:-1]
    e2t = da.squeeze(mask_list[var_dict['e2t']].data)[...,1:-1]
    e3t = da.squeeze(mask_list[var_dict['e3t']].data)[...,1:-1]
    e1f = da.squeeze(mask_list[var_dict['e1f']].data)[...,1:-1]
    e2f = da.squeeze(mask_list[var_dict['e2f']].data)[...,1:-1]
    e3w = da.squeeze(mask_list[var_dict['e3w']].data)[...,1:-1]
    gdept_array = da.squeeze(mask_list[var_dict['deptht']].data)[...,1:-1]
    ff_t = da.squeeze(mask_list[var_dict['ff_t']].data)[...,1:-1]

    depth_array = da.sum(e3t * tmask, axis=-3)
    beta = 2 * (2*np.pi/T) * np.cos( np.pi * lat0 / 180 ) / R

    #Calculate ACC contributions from top and bottom flow

    utop, ubot = first_last_index(u_cube.data, umask, axis=-3)
    vtop, vbot = first_last_index(v_cube.data, vmask, axis=-3)

    utop_zint = utop * ssh_cube.data 
    vtop_zint = vtop * ssh_cube.data
    ubot_zint = ubot * depth_array
    vbot_zint = vbot * depth_array

    integrand_top= -utop_zint * e2u
    integrand_bot = -ubot_zint * e2u

    acc_top = da.sum(-integrand_top, axis=-2)/1e6
    acc_bot = da.sum(-integrand_bot, axis=-2)/1e6

    #Calculate the transport due to thermal wind 
    #Load densities centred on T points
    rhop_array = tmask * rhop_cube.data


    #Calculate vertical density gradient in order to calculate the correction
    drhop_dk = rhop_cube.data - da.roll(rhop_array, 1, axis=-3) # = rho[i,j,k] - rho[i,j,k-1]
    drhop_dk_mask = da.roll(tmask,1,axis=-3)*tmask # No extrapolation should be needed
    roll_mask = np.ones(drhop_dk_mask.shape, dtype=bool)
    roll_mask[0,:,:] = False
    roll_mask = da.from_array(roll_mask)

    drhop_dk_mask = drhop_dk_mask * roll_mask

    drhop_dk = drhop_dk * drhop_dk_mask

    #For densities near partial cells, interpolate the densities such that gradients are purely horizontal

    # Case 1 - jth cell is lower than the partial cell neighbour at j-1
    alpha = e3w / da.roll(e3w, 1, axis=-2)
    roll_mask = np.ones(alpha.shape, dtype=bool)
    roll_mask[:,0,:] = False
    roll_mask = da.from_array(roll_mask)
    alpha = da.ma.masked_where(  ~(tmask * da.roll(tmask,1,axis=-2)*roll_mask), alpha )
    alpha = da.ma.masked_less_equal( alpha, 1 )
    alpha = da.ma.masked_invalid(alpha)

    correction = (1 - 1/alpha)*drhop_dk
    correction = da.ma.filled(correction, fill_value=0)
    rhop_corrected_dy = rhop_array - correction

    # Case 2 - jth cell is lower than the partial cell neighbour at j+1
    alpha =  e3w / da.roll(e3w, -1, axis=-2)
    roll_mask = np.ones(alpha.shape, dtype=bool)
    roll_mask[:,-1,:] = False
    roll_mask = da.from_array(roll_mask)

    alpha = da.ma.masked_where(  ~(tmask * da.roll(tmask,-1,axis=-2)*roll_mask), alpha )
    alpha = da.ma.masked_less_equal( alpha, 1 )
    alpha = da.ma.masked_invalid(alpha)

    correction = (1 - 1/alpha)*drhop_dk
    correction = da.ma.filled(correction, fill_value=0)
    rhop_corrected_dy = rhop_corrected_dy - correction

    # Case 3 - ith cell is lower than the partial cell neighbour at i-1
    alpha = e3w / da.roll(e3w, 1, axis=-1)
    alpha = da.ma.masked_where(  ~(tmask * da.roll(tmask,1,axis=-1)), alpha )
    alpha = da.ma.masked_less_equal( alpha, 1 )
    alpha = da.ma.masked_invalid(alpha)

    correction = (1 - 1/alpha)*drhop_dk
    correction = da.ma.filled(correction, fill_value=0)
    rhop_corrected_dx = rhop_array - correction

    # Case 4 - ith cell is lower than the partial cell neighbour at i+1
    alpha =  e3w / da.roll(e3w, -1, axis=-1)
    alpha = da.ma.masked_where(  ~(tmask * da.roll(tmask,-1,axis=-2)*roll_mask), alpha )
    alpha = da.ma.masked_less_equal( alpha, 1 )
    alpha = da.ma.masked_invalid(alpha)

    correction = (1 - 1/alpha)*drhop_dk
    correction = da.ma.filled(correction, fill_value=0)
    rhop_corrected_dx = rhop_corrected_dx - correction

    #Calculate the horizontal gradients of density using the corrected density fields >>>>
    drhop_dx = umask*(da.roll(rhop_corrected_dx, -1, axis=-1) - rhop_corrected_dx ) / e1u
    drhop_dy = vmask*(da.roll(rhop_corrected_dy, -1, axis=-2) - rhop_corrected_dy ) / e2v

    drhop_dx = da.ma.masked_invalid(drhop_dx)
    drhop_dy = da.ma.masked_invalid(drhop_dy)

    #Extrapolate density gradients that approach the sea floor (free slip condition)
    east_bound_mask  = tmask * (~umask)
    west_bound_mask  = (~tmask) * (~umask)

    drhop_dx = drhop_dx + east_bound_mask * da.roll(drhop_dx, 1, axis=-1)
    drhop_dx = drhop_dx + west_bound_mask * da.roll(drhop_dx,-1, axis=-1)

    north_bound_mask = tmask * (~vmask)
    south_bound_mask = (~tmask) * (~vmask)

    drhop_dy = drhop_dy + north_bound_mask * da.roll(drhop_dy, 1, axis=-2)
    drhop_dy = drhop_dy + south_bound_mask * da.roll(drhop_dy,-1,axis=-2)

    roll_mask = np.ones(tmask2d.shape, dtype=bool)
    roll_mask[-1,:] = False
    roll_mask[0 ,:] = False
    roll_mask = da.from_array(roll_mask)

    drhop_dy = drhop_dy * roll_mask

    # Center drhop_dx on V points and center drhop_dy on U points
    roll_mask = np.ones(tmask2d.shape, dtype=bool)
    roll_mask[-1,:] = False
    roll_mask = da.from_array(roll_mask)

    drhop_dx = vmask * roll_mask * (          drhop_dx*e1u*e2u 
                                    + da.roll(drhop_dx*e1u*e2u,1 ,axis=-1)
                                    + da.roll(drhop_dx*e1u*e2u,-1,axis=-2)
                                    + da.roll(da.roll(drhop_dx*e1u*e2u,-1,axis=-2),1,axis=-1) )/(4*e1v*e2v)

    roll_mask = np.ones(tmask2d.shape, dtype=bool)
    roll_mask[0,:] = False
    roll_mask = da.from_array(roll_mask)

    drhop_dy = umask * roll_mask * (  drhop_dy*e1v*e2v
                            + da.roll(drhop_dy*e1v*e2v,-1,axis=-1)
                            + da.roll(drhop_dy*e1v*e2v,1 ,axis=-2)
                            + da.roll(da.roll(drhop_dy*e1v*e2v,-1,axis=-1),1,axis=-2) )/(4*e1u*e2u)

    #Calculate the depth of u and v points
    gdepu_array = da.cumsum( e3u * umask, axis=-3 ) - umask*e3u/2
    gdepv_array = da.cumsum( e3v * vmask, axis=-3) - vmask*e3v/2

    #Calculate 1/f centred on V and U points
    inv_ff_v = (2*e1v*e2v)/(da.roll(ff_t*e1t*e2t, -1, axis=-2) + ff_t*e1t*e2t )
    ff_v = vmask*(da.roll(ff_t*e1t*e2t, -1, axis=-2) + ff_t*e1t*e2t )/(2*e1v*e2v)
    
    inv_ff_u = (2*e1v*e2v)/(da.roll(ff_t*e1t*e2t, -1, axis=-2) + ff_t*e1t*e2t )
    ff_u = umask * (da.roll(ff_t*e1t*e2t, -1, axis=-2) + ff_t*e1t*e2t )/(2*e1v*e2v)


    #Calculate the thermal wind (depth-integrated) velocity field
    uthw_zint =   ( g*inv_ff_u/rhop_0 ) * da.sum( da.ma.masked_invalid(umask * drhop_dy * gdepu_array *  e3u) , axis=-3)
    vthw_zint =  -( g*inv_ff_v/rhop_0 ) * da.sum( da.ma.masked_invalid(vmask * drhop_dx * gdepv_array *  e3v) , axis=-3)

    uthw_zint = da.ma.masked_array(uthw_zint, mask=da.broadcast_to(~umask2d, uthw_zint.shape))
    vthw_zint = da.ma.masked_array(vthw_zint, mask=da.broadcast_to(~vmask2d, vthw_zint.shape))

    #Combine to calculate the integrand
    integrand = uthw_zint * e2u
    integrand = da.ma.masked_invalid(integrand)

    #Finally, integrate over y for the ACC transport
    acc_thw = da.sum(integrand, axis=-2)/1e6

    #We calculate the density at the Northern and Southern boundaries of the domain
    rhoS, rhoN = first_last_index( rhop_array, tmask, axis=-2 )

    #Calculate the Coriolis parameter at the Northern and Southern boundaries of the domain
    tmp_array = da.broadcast_to(ff_t, rhop_cube.shape[-3:])
    ff_tN, ff_tS = first_last_index( tmp_array, tmask, axis=-2 )

    #Calcaulte the density transformation due to intervening topography
    intbound_array = ( da.roll(rhop_array, -1, axis=-2) - rhop_array )
    intbound_array = intbound_array * ( da.roll(tmask, -1, axis=-2) * tmask )
    roll_mask = np.ones(intbound_array.shape, dtype=bool)
    roll_mask[...,-1,:] = False
    roll_mask = da.from_array(roll_mask)
    intbound_array = intbound_array * roll_mask
    intbound_array = da.sum(intbound_array, axis=-2)
    intbound_array = intbound_array - (rhoN-rhoS) # --> Subtract North-South density contrast and interior boundaries remain


    time_coord = u_cube.coords["time_counter"]

    utop_zint_cube = xr.DataArray( utop_zint.compute(),
                                dims=["time_counter", "y", "x"],
                                coords={"time_counter":time_coord.values},
                                name='utop_zint',
                                attrs={'units':'m2s-2'})

    ubot_zint_cube = xr.DataArray( ubot_zint.compute(),
                                dims=["time_counter", "y", "x"],
                                coords={"time_counter":time_coord.values},
                                name='ubot_zint',
                                attrs={'units':'m2 s-2'})

    uthw_zint_cube = xr.DataArray( uthw_zint.compute(),
                                dims=["time_counter", "y", "x"],
                                coords={"time_counter":time_coord.values},
                                name='uthw_zint',
                                attrs={'units':'m2 s-2'})

    vtop_zint_cube = xr.DataArray( vtop_zint.compute(),
                                dims=["time_counter", "y", "x"],
                                coords={"time_counter":time_coord.values},
                                name='vtop_zint',
                                attrs={'units':'m2 s-2'})

    vbot_zint_cube = xr.DataArray( vbot_zint.compute(),
                                dims=["time_counter", "y", "x"],
                                coords={"time_counter":time_coord.values},
                                name='vbot_zint',
                                attrs={'units':'m2 s-2'})
    
    vthw_zint_cube = xr.DataArray( vthw_zint.compute(),
                                dims=["time_counter", "y", "x"],
                                coords={"time_counter":time_coord.values},
                                name='vthw_zint',
                                attrs={'units':'m2 s-2'})

    
    acc_top_cube = xr.DataArray(acc_top.compute(),
                                dims=["time_counter", "x"],
                                coords={"time_counter":time_coord.values},
                                name='acc_top',
                                attrs={'units':'Sv'})

    acc_bot_cube = xr.DataArray(acc_bot.compute(),
                                dims=["time_counter", "x"],
                                coords={"time_counter":time_coord.values},
                                name='acc_bot',
                                attrs={'units':'Sv'})

    acc_thw_cube = xr.DataArray( acc_thw.compute(),
                                dims=["time_counter", "x"],
                                coords={"time_counter":time_coord.values},
                                name='acc_thw',
                                attrs={'units':'Sv'})    

    rhoS_cube = xr.DataArray( rhoS.compute(),
                                dims=["time_counter", "model_level","x"],
                                coords={"time_counter":time_coord.values},
                                name='rhos',
                                attrs={'units':rhop_cube.attrs['units']})
    
    
    rhoN_cube = xr.DataArray( rhoN.compute(),
                                dims=["time_counter", "model_level","x"],
                                coords={"time_counter":time_coord.values},
                                name='rhon',
                                attrs={'units':rhop_cube.attrs['units']})

    ff_tN_cube = xr.DataArray( ff_tN.compute(),
                                dims=["model_level","x"],
                                coords={},
                                name='ff_tn',
                                attrs={'units':"s-1"})

    ff_tS_cube = xr.DataArray( ff_tS.compute(),
                                dims=["model_level","x"],
                                coords={},
                                name='ff_ts',
                                attrs={'units':"s-1"})

    intbound_cube = xr.DataArray( intbound_array.compute(),
                                dims=["time_counter", "model_level", "x"],
                                coords={"time_counter":time_coord.values},
                                name='intbnd_rhop',
                                attrs={'units':rhop_cube.attrs['units']})    

    return [utop_zint_cube, ubot_zint_cube, uthw_zint_cube, 
            vtop_zint_cube, vbot_zint_cube, vthw_zint_cube,
            acc_top_cube  , acc_bot_cube  , acc_thw_cube  ,
            rhoS_cube     , rhoN_cube     ,
            ff_tN_cube    , ff_tS_cube    , intbound_cube ]      


def first_last_index(array, mask, axis):
    array = da.swapaxes(array, axis, -1)
    mask = da.swapaxes(mask, axis, -1)

    inds = da.squeeze(da.indices(mask.shape))

    ind_target = inds[-1]
    ind_target = da.ma.masked_array(ind_target, mask=~mask)

    ind_1_array = da.broadcast_to(ind_target.min(axis=-1, keepdims=True), mask.shape)
    ind_2_array = da.broadcast_to(ind_target.max(axis=-1, keepdims=True), mask.shape)

    output_1 = da.ma.masked_where(da.broadcast_to(ind_target != ind_1_array, array.shape), array)
    output_1 = output_1.sum(axis=-1, keepdims=True)

    output_2 = da.ma.masked_where(da.broadcast_to(ind_target != ind_2_array, array.shape), array)
    output_2 = output_2.sum(axis=-1, keepdims=True)

    output_1 = da.squeeze(da.swapaxes(output_1, axis, -1),axis=axis)
    output_2 = da.squeeze(da.swapaxes(output_2, axis, -1),axis=axis)

    return output_1, output_2



def WG_decomp( data_list, acc_decomp_cube_list, mask_list, var_dict):
    """
    Calculate the velocity potential associated with the compressible component of the 
    full, baroclinic, and barotropic flow.
    Use this to calculate a corrected stream function

    INPUT variables
    data_list - Dataset from NEMO run
    acc_decomp_cube_list - Dataset from acc_decomp
    mask_list - Dataset from mesh_mask
    var_dict  - Dictionary of variable names for needed data


    OUTPUT variables
    Dataset including the velocity potential from the following depth-integrated flows:

    Variable name | Description
    ______________|___________________________________________________________________________________
    phi_full_cube | Velocity potential of the depth-integrated flow
    phi_thw_cube  | Velocity potential of the depth-intregrated flow calculated from thermal wind (relative to the sea floor)
    phi_top_cube  | Velocity potential of the surface velocity multiplied by the free surface height
    phi_bot_cube  | Velocity potential of the sea floor velocity multiplied by the depth of the sea floor
    
    """
    from SOR import sor

    #Load the time-averaged velocity fields
    u_cube = data_list[var_dict['u']].mean("time_counter")[...,1:-1]
    v_cube = data_list[var_dict['v']].mean("time_counter")[...,1:-1]

    uthw_cube = acc_decomp_cube_list['uthw_zint'].mean("time_counter")
    vthw_cube = acc_decomp_cube_list['vthw_zint'].mean("time_counter")

    utop_cube = acc_decomp_cube_list['utop_zint'].mean("time_counter")
    vtop_cube = acc_decomp_cube_list['vtop_zint'].mean("time_counter")

    ubot_cube = acc_decomp_cube_list['ubot_zint'].mean("time_counter")
    vbot_cube = acc_decomp_cube_list['vbot_zint'].mean("time_counter")

    #Load the grid information
    umask = da.squeeze(mask_list[var_dict['umask']].data).astype(bool)[...,1:-1]
    vmask = da.squeeze(mask_list[var_dict['vmask']].data).astype(bool)[...,1:-1]

    tmask2d = da.squeeze(mask_list[var_dict['tmask2d']].data).astype(bool)[...,1:-1]
    umask2d = da.squeeze(mask_list[var_dict['umask2d']].data).astype(bool)[...,1:-1]
    vmask2d = da.squeeze(mask_list[var_dict['vmask2d']].data).astype(bool)[...,1:-1]

    e1u = da.squeeze(mask_list[var_dict['e1u']].data)[...,1:-1]
    e2u = da.squeeze(mask_list[var_dict['e2u']].data)[...,1:-1]
    e3u = da.squeeze(mask_list[var_dict['e3u']].data)[...,1:-1]

    e1v = da.squeeze(mask_list[var_dict['e1v']].data)[...,1:-1]
    e2v = da.squeeze(mask_list[var_dict['e2v']].data)[...,1:-1]
    e3v = da.squeeze(mask_list[var_dict['e3v']].data)[...,1:-1]

    e1t = da.squeeze(mask_list[var_dict['e1t']].data)[...,1:-1]
    e2t = da.squeeze(mask_list[var_dict['e2t']].data)[...,1:-1]

    #Depth integrate the full velocity fields
    u_array = da.sum( u_cube.data * umask * e3u , axis=-3)
    v_array = da.sum( v_cube.data * vmask * e3v , axis=-3)

    #Set masked velocities to zero
    uthw_array = da.ma.filled(da.ma.masked_invalid(uthw_cube.data), fill_value=0.)
    vthw_array = da.ma.filled(da.ma.masked_invalid(vthw_cube.data), fill_value=0.)

    utop_array = da.ma.filled(da.ma.masked_invalid(utop_cube.data), fill_value=0.)
    vtop_array = da.ma.filled(da.ma.masked_invalid(vtop_cube.data), fill_value=0.)

    ubot_array = da.ma.filled(da.ma.masked_invalid(ubot_cube.data), fill_value=0.)
    vbot_array = da.ma.filled(da.ma.masked_invalid(vbot_cube.data), fill_value=0.)

    #Calculate the compressible part of the full flow >>>>>

    ny = umask2d.shape[-2]
    nx = umask2d.shape[-1]

    bd_xminvals = np.zeros(nx)
    bd_xmaxvals = np.zeros(nx)
    bd_yminvals = np.zeros(ny)
    bd_ymaxvals = np.zeros(ny)
    bd_xminnm = True   #        Note! - Currently sor_cgrid does not work for non-zero Neumann boundary conditions
    bd_xmaxnm = True   #        i.e. (Non-zero volume fluxes on the boundaries). This is fine for the CANAL model
    bd_yminnm = True   #        as it is closed and periodic but annoying for any future work.    
    bd_ymaxnm = True   #        
    bd_xpd = True
    bd_ypd = False

    omega = 1.3
    thresh = 1e-2
    niterations_max = int(1e6)


    #Divergence of the flow
    roll_mask = np.ones(umask2d.shape, dtype=bool)
    roll_mask[0,:] = False
    roll_mask = da.from_array(roll_mask)

    
    #Calculate the velocity potential using Successive Over Relaxation (SOR) >>>>

    #Divergence of the full flow
    print("Solving for full flow")
    D = tmask2d*roll_mask*( u_array * e2u - da.roll(u_array * e2u, 1, axis=-1)
                          + v_array * e1v - da.roll(v_array * e1v, 1, axis=-2) )/(2*e1t*e2t)
    D = da.ma.masked_invalid(D)
    D = da.ma.filled(D, fill_value=0.)

    phi_zero = np.zeros(D.shape)

    phi_full = sor.sor_cgrid(D.compute(), e1u.compute(), e2u.compute(), e1v.compute(),
                    e2v.compute(), e1t.compute(), e2t.compute(), phi_zero,
                    bd_xminvals, bd_xmaxvals, bd_yminvals, bd_ymaxvals,
                    bd_xminnm , bd_xmaxnm, bd_yminnm, bd_ymaxnm,
                    bd_xpd, bd_ypd, (~tmask2d).compute(), omega, thresh, niterations_max   )
    phi_full = da.from_array(phi_full)

    #Divergence of the thermal wind component
    print("Solving for thermal wind component")
    D = tmask2d*roll_mask*( uthw_array * e2u - da.roll(uthw_array * e2u, 1, axis=-1)
                          + vthw_array * e1v - da.roll(vthw_array * e1v, 1, axis=-2) )/(2*e1t*e2t)
    D = da.ma.masked_invalid(D)
    D = da.ma.filled(D, fill_value=0.)

    phi_zero = np.zeros(D.shape)
    phi_thw = sor.sor_cgrid(D.compute(), e1u.compute(), e2u.compute(), e1v.compute(),
                    e2v.compute(), e1t.compute(), e2t.compute(), phi_zero,
                    bd_xminvals, bd_xmaxvals, bd_yminvals, bd_ymaxvals,
                    bd_xminnm , bd_xmaxnm, bd_yminnm, bd_ymaxnm,
                    bd_xpd, bd_ypd, (~tmask2d).compute(), omega, thresh, niterations_max   )
    phi_thw = da.from_array(phi_thw)

    #Divergence of the top flow
    print("Solving for top flow")

    D = tmask2d*roll_mask*( utop_array * e2u - da.roll(utop_array * e2u, 1, axis=-1)
                          + vtop_array * e1v - da.roll(vtop_array * e1v, 1, axis=-2) )/(2*e1t*e2t)
    D = da.ma.masked_invalid(D)
    D = da.ma.filled(D, fill_value=0.)

    phi_zero = np.zeros(D.shape)
    phi_top = sor.sor_cgrid(D.compute(), e1u.compute(), e2u.compute(), e1v.compute(),
                    e2v.compute(), e1t.compute(), e2t.compute(), phi_zero,
                    bd_xminvals, bd_xmaxvals, bd_yminvals, bd_ymaxvals,
                    bd_xminnm , bd_xmaxnm, bd_yminnm, bd_ymaxnm,
                    bd_xpd, bd_ypd, (~tmask2d).compute(), omega, thresh, niterations_max   )
    phi_top = da.from_array(phi_top)

    #Divergence of the bottom flow
    print("Solving for bottom flow")
    D = tmask2d*roll_mask*( ubot_array * e2u - da.roll(ubot_array * e2u, 1, axis=-1)
                          + vbot_array * e1v - da.roll(vbot_array * e1v, 1, axis=-2) )/(2*e1t*e2t)
    D = da.ma.masked_invalid(D)
    D = da.ma.filled(D, fill_value=0.)

    phi_zero = np.zeros(D.shape)
    phi_bot = sor.sor_cgrid(D.compute(), e1u.compute(), e2u.compute(), e1v.compute(),
                    e2v.compute(), e1t.compute(), e2t.compute(), phi_zero,
                    bd_xminvals, bd_xmaxvals, bd_yminvals, bd_ymaxvals,
                    bd_xminnm , bd_xmaxnm, bd_yminnm, bd_ymaxnm,
                    bd_xpd, bd_ypd, (~tmask2d).compute(), omega, thresh, niterations_max   )
    phi_bot = da.from_array(phi_bot)

    phi_full_cube = xr.DataArray( da.ma.masked_array(phi_full, mask=~tmask2d ).compute(),
                                 dims=["y", "x"],
                                 coords={},
                                 name="phi_full",
                                 attrs={'units':'m3 s-1'}  )

    phi_thw_cube = xr.DataArray( da.ma.masked_array(phi_thw, mask=~tmask2d ).compute(),
                                 dims=["y", "x"],
                                 coords={},
                                 name="phi_thw",
                                 attrs={'units':'m3 s-1'}  )

    phi_top_cube = xr.DataArray( da.ma.masked_array(phi_top, mask=~tmask2d ).compute(),
                                 dims=["y", "x"],
                                 coords={},
                                 name="phi_top",
                                 attrs={'units':'m3 s-1'}  )

    phi_bot_cube = xr.DataArray( da.ma.masked_array(phi_bot, mask=~tmask2d ).compute(),
                                 dims=["y", "x"],
                                 coords={},
                                 name="phi_bot",
                                 attrs={'units':'m3 s-1'}  )

    return [phi_full_cube, phi_thw_cube, phi_top_cube, phi_bot_cube]

def u_corr( data_list, acc_decomp_cube_list, phi_cube_list, sf_zint_cube, WG_bounds, mask_list, var_dict ):
    """
    Calculate the velocity of the compressible field using the potential calculated by WG_decomp.
    Use this corrected velocity to find the corrected stream function.
    """

    e1u = da.squeeze(mask_list[var_dict['e1u']].data)[...,1:-1]
    e2u = da.squeeze(mask_list[var_dict['e2u']].data)[...,1:-1]
    umask2d = da.squeeze(mask_list[var_dict['umask2d']].data)[...,1:-1]

    lat = da.squeeze(mask_list[var_dict['y']].data)[...,1:-1]
    lon = da.squeeze(mask_list[var_dict['x']].data)[...,1:-1]

    cube_list = []

    for lab in ['full', 'thw', 'top', 'bot']:

        print("lab = ", lab)

        phi      = da.ma.masked_invalid(da.squeeze(phi_cube_list['phi_'+lab].data))
        ucompr = (da.roll(phi, -1, axis=-1) - phi )/e1u

        if lab =='full':
            integrand = -ucompr * e2u
            sf_corr = da.ma.masked_invalid(sf_zint_cube.mean("time_counter").data[...,1:-1]) - integrand.cumsum(axis=-2)/1e6
        
        else:
            u_zint   = acc_decomp_cube_list['u'+lab+'_zint'].mean("time_counter")
            u_zint   = da.squeeze(u_zint.data)
            
            #Calculate rotational part of the flow
            urot = u_zint - ucompr 
            integrand = -urot * e2u
            sf_corr = integrand.cumsum(axis=-2)/1e6

        #Calculate the Weddell Gyre and ACC transport for each time interval
        if isinstance(WG_bounds,tuple):
            xmin, xmax, ymin, ymax = WG_bounds

            WG_mask = da.ma.getmaskarray(da.ma.masked_less(lat,ymax)
                                    *da.ma.masked_greater(lat,ymin)
                                    *da.ma.masked_less(lon,xmax)
                                    *da.ma.masked_greater(lat,xmin)).astype(bool)

            WG_transport = (sf_corr * WG_mask).max(axis=(-1,-2))


        else:
            WG_mask = True
            WG_transport = sf_corr.max(axis=(-1,-2))



        print(WG_transport.compute())

        sf_corr_cube = xr.DataArray(sf_corr.compute(), 
                                dims=[ "y", "x"], 
                                coords={},
                                name="sf_corr_"+lab, 
                                attrs={'long_name':'Stream function of corrected flow ('+lab+')',
                                       'standard_name' : 'ocean_barotropic_streamfunction',
                                       'units' : 'Sv'})    
    

        WG_corr_cube = xr.DataArray(WG_transport.compute(),
                                    name="WG_corr_transport_"+lab,
                                    attrs={'units': 'Sv', 'WG_bounds': str(WG_bounds)})

        cube_list.append(sf_corr_cube)
        cube_list.append(WG_corr_cube)    

    return cube_list  

def rhop_profile(z, T_top=10., delta_z=1500., a0=0.28e-3, T0=10., S0=35., rau0=1026.):

    Tprofile = T_top * np.exp(z/delta_z)

    output = rau0 * ( 1 - a0 * (Tprofile - T0))

    return output

def jp1(M): 
    output = np.roll(M,-1,axis=-2)
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


# Unused code below  VVVV
"""
    # Calculate the transport due to f-plane thermal wind
    ni = e3u.shape[-1]
    nk = e3u.shape[-3]
    canyon_labels = umask * np.cumsum(umask, axis=-2)
    canyon_labels = canyon_labels.astype(int)

    # tmp_array =  -g * gdept_array * rhop_cube.data / (rhop_0 * ff_t)
    # tmp_array = e3u*( np.roll(tmp_array * e1t * e2t, -1, axis=-1) + tmp_array * e1t * e2t ) / (2 * e1u * e2u) #Centre the array on u points

    # tmp_array = (0*rhop_cube.data - rhop_ref)
    # tmp_array = ( np.roll(tmp_array * e1t * e2t, -1, axis=-1) + tmp_array * e1t * e2t ) / (2 * e1u * e2u)
    # tmpS, tmpN = first_last_index( tmp_array , umask, axis=-2)

    # print("Density anomalies at N and S >>>")
    # print("N: ", np.mean(tmpN[:,:,10], axis=0))
    # print("S: ", np.mean(tmpS[:,:,10], axis=0))
    # print("N-S: ", np.mean(tmpN[:,:,10] - tmpS[:,:,10], axis=0))


    # tmp_array = np.broadcast_to(ff_t, rhop_cube.shape)
    # tmp_array = ( np.roll(tmp_array * e1t * e2t, -1, axis=-1) + tmp_array * e1t * e2t ) / (2 * e1u * e2u)
    # tmpS, tmpN = first_last_index( tmp_array , umask, axis=-2)

    # print("f at N and S >>>")
    # print("N: ", np.mean(tmpN[:,:,10], axis=0))
    # print("S: ", np.mean(tmpS[:,:,10], axis=0))
    # print("N-S: ", np.mean(tmpN[:,:,10] - tmpS[:,:,10], axis=0))

    # tmp_array = (0*rhop_cube.data - rhop_ref)/ff_t
    # tmp_array = ( np.roll(tmp_array * e1t * e2t, -1, axis=-1) + tmp_array * e1t * e2t ) / (2 * e1u * e2u)
    # tmpS, tmpN = first_last_index( tmp_array , umask, axis=-2)

    # print("Density / f anomalies at N and S >>>")
    # print("N: ", np.mean(tmpN[:,:,10], axis=0))
    # print("S: ", np.mean(tmpS[:,:,10], axis=0))
    # print("N-S: ", np.mean(tmpN[:,:,10] - tmpS[:,:,10], axis=0))

    # tmp_array = - np.broadcast_to(gdept_array, rhop_cube.shape)
    # tmp_array = ( np.roll(tmp_array * e1t * e2t, -1, axis=-1) + tmp_array * e1t * e2t ) / (2 * e1u * e2u)
    # tmpS, tmpN = first_last_index( tmp_array , umask, axis=-2)

    # print("Cell depths at N and S >>>")
    # print("N: ", np.mean(tmpN[:,:,10], axis=0))
    # print("S: ", np.mean(tmpS[:,:,10], axis=0))
    # print("N-S: ", np.mean(tmpN[:,:,10] - tmpS[:,:,10], axis=0))



    tmp_array =  g * gdept_array * (  rhop_cube.data - rhop_ref) / (rhop_0 * ff_t) 
    tmp_array = e3u * ( np.roll(tmp_array * e1t * e2t, -1, axis=-1) + tmp_array * e1t * e2t ) / (2 * e1u * e2u) #Centre the array on u points
    print("tmp_array.shape = ", tmp_array.shape)

    # Sbound_array, Nbound_array = first_last_index( tmp_array , umask, axis=-2)
    
    # print("gdept z variation")
    # print(gdept_array[:,35,10])
    
    # print("density variation")
    # print(rhop_cube.data[-1,:,35,10])

    # print("Nbound z variation")
    # print(Nbound_array[-1,:,10]/1e6)

    # print("Sbound z variation")
    # print(Sbound_array[-1,:,10]/1e6)
    
    # Nbound_array = np.sum(Nbound_array, axis=-2)/1e6 #Vertical sum
    # Sbound_array = -np.sum(Sbound_array, axis=-2)/1e6 #Vertical sum
    




    Calculate the effect of intermediate boundaries. This is done by summing the differenced array and
    comparing the value to Nbound - Sbound. In cases where there are no intermediate boundaries, these
    two values will be the same.

    intbound_array = ( np.roll(tmp_array, -1, axis=-2) - tmp_array )
    intbound_array = intbound_array * ( np.roll(umask, -1, axis=-2) * umask )
    intbound_array[...,-1,:] = 0.

    intbound_array = np.sum(intbound_array, axis=(-3,-2)) / 1e6

    #Identify cells that have been overshadowed by neighbouring partial cells >>>>>>>>>>>>
    
    #Northern boundaries (smaller partial cell at j = j+1)
    alpha = umask * np.roll(e3u, -1, axis=-2) / e3u
    alpha[:,-1,:] = 0.

    missed_cell_mask = (alpha < 1)*(alpha > 0)
    missed_contribution = np.sum((1-alpha)*tmp_array*missed_cell_mask, axis=(-3,-2)) / 1e6
    print("Missed contribution (1) = ", np.mean(missed_contribution[:,12]))

    #Southern boundaries (smaller partial cell at j = j+1)
    alpha = umask * np.roll( umask*e3u, 1, axis=-2) / e3u
    alpha[:,0,:] = 0.

    missed_cell_mask = (alpha < 1)*(alpha > 0)
    missed_contribution = missed_contribution - np.sum((1-alpha)*tmp_array*missed_cell_mask, axis=(-3,-2)) / 1e6


    print("Missed contribution (2) = ", np.mean(missed_contribution[:,12]))

    intbound_array = intbound_array + missed_contribution

    #Redefine the tmp array so that it contains the compensating term for beta effects
    tmp_array =  g * gdept_array * ( rhop_cube.data - rhop_ref) * beta / (rhop_0 * (ff_t **2) ) 

    #Note this is the geometric f cell thickness NOT the f cell thickness used to calculate the Coriolis acceleration in dynvor.F90 (NEMO src)
    e3f = ( np.roll(e3u*e1u*e2u, -1, axis=-2) + e3u*e1u*e2u )/(2*e1f*e2f) #e3f outputs in NEMO can be dodgy, generate our own
    fmask = np.roll(umask, -1, axis=-2) * umask  #Conservative f mask. If either u cell is masked, ignore.
    tmp_array = e3f*e2f*( np.roll(tmp_array * e1u*e2u, -1, axis=-2) + tmp_array*e1u*e2u )/(2*e1f*e2f)

    beta_array = np.sum( tmp_array * fmask, axis=(-3,-2) )/1e6


    # tmp_array = e3u * e2u * ( np.roll(tmp_array * e1t * e2t, -1, axis=-1) + tmp_array * e1t * e2t ) / (2 * e1u * e2u) #Centre the array on u points 
    # beta_array = np.sum( tmp_array * umask, axis=(-3,-2))/1e6

    print("SNI diagnostics >>")
    print("ACC_SSH: ", np.mean(acc_ssh[:,12]), " Sv")
    print("ACC_BOT: ", np.mean(acc_bot[:,12]), " Sv")
    # print("Sbound_array: ", np.mean(Sbound_array[:,12]), " Sv")
    # print("Nbound_array: ", np.mean(Nbound_array[:,12]), " Sv")
    print("Intbound_array: ", np.mean(intbound_array[:,12]), " Sv")
    print("Beta_array: ", np.mean(beta_array[:,12]), " Sv")
    
    try: 
        time_coord = u_cube.coord("time")
    except:
        aux_time = u_cube.aux_coords[0]
        aux_time.rename("aux_time")
        time_coord = u_cube.coord("time")

    time_coord = u_cube.coord("time")

    acc_ssh_cube = Cube(acc_ssh, dim_coords_and_dims=[(time_coord,0)])
    acc_ssh_cube.var_name = 'acc_ssh'
    acc_ssh_cube.units = 'Sv'

    acc_bot_cube = Cube(acc_bot, dim_coords_and_dims=[(time_coord,0)])
    acc_bot_cube.var_name = 'acc_bot'
    acc_bot_cube.units = 'Sv'

    # Sbound_cube = Cube(Sbound_array, dim_coords_and_dims=[(time_coord,0)])
    # Sbound_cube.var_name = 'sbound'
    # Sbound_cube.units = 'Sv'

    # Nbound_cube = Cube(Nbound_array, dim_coords_and_dims=[(time_coord,0)])
    # Nbound_cube.var_name = 'nbound'
    # Nbound_cube.units = 'Sv'

    intbound_cube = Cube(intbound_array, dim_coords_and_dims=[(time_coord,0)])
    intbound_cube.var_name = 'intbound'
    intbound_cube.units = 'Sv'

    beta_cube = Cube(beta_array, dim_coords_and_dims=[(time_coord,0)])
    beta_cube.var_name = 'beta_comp'
    beta_cube.units = 'Sv'

    # return [acc_ssh_cube, acc_bot_cube, Sbound_cube, Nbound_cube, intbound_cube, beta_cube]

    return [acc_ssh_cube, acc_bot_cube, intbound_cube, beta_cube]
"""

