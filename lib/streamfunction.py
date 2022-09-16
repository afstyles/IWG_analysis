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
from xgcm import Grid
import numpy as np
from fortran_lib import fortran_lib
import warnings
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

    #Calculate density centred on v points
    rhop_v = (tmask * rhop_cube.data * e1t * e2t) + da.roll(tmask * rhop_cube.data * e1t * e2t, -1, axis=-2)
    rhop_v = rhop_v / ((tmask.astype(int) + da.roll(tmask,-1,axis=-2).astype(int))*e1v*e2v)
    roll_mask = np.ones(vmask.shape, dtype=bool)
    roll_mask[...,-1,:] = False
    roll_mask = da.from_array(roll_mask, chunks='auto')
    rhop_v = rhop_v * roll_mask * vmask
    rhop_v = da.ma.masked_array(rhop_v, mask= da.broadcast_to(~(vmask*roll_mask),rhop_v.shape) )
    rhop_v = da.ma.masked_invalid(rhop_v)

    rhop_min = da.min(rhop_v).compute()
    rhop_max = da.max(rhop_v).compute()
    rhop_coord = np.linspace(rhop_min-0.2, rhop_max+0.2, num=nn_rhop)

    drhop = rhop_coord[1] - rhop_coord[0]

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
    rhop_depth = da.ma.masked_greater(da.ma.masked_invalid(rhop_depth_cube.data), 1e10)
    rhop_depth_xmean = da.mean(da.mean(rhop_depth, axis=-1),axis=0)

    tm_res_ov = da.mean( da.sum(da.ma.masked_invalid(res_ov_cube.data),axis=-1), axis=0 )

    
    tm_res_ov = da.ma.masked_greater(tm_res_ov, 1e10)
    tm_res_ov = da.ma.masked_less(tm_res_ov, -1e10)

    lat_1d = da.mean( da.squeeze(mask_list[var_dict['y']].data), axis=-1 )
    lat = da.broadcast_to(lat_1d, tm_res_ov.shape)


    depth_points = rhop_depth_xmean[~da.ma.getmaskarray(tm_res_ov)]
    lat_points = lat[~da.ma.getmaskarray(tm_res_ov)]
    res_ov_points = tm_res_ov[~da.ma.getmaskarray(tm_res_ov)]

    depth_coord = da.linspace(0, da.max(depth_points).compute(), num=nn_z)
    lat_grid, depth_grid = da.meshgrid(lat_1d, depth_coord)

    out_grid = griddata( (lat_points.compute(), depth_points.compute()), res_ov_points.compute(),
                         (lat_grid.compute(),depth_grid.compute()), method='linear')
    out_grid = np.ma.masked_invalid(out_grid) #Values that are extrapolated from the hull of the data is masked

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
    data_list - CubeList object of data from NEMO run
    mask_list - CubeList object of data from mesh_mask
    var_dict  - Dictionary of variable names for needed data


    OUTPUT variables

    """
    u_cube = data_list[var_dict['u']]
    rhop_cube = data_list[var_dict['rho']]
    ssh_cube = data_list[var_dict['ssh']]

    tmask = da.squeeze(mask_list[var_dict['tmask']].data).astype(bool)
    tmask2d = da.squeeze(mask_list[var_dict['tmask2d']].data).astype(bool)
    umask = da.squeeze(mask_list[var_dict['umask']].data).astype(bool)
    umask2d = da.squeeze(mask_list[var_dict['umask2d']].data).astype(bool)
    vmask = da.squeeze(mask_list[var_dict['vmask']].data).astype(bool)
    vmask2d = da.squeeze(mask_list[var_dict['vmask2d']].data).astype(bool)

    e1u = da.squeeze(mask_list[var_dict['e1u']].data)
    e2u = da.squeeze(mask_list[var_dict['e2u']].data)
    e3u = da.squeeze(mask_list[var_dict['e3u']].data)
    e1v = da.squeeze(mask_list[var_dict['e1v']].data)
    e2v = da.squeeze(mask_list[var_dict['e2v']].data)
    e3v = da.squeeze(mask_list[var_dict['e3v']].data)
    e1t = da.squeeze(mask_list[var_dict['e1t']].data)
    e2t = da.squeeze(mask_list[var_dict['e2t']].data)
    e3t = da.squeeze(mask_list[var_dict['e3t']].data)
    e1f = da.squeeze(mask_list[var_dict['e1f']].data)
    e2f = da.squeeze(mask_list[var_dict['e2f']].data)
    e3w = da.squeeze(mask_list[var_dict['e3w']].data)
    gdept_array = da.squeeze(mask_list[var_dict['deptht']].data)
    ff_t = da.squeeze(mask_list[var_dict['ff_t']].data)

    depth_array = da.sum(e3t * tmask, axis=-3)
    beta = 2 * (2*np.pi/T) * np.cos( np.pi * lat0 / 180 ) / R

    #Calculate ACC contributions from top and bottom flow

    u_top, u_bot = first_last_index(u_cube.data, umask, axis=-3)
    
    acc_ssh = da.sum(u_top * ssh_cube.data * e2u, axis=-2)/1e6
    acc_bot = da.sum(u_bot * depth_array * e2u, axis=-2)/1e6

    #Calculate the transport due to thermal wind 
    #Densities centred on v points
    rhop_array = tmask * rhop_cube.data

    drhop_dk = rhop_cube.data - da.roll(rhop_array, 1, axis=-3) # = rho[i,j,k] - rho[i,j,k-1]
    drhop_dk_mask = da.roll(tmask,1,axis=-3)*tmask # No extrapolation should be needed
    roll_mask = np.ones(drhop_dk_mask.shape, dtype=bool)
    roll_mask[0,:,:] = False
    roll_mask = da.from_array(roll_mask)

    drhop_dk_mask = drhop_dk_mask * roll_mask

    drhop_dk = drhop_dk * drhop_dk_mask

    #For densities near partial cells, interpolate the densities such that gradients are purely horizontal

    #Case 1 - jth cell is lower than the partial cell neighbour at j-1
    alpha = da.roll(tmask,1,axis=-2) * tmask * e3w / da.roll(e3w, 1, axis=-2)
    roll_mask = np.ones(alpha.shape, dtype=bool)
    roll_mask[:,0,:] = False
    roll_mask = da.from_array(roll_mask)
    alpha = alpha * roll_mask

    tmp_mask = alpha > 1
    correction = tmp_mask * (1 - 1/alpha)*drhop_dk


    #Case 2 - jth cell is lower than the partial cell neighbour at j+1
    alpha = da.roll(tmask,-1,axis=-2) * tmask * e3w / da.roll(e3w, -1, axis=-2)
    roll_mask = np.ones(alpha.shape, dtype=bool)
    roll_mask[:,-1,:] = False
    roll_mask = da.from_array(roll_mask)
    alpha = alpha * roll_mask

    tmp_mask = alpha > 1
    correction = correction  + tmp_mask * (1 - 1/alpha)*drhop_dk

    print(correction)

    #Calculate horizontal gradients of the denisty
    rhop_corrected = rhop_array - correction
    drhop_dy = vmask * (da.roll(rhop_corrected, -1, axis=-2) - rhop_corrected)/e2v

    #Calculate the depth's of v points
    gdepv_array = vmask*(da.roll(gdept_array*e1t*e2t, -1, axis=-2) + gdept_array*e1t*e2t )/(2*e1v*e2v)
    gdepv_array = gdepv_array + vmask * e3v / 2

    #Calculate 1/f centred on V points
    inv_ff_v = vmask * (2*e1v*e2v)/(da.roll(ff_t*e1t*e2t, -1, axis=-2) + ff_t*e1t*e2t )
    ff_v = vmask * (da.roll(ff_t*e1t*e2t, -1, axis=-2) + ff_t*e1t*e2t )/(2*e1v*e2v)

    #Combine to calculate the integrand in Sverdrups
    integrand = vmask * g * drhop_dy * gdepv_array * inv_ff_v * e3v * e2v / (1e6 * rhop_0)
    integrand = da.ma.masked_invalid(integrand)

    #Finally, integrate over depth and y
    integral = da.sum(integrand, axis=(-3,-2))

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

    acc_ssh_cube = xr.DataArray(acc_ssh.compute(),
                                dims=["time_counter", "x"],
                                coords={"time_counter":time_coord.values},
                                name='acc_ssh',
                                attrs={'units':'Sv'})

    acc_bot_cube = xr.DataArray(acc_bot.compute(),
                                dims=["time_counter", "x"],
                                coords={"time_counter":time_coord.values},
                                name='acc_bot',
                                attrs={'units':'Sv'})

    print(integrand.shape)

    integrand_cube = xr.DataArray(integrand.compute(),
                                dims=["time_counter", "model_level", "y", "x"],
                                coords={"time_counter":time_coord.values},
                                name='integrand_thw',
                                attrs={'units':'Sv'})

    integrand_top_cube = xr.DataArray( (u_top * ssh_cube.data * e2u /1e6).compute(),
                                dims=["time_counter", "y", "x"],
                                coords={"time_counter":time_coord.values},
                                name='integrand_top',
                                attrs={'units':'Sv'})

    integrand_bot_cube = xr.DataArray( (u_bot * depth_array * e2u /1e6).compute(),
                                dims=["time_counter", "y", "x"],
                                coords={"time_counter":time_coord.values},
                                name='integrand_bot',
                                attrs={'units':'Sv'})

    thw_cube = xr.DataArray( integral.compute(),
                                dims=["time_counter", "x"],
                                coords={"time_counter":time_coord.values},
                                name='thw',
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

    return [acc_ssh_cube, acc_bot_cube, thw_cube, 
            integrand_cube, integrand_top_cube,
            integrand_bot_cube, rhoS_cube, rhoN_cube,
            ff_tN_cube, ff_tS_cube, intbound_cube]


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

# def first_last_index(array, mask, axis):
#     array = da.swapaxes(array, axis, -1)
#     mask = da.swapaxes(mask, axis, -1) #Put operating axis at end of the array
#     n1, n2 = mask.shape[0], mask.shape[1]

#     inds = np.indices([mask.shape[-1]])

#     inds = np.broadcast_to(inds, mask.shape)
#     inds = np.ma.masked_array(inds, mask=~np.ma.make_mask(mask))

#     ind_1_array = np.min( inds, axis=-1)
#     ind_2_array = np.max( inds, axis=-1)

#     output_1 = np.expand_dims(np.zeros(list(array.shape[:-1])), axis=-1)
#     output_2 = np.zeros(output_1.shape)

#     for i1 in range(n1):
#         for i2 in range(n2):
#             ind_1 = ind_1_array[i1, i2]
#             ind_2 = ind_2_array[i1, i2]

#             if (np.ma.is_masked(ind_1)) or (np.ma.is_masked(ind_2)):
#                 output_1[...,i1,i2,0] = np.nan
#                 output_2[...,i1,i2,0] = np.nan
#             else: 
#                 output_1[...,i1, i2,0] = array[...,i1,i2,ind_1]
#                 output_2[...,i1, i2,0] = array[...,i1,i2,ind_2]

#     output_1 = np.ma.masked_invalid(output_1)
#     output_2 = np.ma.masked_invalid(output_2)

#     output_1 = np.swapaxes(output_1, axis, -1)
#     output_2 = np.swapaxes(output_2, axis, -1)

#     output_1 = np.squeeze(output_1, axis)
#     output_2 = np.squeeze(output_2, axis)

#     return output_1, output_2




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

