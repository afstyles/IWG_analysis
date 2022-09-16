"""
 _                                             
| |                                            
| |_ _ __ __ _  ___ ___ _ __ ___   _ __  _   _ 
| __| '__/ _` |/ __/ _ \ '__/ __| | '_ \| | | |
| |_| | | (_| | (_|  __/ |  \__ \_| |_) | |_| |
 \__|_|  \__,_|\___\___|_|  |___(_) .__/ \__, |
                                  | |     __/ |
                                  |_|    |___/ 
tracers.py

Manipulation and integration of tracer fields

Contains methods:
    tracer_xint --> Calculates the zonal mean of temperature, salinity, and density over a list of ranges.
"""

from cubeprep import CubeListExtract as CLE
# import iris
# from iris.cube import Cube
# from iris.coords import AuxCoord, DimCoord
import dask
import dask.array as da
import xarray as xr
import numpy as np

def tracer_xint(data_list, mask_list, var_dict, xmin_list, xmax_list, range_labels):
    """
    Calculates the zonal mean of temperature, salinity, and density over a list of ranges.

    INPUT variables
    data_list - CubeList object of data from NEMO run
    mask_list - CubeList object of data from mesh_mask
    var_dict - Dictionary of variable names for needed data
    xmin_list - List or tuple of floats: minimum x values for n ranges (n) [km]
    xmax_list - List or tuple of floats: maximum x values for n ranges (n) [km]
    range_labels - Tuple of strings to label the n ranges              (n) 

    OUTPUT variables
    temp_xint_cube_list - CubeList object: Zonal mean of temperature over the specified ranges and the whole domain (t,z,y)
    sal_xint_cube_list - CubeList object: Zonal mean of salinity over the specified ranges and the whole domain (t,z,y)
    rhop_xint_cube_list - CubeList object: Zonal mean of density over the specified ranges and the whole domain (t,z,y)
    """

    nranges = len(xmin_list) + 1

    temp_cube = data_list[var_dict['temp']]
    sal_cube = data_list[var_dict['sal']]
    rhop_cube = data_list[var_dict['rho']]

    time_coord = temp_cube.coords["time_counter"]

    x_t = da.squeeze(mask_list[var_dict['x']].data).mean(axis=-2)
    y_t = da.squeeze(mask_list[var_dict['y']].data).mean(axis=-1)
    e1t = da.squeeze(mask_list[var_dict['e1t']].data)
    e2t = da.squeeze(mask_list[var_dict['e1t']].data)
    e3t = da.squeeze(mask_list[var_dict['e1t']].data)

    tmask = da.squeeze(mask_list[var_dict['tmask']].data).astype(bool)

    volume_array = e1t * e2t * e3t * tmask

    volume_array = da.ma.masked_array(volume_array, mask=da.broadcast_to(~tmask, volume_array.shape))

    temp_array = temp_cube.data * volume_array
    sal_array = sal_cube.data * volume_array
    rhop_array = rhop_cube.data * volume_array

    temp_xint_cube_list = []
    sal_xint_cube_list  = []
    rhop_xint_cube_list = []


    for ir in range(nranges):
        
        if ir <= len(xmin_list) - 1:
            xmin = xmin_list[ir]
            xmax = xmax_list[ir]
            label = range_labels[ir]

            tmp_mask = ( da.ma.getmaskarray(da.ma.masked_greater(x_t, xmin))
                       * da.ma.getmaskarray(da.ma.masked_less(x_t, xmax)))

            xind_min = da.argmin(np.abs(x_t - xmin)).compute()
            xind_max = da.argmin(np.abs(x_t - xmax)).compute()

            if ( x_t[xind_min] - xmin < 0 ): xind_min = xind_min + 1
            if ( x_t[xind_max] - xmax > 0 ): xind_max = xind_max - 1

        else:
            xind_min = 0
            xind_max = len(x_t) + 1
            xmin = x_t.compute()[xind_min]
            xmax = x_t.compute()[-1]
            label = "fullwidth"
            tmp_mask = True

        vol_xint  = da.sum( volume_array * tmp_mask, axis=-1 )
        vol_xint = da.ma.masked_less_equal(vol_xint, 0)
        temp_xint = da.sum( temp_array   * tmp_mask, axis=-1  ) / vol_xint
        sal_xint  = da.sum(  sal_array   * tmp_mask, axis=-1  ) / vol_xint
        rhop_xint = da.sum( rhop_array   * tmp_mask, axis=-1  ) / vol_xint

        temp_xint_cube = xr.DataArray(temp_xint.compute(),
                                      dims=["time_counter", "model_level", "y"],
                                      coords={"time_counter":time_coord.values, "y":y_t.compute()},
                                      name="temp_xint_" + label,
                                      attrs={'units':temp_cube.attrs['units'],
                                             'xmin':xmin, 'xmax':xmax,
                                             'xind_min':xind_min, 'xind_max':xind_max})

        sal_xint_cube = xr.DataArray(sal_xint.compute(),
                                      dims=["time_counter", "model_level", "y"],
                                      coords={"time_counter":time_coord.values, "y":y_t.compute()},
                                      name="sal_xint_" + label,
                                      attrs={'units':sal_cube.attrs['units'],
                                             'xmin':xmin, 'xmax':xmax,
                                             'xind_min':xind_min, 'xind_max':xind_max})

        rhop_xint_cube = xr.DataArray(rhop_xint.compute(),
                                      dims=["time_counter", "model_level", "y"],
                                      coords={"time_counter":time_coord.values, "y":y_t.compute()},
                                      name="rhop_xint_" + label,
                                      attrs={'units':rhop_cube.attrs['units'],
                                             'xmin':xmin, 'xmax':xmax,
                                             'xind_min':xind_min, 'xind_max':xind_max})

        temp_xint_cube_list.append(temp_xint_cube)
        sal_xint_cube_list.append(sal_xint_cube)
        rhop_xint_cube_list.append(rhop_xint_cube)

    return temp_xint_cube_list, sal_xint_cube_list, rhop_xint_cube_list

def uvar_xint(data_list, mask_list, var_dict, xmin_list, xmax_list, range_labels):
    """
    Calculates the zonal mean of x-velocity over a list of ranges.

    INPUT variables
    data_list - CubeList object of data from NEMO run
    mask_list - CubeList object of data from mesh_mask
    var_dict - Dictionary of variable names for needed data
    xmin_list - List or tuple of floats: minimum x values for n ranges (n) [km]
    xmax_list - List or tuple of floats: maximum x values for n ranges (n) [km]
    range_labels - Tuple of strings to label the n ranges              (n) 

    OUTPUT variables
    u_xint_cube_list - CubeList object: Zonal mean of x-velocity over the specified ranges and the whole domain (t,z,y)
    """

    nranges = len(xmin_list) + 1

    u_cube = data_list[var_dict['u']]
    time_coord = u_cube.coords["time_counter"]

    x_t = da.squeeze(mask_list[var_dict['x']].data).mean(axis=-2)
    y_t = da.squeeze(mask_list[var_dict['y']].data).mean(axis=-1)
    e1u = da.squeeze(mask_list[var_dict['e1u']].data)
    e2u = da.squeeze(mask_list[var_dict['e2u']].data)
    e3u = da.squeeze(mask_list[var_dict['e3u']].data)
    umask = da.squeeze(mask_list[var_dict['umask']].data).astype(bool)

    volume_array = e1u * e2u * e3u * umask
    volume_array = da.ma.masked_array(volume_array, mask=da.broadcast_to(~umask, volume_array.shape))

    u_array = u_cube.data * volume_array
    u_xint_cube_list = []

    for ir in range(nranges):
        
        if ir <= len(xmin_list) - 1:
            xmin = xmin_list[ir]
            xmax = xmax_list[ir]
            label = range_labels[ir]

            tmp_mask = ( da.ma.getmaskarray(da.ma.masked_greater(x_t, xmin))
                       * da.ma.getmaskarray(da.ma.masked_less(x_t, xmax)))

            xind_min = da.argmin(np.abs(x_t - xmin)).compute()
            xind_max = da.argmin(np.abs(x_t - xmax)).compute()

            if ( x_t[xind_min] - xmin < 0 ): xind_min = xind_min + 1
            if ( x_t[xind_max] - xmax > 0 ): xind_max = xind_max - 1

        else:
            xind_min = 0
            xind_max = len(x_t) + 1
            xmin = x_t.compute()[xind_min]
            xmax = x_t.compute()[-1]
            label = "fullwidth"
            tmp_mask = True

        vol_xint  = da.sum( volume_array * tmp_mask, axis=-1 )
        vol_xint = da.ma.masked_less_equal(vol_xint, 0)
        u_xint = da.sum( u_array   * tmp_mask, axis=-1  ) / vol_xint

        u_xint_cube = xr.DataArray(u_xint.compute(),
                                      dims=["time_counter", "model_level", "y"],
                                      coords={"time_counter":time_coord.values, "y":y_t.compute()},
                                      name="u_xint_" + label,
                                      attrs={'units':u_cube.attrs['units'],
                                             'xmin':xmin, 'xmax':xmax,
                                             'xind_min':xind_min, 'xind_max':xind_max})

        u_xint_cube_list.append(u_xint_cube)

    return u_xint_cube_list

def vvar_xint(data_list, mask_list, var_dict, xmin_list, xmax_list, range_labels):
    """
    Calculates the zonal mean of y-velocity over a list of ranges.

    INPUT variables
    data_list - CubeList object of data from NEMO run
    mask_list - CubeList object of data from mesh_mask
    var_dict - Dictionary of variable names for needed data
    xmin_list - List or tuple of floats: minimum x values for n ranges (n) [km]
    xmax_list - List or tuple of floats: maximum x values for n ranges (n) [km]
    range_labels - Tuple of strings to label the n ranges              (n) 

    OUTPUT variables
    v_xint_cube_list - CubeList object: Zonal mean of x-velocity over the specified ranges and the whole domain (t,z,y)
    """

    nranges = len(xmin_list) + 1

    v_cube = data_list[var_dict['v']]
    time_coord = v_cube.coords["time_counter"]

    x_t = da.squeeze(mask_list[var_dict['x']].data).mean(axis=-2)
    y_t = da.squeeze(mask_list[var_dict['y']].data).mean(axis=-1)
    e1v = da.squeeze(mask_list[var_dict['e1v']].data)
    e2v = da.squeeze(mask_list[var_dict['e2v']].data)
    e3v = da.squeeze(mask_list[var_dict['e3v']].data)
    vmask = da.squeeze(mask_list[var_dict['vmask']].data).astype(bool)

    volume_array = e1v * e2v * e3v * vmask
    volume_array = da.ma.masked_array(volume_array, mask=da.broadcast_to(~vmask, volume_array.shape))

    v_array = v_cube.data * volume_array
    v_xint_cube_list = []

    for ir in range(nranges):
        
        if ir <= len(xmin_list) - 1:
            xmin = xmin_list[ir]
            xmax = xmax_list[ir]
            label = range_labels[ir]

            tmp_mask = ( da.ma.getmaskarray(da.ma.masked_greater(x_t, xmin))
                       * da.ma.getmaskarray(da.ma.masked_less(x_t, xmax)))

            xind_min = da.argmin(np.abs(x_t - xmin)).compute()
            xind_max = da.argmin(np.abs(x_t - xmax)).compute()

            if ( x_t[xind_min] - xmin < 0 ): xind_min = xind_min + 1
            if ( x_t[xind_max] - xmax > 0 ): xind_max = xind_max - 1

        else:
            xind_min = 0
            xind_max = len(x_t) + 1
            xmin = x_t.compute()[xind_min]
            xmax = x_t.compute()[-1]
            label = "fullwidth"
            tmp_mask = True

        
        vol_xint  = da.sum( volume_array * tmp_mask, axis=-1 )
        vol_xint = da.ma.masked_less_equal(vol_xint, 0)
        v_xint = da.sum( v_array   * tmp_mask, axis=-1  ) / vol_xint

        v_xint_cube = xr.DataArray(v_xint.compute(),
                                      dims=["time_counter", "model_level", "y"],
                                      coords={"time_counter":time_coord.values, "y":y_t.compute()},
                                      name="v_xint_" + label,
                                      attrs={'units':v_cube.attrs['units'],
                                             'xmin':xmin, 'xmax':xmax,
                                             'xind_min':xind_min, 'xind_max':xind_max})

        v_xint_cube_list.append(v_xint_cube)

    return v_xint_cube_list

def eke_xint(eke_cube, mask_list, var_dict, xmin_list, xmax_list, range_labels):
    """
    Calculates the zonal mean of eddy kinetic energy over a list of ranges.

    INPUT variables
    eke_cube - IRIS cube of eddy kinetic energy, calculated by eddy_energy.eddy_kinetic_energy()
    var_dict - Dictionary of variable names for needed data
    xmin_list - List or tuple of floats: minimum x values for n ranges (n) [km]
    xmax_list - List or tuple of floats: maximum x values for n ranges (n) [km]
    range_labels - Tuple of strings to label the n ranges              (n) 

    OUTPUT variables
    v_xint_cube_list - CubeList object: Zonal mean of x-velocity over the specified ranges and the whole domain (t,z,y)
    """

    nranges = len(xmin_list) + 1

    x_t = da.squeeze(mask_list[var_dict['x']].data).mean(axis=-2)
    y_t = da.squeeze(mask_list[var_dict['y']].data).mean(axis=-1)
    e1t = da.squeeze(mask_list[var_dict['e1t']].data)
    e2t = da.squeeze(mask_list[var_dict['e1t']].data)
    e3t = da.squeeze(mask_list[var_dict['e1t']].data)

    tmask = da.squeeze(mask_list[var_dict['tmask']].data).astype(bool)

    volume_array = e1t * e2t * e3t * tmask

    volume_array = da.ma.masked_array(volume_array, mask=da.broadcast_to(~tmask, volume_array.shape))

    eke_array = eke_cube.data * volume_array
    eke_xint_cube_list = []

    for ir in range(nranges):
        
        if ir <= len(xmin_list) - 1:
            xmin = xmin_list[ir]
            xmax = xmax_list[ir]
            label = range_labels[ir]

            tmp_mask = ( da.ma.getmaskarray(da.ma.masked_greater(x_t, xmin))
                       * da.ma.getmaskarray(da.ma.masked_less(x_t, xmax)))

            xind_min = da.argmin(np.abs(x_t - xmin)).compute()
            xind_max = da.argmin(np.abs(x_t - xmax)).compute()

            if ( x_t[xind_min] - xmin < 0 ): xind_min = xind_min + 1
            if ( x_t[xind_max] - xmax > 0 ): xind_max = xind_max - 1

        else:
            xind_min = 0
            xind_max = len(x_t) + 1
            xmin = x_t.compute()[xind_min]
            xmax = x_t.compute()[-1]
            label = "fullwidth"
            tmp_mask = True

        vol_xint  = da.sum( volume_array * tmp_mask, axis=-1 )
        vol_xint = da.ma.masked_less_equal(vol_xint, 0)
        eke_xint = da.sum( eke_array   * tmp_mask, axis=-1  ) / vol_xint

        eke_xint_cube = xr.DataArray(eke_xint.compute(),
                                      dims=["model_level", "y"],
                                      coords={"y":y_t.compute()},
                                      name="eke_xint_" + label,
                                      attrs={'units':eke_cube.attrs['units'],
                                             'xmin':xmin, 'xmax':xmax,
                                             'xind_min':xind_min, 'xind_max':xind_max})


        eke_xint_cube_list.append(eke_xint_cube)

    return eke_xint_cube_list

