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
import iris
from iris.cube import Cube
from iris.coords import AuxCoord, DimCoord
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

    temp_cube = CLE(data_list, var_dict['temp'])
    sal_cube = CLE(data_list, var_dict['sal'])
    rhop_cube = CLE(data_list, var_dict['rho'])

    x_t = np.mean(CLE(mask_list, var_dict['x']).data, axis=-2)
    y_t = np.mean(CLE(mask_list, var_dict['y']).data, axis=-1)
    e1t = iris.util.squeeze(CLE(mask_list, var_dict['e1v'])).data
    e2t = iris.util.squeeze(CLE(mask_list, var_dict['e2v'])).data
    e3t = iris.util.squeeze(CLE(mask_list, var_dict['e3v'])).data
    tmask = iris.util.squeeze(CLE(mask_list, var_dict['tmask'])).data

    volume_array = e1t * e2t * e3t * tmask
    volume_array = np.ma.masked_array(volume_array, mask=~np.ma.make_mask(tmask))

    temp_array = temp_cube.data * volume_array
    sal_array = sal_cube.data * volume_array
    rhop_array = rhop_cube.data * volume_array

    try: 
        time_coord = temp_cube.coord("time")
    except:
        aux_time = temp_cube.aux_coords[0]
        aux_time.rename("aux_time")
        time_coord = temp_cube.coord("time")


    temp_xint_cube_list = []
    sal_xint_cube_list  = []
    rhop_xint_cube_list = []


    for ir in range(nranges):
        
        if ir <= len(xmin_list) - 1:
            xmin = xmin_list[ir]
            xmax = xmax_list[ir]
            label = range_labels[ir]

            xind_min = np.argmin(np.abs(x_t - xmin))
            xind_max = np.argmin(np.abs(x_t - xmax))

            if ( x_t[xind_min] - xmin < 0 ): xind_min = xind_min + 1
            if ( x_t[xind_max] - xmax > 0 ): xind_max = xind_max - 1

        else:
            xind_min = 0
            xind_max = len(x_t) + 1
            xmin = x_t[xind_min]
            xmax = x_t[-1]
            label = "fullwidth"

        temp_xint = np.sum( temp_array[...,xind_min:xind_max+1], axis=-1) / np.sum(volume_array[...,xind_min:xind_max+1], axis=-1)
        sal_xint = np.sum( sal_array[...,xind_min:xind_max+1], axis=-1) / np.sum(volume_array[...,xind_min:xind_max+1], axis=-1)
        rhop_xint = np.sum( rhop_array[...,xind_min:xind_max+1], axis=-1) / np.sum(volume_array[...,xind_min:xind_max+1], axis=-1)


        temp_xint_cube = Cube(temp_xint, dim_coords_and_dims=[(time_coord,0)])
        temp_xint_cube.long_name = "Mean Temperature over Range " + label
        temp_xint_cube.var_name = "temp_xint_" + label
        temp_xint_cube.units = temp_cube.units
        temp_xint_cube.attributes = {'xmin':xmin, 'xmax':xmax, 'xind_min':xind_min, 'xind_max':xind_max}
        
        sal_xint_cube = Cube(sal_xint, dim_coords_and_dims=[(time_coord,0)])
        sal_xint_cube.long_name = "Mean Salinity over Range: "+ label
        sal_xint_cube.var_name = "sal_xint_" +  label
        sal_xint_cube.units = sal_cube.units
        sal_xint_cube.attributes = {'xmin':xmin, 'xmax':xmax, 'xind_min':xind_min, 'xind_max':xind_max}

        rhop_xint_cube = Cube(rhop_xint, dim_coords_and_dims=[(time_coord,0)])
        rhop_xint_cube.long_name = "Mean Density over Range: "+ label
        rhop_xint_cube.var_name = "rhop_xint_"+  label 
        rhop_xint_cube.units = rhop_cube.units
        rhop_xint_cube.attributes = {'xmin':xmin, 'xmax':xmax, 'xind_min':xind_min, 'xind_max':xind_max}

        temp_xint_cube_list.append(temp_xint_cube)
        sal_xint_cube_list.append(sal_xint_cube)
        rhop_xint_cube_list.append(rhop_xint_cube)

    return temp_xint_cube_list, sal_xint_cube_list, rhop_xint_cube_list