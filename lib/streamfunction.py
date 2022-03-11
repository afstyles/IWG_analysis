#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# streamfunction.py
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

def DepIntSf(data_list, mask_list, var_dict, WG_bounds=None ):
    """
    Calculate the stream function for the depth-integrated flow
    """
    from cubeprep import CubeListExtract as CLE
    import iris
    from iris.cube import Cube
    from iris.coords import AuxCoord
    import numpy as np

    u_cube = CLE(data_list, var_dict['u']) 

    umask = iris.util.squeeze(CLE(mask_list, var_dict['umask'])).data
    e1u = iris.util.squeeze(CLE(mask_list, var_dict['e1u'])).data
    e2u = iris.util.squeeze(CLE(mask_list, var_dict['e2u'])).data
    e3u = iris.util.squeeze(CLE(mask_list, var_dict['e3u'])).data


    u = np.ma.masked_array(u_cube.data, mask=np.broadcast_to(~np.ma.make_mask(umask),u_cube.shape))
    u_zint = np.sum(u*e3u, axis=-3)

    integrand = -u_zint*e2u

    sf_zint = np.cumsum(integrand, axis=-2)


    #We need to remove the non-zero values at [j=0]
    # sf_zint = np.swapaxes(sf_zint, -2, 0) #Temporarily put the j axis at the end for easier broadcasting
    # sf_zint = sf_zint - np.broadcast_to(sf_zint[0,...],sf_zint.shape)
    # sf_zint = np.swapaxes(sf_zint, -2, 0) # Restore the original axis order


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


def ZonIntSF(data_list, mask_list, var_dict):
    """
    Calculate the stream function of the zonally integrated flow
    """
    import numpy as np
    import iris
    from iris.cube import Cube
    from iris.coords import AuxCoord, DimCoord
    from cubeprep import CubeListExtract as CLE

    v_cube = CLE(data_list, var_dict['v']) 

    vmask = iris.util.squeeze(CLE(mask_list, var_dict['vmask'])).data
    vmask2d = iris.util.squeeze(CLE(mask_list, var_dict['vmask2d'])).data

    e1v = iris.util.squeeze(CLE(mask_list, var_dict['e1v'])).data
    e2v = iris.util.squeeze(CLE(mask_list, var_dict['e2v'])).data
    e3v = iris.util.squeeze(CLE(mask_list, var_dict['e3v'])).data

    v = np.ma.masked_array(v_cube.data, mask=np.broadcast_to(~np.ma.make_mask(vmask),v_cube.shape))
    vflux_xint = np.sum(v*e3v*e1v, axis=-1)

    sf_xint = np.cumsum(np.flip(vflux_xint,axis=-2), axis=-2) #We flip the depth axis so we can integrate from the bottom 

    #We need to remove the non-zero values at [k=0]
    # sf_xint = np.swapaxes(sf_xint, -2, 0) #Temporarily put the j axis at the end for easier broadcasting
    # sf_xint = sf_xint - np.broadcast_to(sf_xint[0,...],sf_xint.shape)
    # sf_xint = np.swapaxes(sf_xint, -2, 0) # Restore the original axis order

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

def ResidualOverturning(data_list, mask_list, nn_rhop, var_dict):
    """
    Calculate the stream function of the zonally integrated flow
    """
    import numpy as np
    import iris
    from iris.cube import Cube
    from iris.coords import AuxCoord, DimCoord
    from cubeprep import CubeListExtract as CLE
    from fortran_lib import fortran_lib



    v_cube = CLE(data_list, var_dict['v'])
    rhop_cube = CLE(data_list, var_dict['rho'] )

    tmask = iris.util.squeeze(CLE(mask_list, var_dict['tmask'])).data
    vmask = iris.util.squeeze(CLE(mask_list, var_dict['vmask'])).data
    # vmask2d = iris.util.squeeze(CLE(mask_list, var_dict['vmask2d'])).data

    e1v = iris.util.squeeze(CLE(mask_list, var_dict['e1v'])).data
    e2v = iris.util.squeeze(CLE(mask_list, var_dict['e2v'])).data
    e3v = iris.util.squeeze(CLE(mask_list, var_dict['e3v'])).data

    e1t = iris.util.squeeze(CLE(mask_list, var_dict['e1v'])).data
    e2t = iris.util.squeeze(CLE(mask_list, var_dict['e2v'])).data
    e3t = iris.util.squeeze(CLE(mask_list, var_dict['e3v'])).data
    # deptht_cube = iris.util.squeeze(CLE(mask_list, var_dict['deptht']))

    # #We need to calculate the density and depth centred on the V points
    rhop_v = (tmask*rhop_cube.data*e1t*e2t + jp1(tmask*rhop_cube.data*e1t*e2t))/((tmask+jp1(tmask))*e1v*e2v)  
    rhop_v = np.ma.masked_array( rhop_v, mask = np.broadcast_to(~np.ma.make_mask(vmask), rhop_v.shape))

    # depthv = (tmask*deptht_cube.data*e1t*e2t + jp1(tmask*deptht_cube.data*e1t*e2t))/((tmask+jp1(tmask))*e1v*e2v)  
    # depthv = np.ma.masked_array( depthv, mask = ~np.ma.make_mask(vmask))

    rhop_min = np.min(rhop_v)
    rhop_max = np.max(rhop_v)
    rhop_coord = np.linspace(rhop_min, rhop_max, num=nn_rhop)
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



def jp1(M): 
    import numpy as np
    output = np.roll(M, -1, axis=-2)
    output[...,-1,:] = 0.

    return output





