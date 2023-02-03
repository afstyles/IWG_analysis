"""
            _                                            
           | |                                           
  ___ _   _| |__   ___ _ __  _ __ ___ _ __   _ __  _   _ 
 / __| | | | '_ \ / _ \ '_ \| '__/ _ \ '_ \ | '_ \| | | |
| (__| |_| | |_) |  __/ |_) | | |  __/ |_) || |_) | |_| |
 \___|\__,_|_.__/ \___| .__/|_|  \___| .__(_) .__/ \__, |
                      | |            | |    | |     __/ |
                      |_|            |_|    |_|    |___/ 
cubeprep.py

Routines and information for extracting data from IRIS cubelists

Contains methods:
    CubeListExtract --> Extract a cube from a CubeList object using the variable name
    var_names --> Returns a dictionary relating generic names of variables to specific variable names for the data.
"""

def var_names():
    """
    Load the variable names for cubes that might be used in the analysis
    """
    varname_dict = { 'u'  : 'uoce' ,     # x velocity
                 'v'  : 'voce' ,         # y velocity
                 'w'  : 'woce' ,         # z velocity
                 'rho': 'swsigthet' ,    # Density
                 'temp': 'toce',         # Temperature
                 'sal' : 'soce',         # Salinity
                 'ssh': 'ssh'  ,         # Sea surface height 
                 'e1u': 'e1u' ,          # U cell width     (i direction)
                 'e2u': 'e2u' ,          # U cell width     (j direction)
                 'e3u': 'e3u_0',         # U cell thickness (k direction)
                 'e1v': 'e1v' ,          # V cell width     (i direction)
                 'e2v': 'e2v' ,          # V cell width     (j direction)
                 'e3v': 'e3v_0',         # V cell thickness (k direction)
                 'e1t': 'e1t' ,          # T cell width     (i direction)
                 'e2t': 'e2t' ,          # T cell width     (j direction)
                 'e3t': 'e3t_0',         # T cell thickness (k direction)
                 'e1f': 'e1f' ,          # F cell width     (i direction) 
                 'e2f': 'e2f' ,          # F cell width     (j direction)
                 'e3w': 'e3w_0' ,        # W cell thickness (k direction)
                 'deptht': 'gdept_0',    # T cell depth 
                 'umask': 'umask',       # U cell mask 
                 'vmask': 'vmask',       # V cell mask
                 'tmask': 'tmask',       # T cell mask
                 'umask2d': 'umaskutil', # U cell mask (2d fields)
                 'vmask2d': 'vmaskutil', # V cell mask (2d fields)
                 'tmask2d': 'tmaskutil', # T cell mask (2d fields)
                 'x': 'nav_lon',       # X coordinates of T cell (km)
                 'y': 'nav_lat',       # Y coordinates of T cell (km)
                 'ff_f': 'ff_f',          # Coriolis parameter centred on F point (1/s)
                 'ff_t': 'ff_t',          # Coriolis parameter centred on T point (1/s)
                }

    return varname_dict

def CubeListExtract(cube_list, var_name):
    """
    Extracts a cube with a specific variable name from a cube list. 
    If two cubes have the same variable name in the list then the first occurence in the list is extracted
    cube_list - CubeList object
    var_name - String, variable name of cube to extract
    Returns
    cube - IRIS cube with variable name matching var_name
    """
    VarNames = [cube_list[i].var_name for i in range(len(cube_list))]
    
    try:
    
        index = VarNames.index(var_name)
    
        return cube_list[index]
    
    except:
        
        print('Variable name not found: ', var_name)
        
        return 
