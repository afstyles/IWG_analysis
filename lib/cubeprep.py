# >>>>>>>>>>>>>>>>>>>>>>>>>>>>
# cubeprep.py
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>

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

def var_names():
    """
    Load the variable names for cubes that might be used in the analysis
    """
    varname_dict = { 'u'  : 'uoce' ,     #x velocity
                 'v'  : 'voce' ,     #y velocity
                 'w'  : 'woce' ,     #z velocity
                 'rho': 'swsigthet' ,#Density
                 'ssh': 'ssh'  ,     # ssh 
                 'e1u': 'e1u' ,      
                 'e2u': 'e2u' ,  
                 'e3u': 'e3u_0',     
                 'e1v': 'e1v' ,       
                 'e2v': 'e2v' ,
                 'e3v': 'e3v_0',       
                 'e1t': 'e1t' ,       
                 'e2t': 'e2t' ,
                 'e3t': 'e3t_0',
                 'e1f': 'e1f' ,      
                 'e2f': 'e2f' ,
                 'e3w': 'e3w_0' ,
                 'deptht': 'gdept_0',
                 'umask': 'umask',    
                 'vmask': 'vmask',
                 'tmask': 'tmask',  
                 'umask2d': 'umaskutil',    
                 'vmask2d': 'vmaskutil',
                 'tmask2d': 'tmaskutil', 
                 'lon': 'nav_lon',
                 'lat': 'nav_lat', 
                }

    return varname_dict
