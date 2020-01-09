import xarray as xr
import numpy as np
from utilities import stretching, set_depth
 
def retrieve_ROMS_data(thredds_url: str):
    """Return an xarray DataSet from TDS server
    
    Arguments:
        thredds_url {str} -- url to the THREDDS server where model data is located   
    
    Return:
        roms_ds {xarray.DataSet} -- xarray.DataSet object that is (likely) lazy loaded
    """
    
    try:
        roms_ds = xr.open_dataset(thredds_url, decode_cf=False)
    except OSError:
        print(f"Could Not Access {thredds_url}") # Switch to log
        roms_ds = None
    return(roms_ds)

def calculate_discrete_depths(roms_data: xr.Dataset):
    """Calculate the discrete depths from the sigma-coordinates. This relies on code from Fred B. See simga_coord_utilities
    
    Arguments:
        roms_data {xr.Dataset} -- xarray datasets lazy loaded from TDS
    
    Returns:
        [np.ndarray] -- discrete depths for every locations (lat, lon, z)
    """
    vertical_stretching_function = roms_data['Vstretching'] # vertical terrain following stretching function
    vertical_transformation_function = roms_data['Vtransform'] # vertical terrain following tranformation equation
    h = roms_data['h'] # bathymetry at rho points
    h_critical = roms_data['hc'] # S-coordinate parameter, critical depth
    theta_bottom = roms_data['theta_b'].values # S-coordinate bottom control parameter
    theta_surface = roms_data['theta_s'].values # S-coordinate surface control parameter
    zeta = roms_data['zeta'] # free surface
    N_levels=42 # number of levels, note w has 43 levels
    zeta = np.squeeze(zeta[0,:,:]) # reduce dimensions to

    igrid=1 # density point grid, for (T,S) use igrid=3, for u and igrid=4 for v # use igrid=5 for w
    z_rho = set_depth.set_depth(vertical_transformation_function, vertical_stretching_function, theta_surface, theta_bottom, h_critical, N_levels, igrid, h, zeta, report=0)
    return z_rho

def trim_dataset(roms_data: xr.Dataset, z_rho: np.ndarray):
    """ Reduce the number of variables in the dataset to just track temperature and make the discrete depths a dimension.
    
    Arguments:
        roms_data {xr.Dataset}
        z_rho {np.ndarray}
    
    Returns:
        xr.Dataset -- Dataset with only Temperature data of dimensions (time, latitude, longitude, s_rho)
    """    
    roms_temp = roms_data['temp']
    roms_temp_latitude = roms_data['lat_rho'].values[:,0]
    roms_temp_longitude = roms_data['lon_rho'].values[0,:]
    roms_temp = roms_temp.assign_coords({"eta_rho":roms_temp_latitude, "xi_rho":roms_temp_longitude})
    roms_temp = roms_temp.rename({"eta_rho":"latitude", "xi_rho":"longitude"})
    roms_temp = roms_temp.drop('lat_rho', errors='ignore')
    roms_temp = roms_temp.drop('lon_rho', errors='ignore')
    roms_temp = roms_temp.drop('time_run', errors='ignore')
    roms_temp = roms_temp.transpose('time','latitude','longitude','s_rho')
    roms_temp['z_rho'] = (('latitude','longitude','s_rho'),z_rho) # Add the depth values back to the Dataset
    return roms_temp


def get_roms():
    """ Use this as the main function, should probably be refactored """
    url = "http://oceanmodeling.pmc.ucsc.edu:8080/thredds/dodsC/ccsra_2016a_phys_agg_slevs/fmrc/CCSRA_2016a_Phys_ROMS_Sigma-level_Aggregation_best.ncd"
    ds = retrieve_ROMS_data(url)
    z_rho = calculate_discrete_depths(ds)
    ds = trim_dataset(ds, z_rho)
    return ds