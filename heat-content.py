import numpy as np
import xarray as xr
import get_4DVAR_roms
from scipy import integrate
import time
import datetime as dt

def calculate_upper_heat_content(temperature: xr.Dataset, upper_depth=-100):
    """ 
    Calculate upper ocean heat content from ROMS
    
    Calculate the upper ocean heat (J/m^2) content from the CA ROMS model output
    Density is assumed constant at 1025 kg/m^3
    specific heat is 3850 J/(kg C)
    Arguments:
        temperature {xr.Dataset} -- xarray dataset with dimensions (latitude,longitude,s_sho) 
        upper_depth {int} : depth to integrate to from the surface, values below the surface are negative {default = -100}

    Returns:
        upper_heat_content: array of upper heat content values in the shape of the temperature input dimensions (latitude, longitude)
        upper_heat_content_b_bins: array of the number of depth bins used for integrated upper heat content in the shape of the temperature input dimensions (latitude, longitude)
    """
    z_rho = temperature['z_rho'].values
    cp = 3850
    density = 1025
    
    if len(z_rho.shape) == 1:
        '''Calculating at single point'''
        upper_heat_content = np.zeros(shape=(1))
        roms_temperature_flattened = temperature.values.reshape(-1, temperature.shape[-1])
        
    else:
        upper_heat_content = np.zeros(shape=(z_rho.shape[0]*z_rho.shape[1]))
        roms_temperature_flattened = temperature.values.reshape(-1, temperature.shape[-1])
        
    for i, grid_point in enumerate(z_rho.reshape(-1, z_rho.shape[-1])):
        #Find index of depth less than 100 meters
        upper_ix = np.where(grid_point >= upper_depth)[0]
        # Integrate
        if not upper_ix.size == 0:
            
            interp_top = np.interp(0,grid_point,roms_temperature_flattened[i,:])
            interp_bottom = np.interp(-100, grid_point,roms_temperature_flattened[i,:])
            # Add interpolated values to temperature array for integrtaion
            upper_values = roms_temperature_flattened[i,upper_ix]
            upper_values = np.append(interp_bottom,upper_values)
            upper_values = np.append(upper_values,interp_top)
            # Add interpolated depths to depth array for integrtaion
            upper_depths = grid_point[upper_ix]
            upper_depths = np.append(-100,upper_depths)
            upper_depths = np.append(upper_depths,0)
            integrated_temp = integrate.simps(upper_values,upper_depths)                    
            upper_heat_content[i] = density * cp * integrated_temp
            
        else:
            upper_heat_content[i] = np.nan
    
    if len(z_rho.shape) == 1:
        return upper_heat_content
    
    else:
        upper_heat_content = upper_heat_content.reshape(z_rho.shape[:2])
        return upper_heat_content


def all_heat_content():
    """Initial run, runs over the entire spatial and temporal domain

    Arguments:
        model_data {xr.Dataset} -- [description]
    """
    model_data = get_4DVAR_roms.get_roms()
    heat_content = np.zeros(shape=model_data[:,:,:,0].shape)
    start_time = time.time()
    # for i in range(len(heat_content))[30]:
    for i in range(20):
        upper_heat_content = calculate_upper_heat_content(model_data.isel(time=i), upper_depth=-100)
        heat_content[i,:,:] = upper_heat_content
        if i%10 == 0:
            print(i,round(time.time()-start_time,1),'s')
            start_time = time.time()
    
    # write data out 
    model_dates = model_data.time.values
    longitude = model_data.longitude.values
    latitude = model_data.latitude.values
    dims = ['time', 'lat', 'lon']
    ds = xr.Dataset({'heat_content_100_meters': (dims, heat_content)},
                coords={"time": model_dates,
                        "lon": longitude,
                        "lat": latitude,
                        })
    ds.attrs['title'] = "West Coast Upper Ocean Heat Content - 0-100 meters"
    ds.attrs['notes'] = "Created on "+dt.datetime.today().strftime("%Y-%m-%d") + " by pdaniel"
    fname = "data/4DVAR-ROMS-heat-content.nc"
    ds.to_netcdf(path=fname)


all_heat_content()