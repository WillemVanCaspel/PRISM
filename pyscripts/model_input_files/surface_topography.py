#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 18 13:16:57 2020

@author: wim
"""

import numpy as np
import pyshtools as sh
from scipy.io import FortranFile
from netCDF4 import Dataset
import pyfes
from datetime import datetime, timedelta
from matplotlib import pyplot as plt

# construct pressure grid by reading in surface pressure and using sigprim-coords
folder_surf = '/home/wim/Desktop/Projects/Model/ERA5/'
surf_f = Dataset(folder_surf + 'era5_psurf_geopsurf_mm_clim.nc')

surf_z = np.mean(surf_f['z'][:,::-2,::2],axis=0)

surf_z = np.swapaxes(surf_z,0,1)

surf_lat = surf_f['latitude'][::-2]
surf_lon = surf_f['longitude'][::2]

surf_z = np.array(surf_z,dtype=np.float64)
surf_lat = np.array(surf_lat,dtype=np.float64)
surf_lon = np.array(surf_lon,dtype=np.float64)

# max_m and max_n, corresponds to dim 60 x 120 era5 input files
max_m = np.float32(90)
max_n = np.float32(91)

# save to file
folder_output = '/home/wim/Desktop/Projects/iSSW/dope/input/'
output_name = 'surface_geopotentialtest.dat'

# write meta-data and z,t grid
f = FortranFile(folder_output + output_name,'w')

# write grid dims as integers. *2 factor to make rec len match complex len
rec_int = np.zeros(int(max_m * max_n * 2),dtype=np.int32)
rec_int[0], rec_int[1] = (max_m,
                          max_n)

f.write_record(rec_int)

rec_fl = np.zeros(int(max_m * max_n),dtype=np.float64)
rec_fl[:int(max_m)] = surf_lon
f.write_record(rec_fl)

rec_fl = np.zeros(int(max_m * max_n),dtype=np.float64)
rec_fl[:int(max_n)] = surf_lat
f.write_record(rec_fl)

f.write_record(surf_z.reshape(int(max_m * max_n),
                              order='F'))

f.close()

#%%

plt.pcolormesh(surf_z.T)
plt.colorbar()

#%%

print(24*60*60*3)

#%%

plt.plot(surf_z[0,:])


#%%
# end and start date
start_date = datetime(2012,12,1)
end_date = datetime(2013,4,1)
dt_model = (end_date - start_date).days

# time resolution (hr)
t_step = 1

# save to file
folder_output = '/home/wim/Desktop/Projects/iSSW/dope/input/'
output_name = 'FES2013_ssw_hrly_grid.dat'

# include surface topography in this file (1) or not (0)
surf_top = 0

#%%

def tide_elevation(date_start,
                   n_steps,
                   dt,
                   radial_tide=False):

    # Create handler
    short_tide = pyfes.Handler("ocean",
                               "memory",
                               ('/home/wim/anaconda3/pkgs/pyfes-2.9.2-py38h6bb024c_0/' +
                                'info/recipe/data/fes2014/ocean_tide_m2n2.ini'))
    
    load_tide = pyfes.Handler("radial",
                               "memory",
                               ('/home/wim/anaconda3/pkgs/pyfes-2.9.2-py38h6bb024c_0/' +
                                'info/recipe/data/fes2014/load_tide_m2n2.ini'))
    
    # Create grid that will be used to interpolate the tide
    latitds = np.arange(-90, 91, 2.0)
    longits = np.arange(0, 360, 4.0)
    
    full_arr = np.empty((n_steps,
                         np.shape(latitds)[0],
                         np.shape(longits)[0]))
    
    lons, lats = np.meshgrid(longits, latitds)
    
    shape = lons.shape
    
    date = date_start
    dates = np.empty(shape, dtype='datetime64[us]')
        
    for t in range(n_steps):
        dates.fill(date)
        
        # Create handler; output diurnal + semidiurnal, and long-period tides
        tide, lp, _ = short_tide.calculate(lons.ravel(), lats.ravel(),
                                           dates.ravel())
        tide, lp = tide.reshape(shape), lp.reshape(shape)
        
        if radial_tide:
            load, load_lp, _ = load_tide.calculate(lons.ravel(), lats.ravel(),
                                                   dates.ravel())
            load, load_lp = load.reshape(shape), load_lp.reshape(shape)
            
        else:
            load = np.zeros(lons.shape)
            load_lp = load
        
        # Creating an image to see the result in meters
        geo_tide = (tide + lp + load) * 0.01
        geo_tide = geo_tide.reshape(lons.shape)
        geo_tide = np.nan_to_num(geo_tide,nan=0.0) * 9.81 # geopotential height
        
        full_arr[t,:,:] = geo_tide
        
        date = date + timedelta(hours=dt)
        
        print('Creating tide day: ', t)
        
    return full_arr, latitds, longits

#% create input file

# max_m and max_n, corresponds to dim 60 x 120 era5 input files
max_m = np.float32(90)
max_n = np.float32(91)
max_l = np.float32(1) 

zero_top = np.float32(0) # dummy value for surface geop field

tdim = (dt_model + 100) * 24 # 80 days of spin-up padding
time_steps = np.zeros(tdim,dtype=np.float32)
time_steps[0] = 0 

for k in range(1,tdim):
    time_steps[k] = time_steps[k - 1] + 1 / 24
    
time_steps = np.round(time_steps,3)

tide_field, lati, longi = tide_elevation(start_date,
                                         dt_model * 24,
                                         t_step,
                                         radial_tide=True)

lati = np.array(lati,dtype=np.float64)
longi = np.array(longi,dtype=np.float64)

#%%

# write meta-data and z,t grid
f = FortranFile(folder_output + output_name,'w')

# write grid dims as integers. *2 factor to make rec len match complex len
rec_int = np.zeros(int(max_m * max_n * max_l * 2),dtype=np.int32)
rec_int[0], rec_int[1], rec_int[2], rec_int[3], rec_int[4] = (tdim,
                                                              max_m,
                                                              max_n,
                                                              max_l,
                                                              zero_top)

f.write_record(rec_int)

# write time and altitude array as floats
rec_fl = np.zeros(int(max_m * max_n * max_l * 2),dtype=np.float32)
rec_fl[:tdim] = time_steps
f.write_record(rec_fl)

rec_fl = np.zeros(int(max_m * max_n * max_l),dtype=np.float64)
rec_fl[:int(max_m)] = longi
f.write_record(rec_fl)

rec_fl = np.zeros(int(max_m * max_n * max_l),dtype=np.float64)
rec_fl[:int(max_n)] = lati
f.write_record(rec_fl)

dt, dy, dx = np.shape(tide_field)

field_ext = np.empty((dt + 100 * 24, dy, dx))

for u in range(80 * 24):
    field_ext[u,:,:] = tide_field[u % 24,:,:]

field_ext[(80 * 24):((80 + dt_model) * 24),:,:] = tide_field

for u in range((80 + dt_model) * 24,(100 + dt_model) * 24):
    field_ext[u,:,:] = tide_field[((dt - 24 + u % 24)),:,:]
    
del tide_field

# switch to sigprim [lon,lat] configuration
field_ext = np.swapaxes(field_ext,1,2)
field_ext = np.array(field_ext,dtype=np.float64)
    
for t in range(np.int(tdim)):
    
    temp_coeffs = field_ext[t,:,:]
        
    f.write_record(temp_coeffs.reshape(int(max_m * max_n * max_l),
                                            order='F'))
    
    print('Writing tide field: ',t)
    
f.close()

del field_ext

#%%