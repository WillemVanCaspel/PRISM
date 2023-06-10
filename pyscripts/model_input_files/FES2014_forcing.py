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

for year in range(2014,2018):
    
    print(year)
    
    # end and start date
    start_date = datetime(year,1,1) - timedelta(days=80)
    end_date = datetime(year + 1,1,1) + timedelta(days=20)
    dt_model = (end_date - start_date).days * 24
    
    # time resolution (hr)
    t_step = 1
    
    # save to file
    folder_output = '/home/wim/Desktop/Projects/Lunar tide model/data/'
    output_name = 'ocean_m2_' + str(year) + '.dat'
    
    # include surface topography in this file (1) or not (0)
    surf_top = 0
    
    #%
    
    def tide_elevation(date_start,
                       n_steps,
                       dt,
                       radial_tide=False):
    
        # Create handler
        short_tide = pyfes.Handler("ocean",
                                   "memory",
                                   ('/home/wim/anaconda3/pkgs/pyfes-2.9.2-py38h6bb024c_0/' +
                                    'info/recipe/data/fes2014/ocean_tide_m2.ini'))
        
        load_tide = pyfes.Handler("radial",
                                   "memory",
                                   ('/home/wim/anaconda3/pkgs/pyfes-2.9.2-py38h6bb024c_0/' +
                                    'info/recipe/data/fes2014/load_tide_m2.ini'))
        
        # Create grid that will be used to interpolate the tide
        lats = np.arange(-90, 90, 1.0)
        lons = np.arange(0, 360, 2.0)
        
        full_arr = np.empty((n_steps,
                             np.shape(lats)[0],
                             np.shape(lons)[0]))
        
        lons, lats = np.meshgrid(lons, lats)
        
        shape = lons.shape
        
        date = date_start
        dates = np.empty(shape, dtype='datetime64[us]')
            
        for t in range(n_steps):
            dates.fill(date)
            
            # Create handler; output diurnal + semidiurnal, and long-period tides
            tide, lp, _ = short_tide.calculate(lons.ravel(), 
                                               lats.ravel(),
                                               dates.ravel())
            
            tide, lp = tide.reshape(shape), lp.reshape(shape)
            
            if radial_tide:
                load, load_lp, _ = load_tide.calculate(lons.ravel(), 
                                                       lats.ravel(),
                                                       dates.ravel())
                
                load, load_lp = load.reshape(shape), load_lp.reshape(shape)
                
            else:
                load = np.zeros(lons.shape)
                load_lp = load
            
            # Creating an image to see the result in meters
            geo_tide = (tide + lp + load) * 0.01 # cm to meters
            geo_tide = geo_tide.reshape(lons.shape)
            geo_tide = np.nan_to_num(geo_tide,nan=0.0) * 9.81 # potential 
            
            full_arr[t,:,:] = geo_tide
            
            date = date + timedelta(hours=dt)
            
        return full_arr
    
    #% create input file
    
    # max_m and max_n, corresponds to dim 60 x 120 era5 input files
    max_m = np.float32(45)
    max_n = np.float32(45)
    max_l = np.float32(1) 
    
    zero_top = np.float32(0) # dummy value for surface geop field
    
    tdim = (dt_model + 100) * 24 # 80 days of spin-up padding
    dt_input = 1/24
    t_zero_input = 0
    
    time_steps = np.zeros(tdim,dtype=np.float32)
    time_steps[0] = 0 # era5 time averaged fields
    
    for k in range(1,tdim):
        time_steps[k] = time_steps[k - 1] + 1 / 24
        
    time_steps = np.round(time_steps,3)
    
    # write meta-data and z,t grid
    f = FortranFile(folder_output + output_name,'w')
    
    # write grid dims as integers. *2 factor to make rec len match complex len
    rec_int = np.zeros(int(max_m * max_n * max_l * 2),dtype=np.int32)
    rec_int[0], rec_int[1], rec_int[2], rec_int[3], rec_int[4], rec_int[5], rec_int[6] = (tdim,
                                                                                          max_m,
                                                                                          max_n,
                                                                                          max_l,
                                                                                          zero_top,
                                                                                          dt_input,
                                                                                          t_zero_input)
    
    f.write_record(rec_int)
    
    # write time and altitude array as floats
    rec_fl = np.zeros(int(max_m * max_n * max_l * 2),dtype=np.float32)
    rec_fl[:100] = time_steps[:100] # 100 dummy values; calculate time_steps in Fortran
    f.write_record(rec_fl)
    
    temp_coeffs = np.zeros((int(max_n),int(max_m)),
                               dtype=np.complex64)
    
    # note: surface geop is constant throughout the year
    data = Dataset('/home/wim/Desktop/Projects/Model/ERA5/era5_psurf_geopsurf_mm_clim.nc')
    mean_f = np.mean(data['z'][:,::-1,:],axis=0)
    
    mean_f = np.delete(mean_f,np.shape(mean_f)[0] - 1,axis=0) 
    mean_f = mean_f[::1,::1]
    
    tide_field = tide_elevation(start_date,
                                dt_model,
                                1,
                                radial_tide=True)
    
    dt, dy, dx = np.shape(tide_field)
        
    for t in range(dt):
        
        # convert fields to spherical harmonic representation
        shs = sh.expand.SHExpandDHC(mean_f * surf_top + tide_field[t,:,:])
                        
        temp_coeff = shs[0,:,:] # even shs
        temp_coeffs[:,:] = temp_coeff[:int(max_n),:int(max_m)]
        
        temp_coeffs = np.swapaxes(temp_coeffs,0,1)
        
        f.write_record(temp_coeffs.reshape(int(max_m * max_n * max_l),
                                                order='F'))
        
        print('Writing tide field: ',t)
        
    f.close()
    
    del tide_field

#%%
