#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 14 11:08:11 2019

@author: wim
"""

import numpy as np
from scipy.optimize import curve_fit
from datetime import datetime, timedelta 
import json
import scipy.stats
from scipy.stats import circmean
from scipy.stats import circstd
import random
from netCDF4 import Dataset
import os.path
import warnings

#%%

def measurement_field(start_date,
                      end_date,
                      var='v'):    
    
    """    
    -- start and end date: datetime objects
    -- var:  "v" or "u"
    
    SIGN OF U: HAS BEEN CHANGED?
    
    """
    
    with open("ss_tides_config.json") as config_file:
        cfg = json.load(config_file)
            
    path = cfg["SD_data"]["path"]
    stations = cfg["SD_data"]["stations"].split()
    
    # number of amplitude measurements (t_steps)
    t_steps = int((end_date - start_date) / timedelta(hours=1))

    # array to hold hourly interval measurements
    measurements = np.zeros((t_steps,np.shape(stations)[0]))
    
    # longitudes (rad) scaled to NAVGEM; corresponding to data measurements
    longitudes = np.array([-158.5,-150.1,-123.2,-105.2,-91.9,-83.3,
                           -60.3,-26.9,-18.0,25.2 - 360])
    longitudes = (longitudes + 360) / 180 * np.pi
    
    if var == 'v':
        var = 7
    elif var == 'u':
        var = 8
        
    # fill measurement field using station data from start to end date
    for index,item in enumerate(stations):
        data = np.load(path + 'Treated_full_data_extended' + item + '.npy')
        
        not_found = True
        c = 0
                
        while not_found:
            if datetime(int(data[c,0]),int(data[c,1]),
                        int(data[c,2]),int(data[c,3])) == start_date:
                measurements[:,index] = data[c:(c + t_steps),var]
                not_found = False
            else:
                c += 1
                
    return measurements, longitudes


def updated_superdarn(start_date,
                      end_date,
                      var='v',
                      optimized_spread=False):    
    
    """
    -- start and end date: datetime objects
    -- var:  "v" or "u"
    
    Sign of U has been corrected. 
    
    The (updated) pre-processing of the datafiles used in this function can be
    traced back to the Python script in the datafile's folder.
    
    Configured 11-12-2019.
    
    """

    path = '/home/wim/Desktop/NTNU/SuperDARN data/'
    stations = ['Salmon','Kodiak','PrinceGeorge','Saskatoon',
                'RankinInlet','Kapuskasing','GooseBay','Stokkseyri',
                'Pykkvibaer','Hankasalmi']
    
    # number of amplitude measurements (t_steps)
    t_steps = int((end_date - start_date) / timedelta(hours=1))

    # array to hold hourly interval measurements
    measurements = np.zeros((t_steps,np.shape(stations)[0]))
    
    # longitudes (rad) scaled to NAVGEM; corresponding to data measurements
    longitudes = np.array([-158.5,-150.1,-123.2,-105.2,-91.9,-83.3,
                           -60.3,-26.9,-18.0,25.2 - 360])
    longitudes = (longitudes + 360) / 180 * np.pi
    
    if var == 'v':
        var = 7
    elif var == 'u':
        var = 8
        
    # fill measurement field using station data from start to end date
    for index,item in enumerate(stations):
        data = np.load(path + 'SD_wind_correct_signs_' + item + '.npy')
        
        not_found = True
        c = 0
                
        while not_found:
            if datetime(int(data[c,0]),int(data[c,1]),
                        int(data[c,2]),int(data[c,3])) == start_date:
                measurements[:,index] = data[c:(c + t_steps),var]
                not_found = False
            else:
                c += 1
                
    # sample stations such as to improve equidistant longitudinal spread
    if optimized_spread:
            for t in range(t_steps):
                # choose KSR over Kod
                if np.isfinite(measurements[t,0]) and np.isfinite(measurements[t,1]):
                    measurements[t,1] = np.nan
                # choose Kap over Rkn
                if np.isfinite(measurements[t,4]) and np.isfinite(measurements[t,5]):
                    measurements[t,4] = np.nan
                # choose Sto over Pyk
                if np.isfinite(measurements[t,7]) and np.isfinite(measurements[t,8]):
                    measurements[t,8] = np.nan
                
    return measurements, longitudes

def updated_superdarn_mc_skim(start_date,
                              end_date,
                              var='v',
                              optimized_spread=False,
                              skim_n='25'):    
    
    """
    -- start and end date: datetime objects
    -- var:  "v" or "u"
    -- mc_skim --> skimmed based on meteor echo count rate,
       set count skim rate in input file (25,50,75,100)
    
    Sign of U has been corrected. 
    
    The (updated) pre-processing of the datafiles used in this function can be
    traced back to the Python script in the datafile's folder.
    
    Configured 11-12-2019.
    
    """

    path = '/home/wim/Desktop/NTNU/SuperDARN data/'
    stations = ['Salmon','Kodiak','PrinceGeorge','Saskatoon',
                'RankinInlet','Kapuskasing','GooseBay','Stokkseyri',
                'Pykkvibaer','Hankasalmi']
    
    # number of amplitude measurements (t_steps)
    t_steps = int((end_date - start_date) / timedelta(hours=1))

    # array to hold hourly interval measurements
    measurements = np.zeros((t_steps,np.shape(stations)[0]))
    
    # longitudes (rad) scaled to NAVGEM; corresponding to data measurements
    longitudes = np.array([-158.5,-150.1,-123.2,-105.2,-91.9,-83.3,
                           -60.3,-26.9,-18.0,25.2 - 360])
    longitudes = (longitudes + 360) / 180 * np.pi
    
    if var == 'v':
        var = 7
    elif var == 'u':
        var = 8
    else:
        var = var # standard deviation merid (14), zonal (15)
        
    # fill measurement field using station data from start to end date
    for index,item in enumerate(stations):
        file_mc = 'SD_wind_correct_signs_mc_' + skim_n
        data = np.load(path + file_mc + item + '.npy')
                        
        # all good in terms of STD and wind nans here. No wind without a STD
        not_found = True
        c = 0
                
        while not_found:
            if datetime(int(data[c,0]),int(data[c,1]),
                        int(data[c,2]),int(data[c,3])) == start_date:
                measurements[:,index] = data[c:(c + t_steps),var]
                not_found = False
            else:
                c += 1
                
    # sample stations such as to improve equidistant longitudinal spread   
    if var == 7 or var == 8:
        if optimized_spread:
                for t in range(t_steps):
                    # choose KSR over Kod
                    if np.isfinite(measurements[t,0]) and np.isfinite(measurements[t,1]):
                        measurements[t,1] = np.nan
                    # choose Kap over Rkn
                    if np.isfinite(measurements[t,4]) and np.isfinite(measurements[t,5]):
                        measurements[t,4] = np.nan
                    # choose Sto over Pyk
                    if np.isfinite(measurements[t,7]) and np.isfinite(measurements[t,8]):
                        measurements[t,8] = np.nan
                
    return measurements, longitudes


def new_superdarn(start_day,
                  end_day,
                  invar,
                  optimized_spread=False,
                  n_skim='25'):
    
    """
    
    invar = 7 = U
    invar = 8 = V
    
    """
    
    stations = ['.ksr','.kod','.pgr','.sas','.rkn','.kap','.gbr','.sto','.pyk','.han']
    longitudes = np.array([-158.5,-150.1,-123.2,-105.2,-91.9,-83.3,
                           -60.3,-26.9,-18.0,25.2 - 360])
    longitudes = (longitudes + 360) / 180 * np.pi
    
    min_values = 5
    
    if invar == 7:
        covar = 14
    elif invar == 8:
        covar = 15
    
    n_stations = np.shape(stations)[0]
        
    n_days = (end_day - start_day).days
    length = n_days * 24
    
    wind = np.zeros((length,n_stations)) * np.nan  
    std = np.zeros((length,n_stations)) * np.nan
    n_meteor = np.zeros((length,n_stations)) * np.nan
        
    f_base = '/home/wim/Desktop/Projects/SSW/data/SuperDARN/data_files/'
    
    loop_date = start_day
    
    for d in range(n_days):
        # set path to folder based on loop_date
        yrs = loop_date.strftime("%Y")
        mon = loop_date.strftime("%b") 
        day = loop_date.strftime("%d")
        mnt = loop_date.strftime("%m")
        path = f_base + yrs + '/' + mnt + '/' + yrs + mon + day
                    
        # read in data for each station
        for index, item in enumerate(stations):
            if os.path.exists(path + item + '.txt'):
                # read in data file
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    data = np.loadtxt(path +  item + '.txt')
                
                # calculate time index based on start date
                t_index = (loop_date - start_day).days
                t_ind = t_index * 24
                
                # if all data is present, store data in array
                if np.shape(data) == (24,16):
                    temp = data[:,invar]   
                    if np.nanmean(np.abs(temp)) > min_values:
                        wind[(t_ind):(t_ind + 24),index] = data[:,invar]                    
                        std[(t_ind):(t_ind + 24),index] = data[:,covar]
                        n_meteor[(t_ind):(t_ind + 24),index] = data[:,4]
                        
                    else:
                        None
    
                else:
                    None
                    
        loop_date = loop_date + timedelta(days=1)
            
    inds_std = np.where(std == 0)
    wind[inds_std] = np.nan
    
    inds_mag = np.where(wind > 150)
    wind[inds_mag] = np.nan
    
    inds_met = np.where(n_meteor < float(n_skim))
    
    wind[inds_met] = np.nan
    
    inds_los = np.where(np.isnan(wind))
   
    # sample stations such as to improve equidistant longitudinal spread   
    
    t_steps = (end_day - start_day).days * 24
    if optimized_spread:
            for t in range(t_steps):
                # choose KSR over Kod
                if np.isfinite(wind[t,0]) and np.isfinite(wind[t,1]):
                    wind[t,1] = np.nan
                # choose Kap over Rkn
                if np.isfinite(wind[t,4]) and np.isfinite(wind[t,5]):
                    wind[t,4] = np.nan
                # choose Sto over Pyk
                if np.isfinite(wind[t,7]) and np.isfinite(wind[t,8]):
                    wind[t,8] = np.nan
    
    return wind, std, longitudes

def optimized_measurement_field(start_date,
                                end_date,
                                var="v"):    
    
    """
    -- start and end date: datetime objects
    -- var:  "v" or "u"
    
    """
    
    with open("ss_tides_config.json") as config_file:
        cfg = json.load(config_file)
            
    path = cfg["SD_data"]["path"]
    stations = cfg["SD_data"]["stations"].split()
    
    # number of amplitude measurements (t_steps)
    t_steps = int((end_date - start_date) / timedelta(hours=1))

    # array to hold hourly interval measurements
    measurements = np.zeros((t_steps,np.shape(stations)[0]))
    
    # longitudes (rad) scaled to NAVGEM; corresponding to data measurements
    longitudes = np.array([-158.5,-150.1,-123.2,-105.2,-91.9,-83.3,
                           -60.3,-26.9,-18.0,25.2 - 360])   
    longitudes = (longitudes + 360) / 180 * np.pi 
    
    if var == "v":
        var = 7
    elif var == "u":
        var = 8
        
    # fill measurement field using station data from start to end date
    for index,item in enumerate(stations):
        data = np.load(path + 'Treated_full_data_extended' + item + '.npy')
        
        not_found = True
        c = 0
                
        while not_found:
            if datetime(int(data[c,0]),int(data[c,1]),
                        int(data[c,2]),int(data[c,3])) == start_date:
                measurements[:,index] = data[c:(c + t_steps),var]
                not_found = False
            else:
                c += 1

    # 1. Optimize equidistant spread
    for t in range(t_steps):
        if np.isfinite(measurements[t,0]) and np.isfinite(measurements[t,1]):
            measurements[t,1] = np.nan
        if np.isfinite(measurements[t,7]) and np.isfinite(measurements[t,8]):
            measurements[t,7] = np.nan
        if np.isfinite(measurements[t,4]) and np.isfinite(measurements[t,5]):
            measurements[t,4] = np.nan
            
    return measurements, longitudes


def optimized_new_SD(start_date,
                     end_date,
                     var="v"):    
    
    """
    -- start and end date: datetime objects
    -- var:  "v" or "u"
    
    """
    
    with open("ss_tides_config.json") as config_file:
        cfg = json.load(config_file)
            
    path = '/home/wim/Desktop/NTNU/SuperDARN data/'
    stations = cfg["SD_data"]["stations"].split()
    
    # number of amplitude measurements (t_steps)
    t_steps = int((end_date - start_date) / timedelta(hours=1))

    # array to hold hourly interval measurements
    measurements = np.zeros((t_steps,np.shape(stations)[0]))
    
    # longitudes (rad) scaled to NAVGEM; corresponding to data measurements
    longitudes = np.array([-158.5,-150.1,-123.2,-105.2,-91.9,-83.3,
                           -60.3,-26.9,-18.0,25.2 - 360])   
    longitudes = (longitudes + 360) / 180 * np.pi 

    if var == "v":
        var = 7
        print('SuperDARN variable ' + var + ' is extracted')
    if var == "u":
        var = 8
        print('SuperDARN variable ' + var + ' is extraxted')
        
    # fill measurement field using station data from start to end date
    for index,item in enumerate(stations):
        data = np.load(path + 'Treated_full_data_extended' + item + '.npy')
        
        not_found = True
        c = 0
                
        while not_found:
            if datetime(int(data[c,0]),int(data[c,1]),
                        int(data[c,2]),int(data[c,3])) == start_date:
                measurements[:,index] = data[c:(c + t_steps),var]
                not_found = False
            else:
                c += 1

    # 1. Optimize equidistant spread
    for t in range(t_steps):
        if np.isfinite(measurements[t,0]) and np.isfinite(measurements[t,1]):
            measurements[t,1] = np.nan
        if np.isfinite(measurements[t,7]) and np.isfinite(measurements[t,8]):
            measurements[t,7] = np.nan
        if np.isfinite(measurements[t,4]) and np.isfinite(measurements[t,5]):
            measurements[t,4] = np.nan
            
    return measurements, longitudes


def interpolated_field_fixlat(height_i,
                              lat_i,
                              nan_map,
                              var="v"):
    
    """
    -- top_i = 00 at 95km; 1.25km increments
    -- lat_i = 09 at 60N
    
    -- function checked 23-09-2019
        
    """
    
    with open("ss_tides_config.json") as config_file:
        cfg = json.load(config_file)
            
    folder = cfg["NAVGEM_data"]["path_" + var]
    months = cfg["NAVGEM_data"]["months"].split()    
        
    for index, month in enumerate(months):                  
        data_file = np.load(folder + month + '.npy')[:,height_i,lat_i,:]

        # remove rare interpolation error giving np.nan entry
        inds = np.isnan(data_file)
        data_file[inds] = 0
            
        if index == 0:
            field = data_file
        else:
            field = np.concatenate((field,data_file),axis=0)
            
    # filter data at nan locations
    length = np.isnan(nan_map).sum()
    nans = np.isnan(nan_map)
    
    filtered_field = np.zeros((np.shape(field)[0],length))
    longitude_field = np.linspace(0,2 * np.pi, 360, endpoint=False)[nans]
    
    for t in range(np.shape(field)[0]):
        filtered_field[t,:] = field[t,:][nans] 
        
    return filtered_field, longitude_field


def interpolated_field_SDlat(height_i,
                             var='v',
                             year='2014'):
    
    """
    Function loads interpolated NAVGEM-HA fields sampled at SuperDARN locs
    
    top_i = 00 at 95km, -1.25 km intervals down to 80 km
    
    fuction checked: 23-05-2019
    function revised to include 2015:  

    """
    
    with open("ss_tides_config.json") as config_file:
        cfg = json.load(config_file)
            
    # import meta data
    folder = cfg["NAVGEM_data"]["path_" + var]
    months = cfg["NAVGEM_data"]["months"].split()    
    nav_lats = cfg["NAVGEM_data"]["lat"].split()   
    sd_lats = cfg["SD_data"]["lat"].split()
    
    # calculate NAVGEM latitude indices corresponding to SD locations
    nav_lats = np.array(nav_lats,dtype='float')
    sd_lats = np.array(sd_lats,dtype='float')
    
    lat_inds = np.empty(10)
        
    for index, item in enumerate(sd_lats):
        idy = (np.abs(nav_lats - item)).argmin()
        lat_inds[index] = idy
          
    # correct for NAVGEM longitude ordering
    lat_inds = np.roll(lat_inds,1)
    
    # calculate NAVGEM longitude nan-map corresponding to SD locations
    stations_lon = np.array([-158.5,-150.1,-123.2,-105.2,-91.9,-83.3,
                             -60.3,-26.9,-18.0,25.2 - 360]) 

    lons = (360 + stations_lon) / 180 * np.pi
    navgem_lons = np.linspace(0,2 * np.pi, 360, endpoint=False)
    nan_map = np.zeros(360)
    
    for index, item in enumerate(lons):
        ind = np.argmin(np.abs(navgem_lons - item))
        nan_map[ind] = np.nan
    
    # populate NAVGEM data field
    for index, month in enumerate(months):                  
        data_file = np.load(folder + month + '.npy')[:,height_i,:,:]

        # set extremely rare interpolation error to zero
        inds = np.isnan(data_file)
        data_file[inds] = 0
        
        if index == 0:
            field = data_file
        else:
            field = np.concatenate((field,data_file),axis=0)
            
    # filter data at nan locations
    length = np.isnan(nan_map).sum()
    nans = np.isnan(nan_map)
    nan_inds = np.where(nans)[0]
    
    filtered_field = np.zeros((np.shape(field)[0],length))
    longitude_field = np.linspace(0,2 * np.pi, 360, endpoint=False)[nans]
    
    for t in range(np.shape(field)[0]): 
        for i in range(length):
            filtered_field[t,i] = field[t,int(lat_inds[i]),int(nan_inds[i])]
       
    return filtered_field, longitude_field


def NAVGEM_field_SD_sampled(height_i,
                            var='v',
                            year='2014'):
    
    """
    Function loads interpolated NAVGEM-HA fields sampled at SuperDARN locs.
    Ordering of the fields is equal to that of SuperDARN data, i.e. index 0 
    corresponds to KSR, index 9 to HAN.
    
    top_i = 00 at 95km, -1.25 km intervals down to 80 km.
    
    function revised to include 2015: 10-12-2019.
    
    Works with NAVGEM-HA interpolated to an index [141:157] latitude band.

    """
    
    # import meta data
    folder = '/home/wim/Desktop/NTNU/NAVGEM_HA_Data/' + year + '_interpolated/' + year + '_interpolated_' + var + '_month_'
    months = '01 02 03 04 05 06 07 08 09 10 11 12'.split()
    nav_lats = '51.35 52.35 53.35 54.34 55.34 56.34 57.34 58.33 59.33 60.3 61.32 62.32 63.32 64.32 65.31 66.31'.split()   
    sd_lats = '60.6 59.5 56.1 54.2 65.0 51.4 55.6 64.7 65.7 64.4'.split()
    
    # print meta data of sampled fields
    print('NAVGEM-HA sampled at superDARN locations. \n')
    print('Height: ',np.arange(95000,79999,-1250, dtype='float')[height_i], ' [m]')
    print('Variable: ',var,' [m/s]')
    print('Year: ',year)
    
    # calculate NAVGEM latitude indices corresponding to SD locations
    nav_lats = np.array(nav_lats,dtype='float')
    sd_lats = np.array(sd_lats,dtype='float')
    lat_inds = np.empty(10) 

    for index, item in enumerate(sd_lats):
        idy = np.abs(nav_lats - item).argmin()
        lat_inds[index] = idy

    # calculate NAVGEM longitude nan-map corresponding to SD locations
    stations_lon = np.array([-158.5,-150.1,-123.2,-105.2,-91.9,-83.3,
                             -60.3,-26.9,-18.0,25.2 - 360]) 

    sd_lons = np.array((360 + stations_lon) / 180 * np.pi)
    nav_lons = np.linspace(0,2 * np.pi,360, endpoint=False)
    lon_inds = np.empty(10)
        
    for index, item in enumerate(sd_lons):
        idx = np.abs(nav_lons - item).argmin()
        lon_inds[index] = idx
        
    # NAVGEM-HA field from interpolated numpy arrays
    for index, month in enumerate(months):                  
        data_file = np.load(folder + month + '.npy')[:,height_i,:,:]

        # set extremely rare interpolation error to zero
        inds = np.isnan(data_file)
        data_file[inds] = 0
        
        if index == 0:
            field = data_file
        else:
            field = np.concatenate((field,data_file),axis=0)
            
    # populate SuperDARN-sampled field
    filtered_field = np.zeros((np.shape(field)[0],10))
    longitude_field = np.empty(10)

    for t in range(np.shape(field)[0]): 
        for i in range(10):
            filtered_field[t,i] = field[t,int(lat_inds[i]),int(lon_inds[i])]
            longitude_field[i] = nav_lons[int(lon_inds[i])]
            
    return filtered_field, longitude_field


def NAVGEM_field_SD_lon_sampled(height_i,
                                var='v',
                                year='2014'):
    
    """
    Function loads interpolated NAVGEM-HA fields sampled at SuperDARN locs only
    in LONGITUDE and along 60 degrees LATITUDE.
    
    Ordering of the fields is equal to that of SuperDARN data, i.e. index 0 
    corresponds to KSR, index 9 to HAN.
    
    top_i = 00 at 95km, -1.25 km intervals down to 80 km.
    
    function revised to include 2015: 10-12-2019.
    
    Works with NAVGEM-HA interpolated to an index [141:157] latitude band.

    """
    
    # import meta data
    folder = '/home/wim/Desktop/NTNU/NAVGEM_HA_Data/' + year + '_interpolated/' + year + '_interpolated_' + var + '_month_'
    months = '01 02 03 04 05 06 07 08 09 10 11 12'.split()
    nav_lats = '51.35 52.35 53.35 54.34 55.34 56.34 57.34 58.33 59.33 60.3 61.32 62.32 63.32 64.32 65.31 66.31'.split()   
    sd_lats = '60.0 60.0 60.0 60.0 60.0 60.0 60.0 60.0 60.0 60.0'.split()
    
    # print meta data of sampled fields
    print('NAVGEM-HA sampled at superDARN longitude locations. \n')
    print('Height: ',np.arange(95000,79999,-1250, dtype='float')[height_i], ' [m]')
    print('Variable: ',var,' [m/s]')
    print('Year: ',year)
    
    # calculate NAVGEM latitude indices corresponding to SD locations
    nav_lats = np.array(nav_lats,dtype='float')
    sd_lats = np.array(sd_lats,dtype='float')
    lat_inds = np.empty(10) 

    for index, item in enumerate(sd_lats):
        idy = np.abs(nav_lats - item).argmin()
        lat_inds[index] = idy

    # calculate NAVGEM longitude nan-map corresponding to SD locations
    stations_lon = np.array([-158.5,-150.1,-123.2,-105.2,-91.9,-83.3,
                             -60.3,-26.9,-18.0,25.2 - 360]) 

    sd_lons = np.array((360 + stations_lon) / 180 * np.pi)
    nav_lons = np.linspace(0,2 * np.pi,360, endpoint=False)
    lon_inds = np.empty(10)
        
    for index, item in enumerate(sd_lons):
        idx = np.abs(nav_lons - item).argmin()
        lon_inds[index] = idx
        
    # NAVGEM-HA field from interpolated numpy arrays
    for index, month in enumerate(months):                  
        data_file = np.load(folder + month + '.npy')[:,height_i,:,:]

        # set extremely rare interpolation error to zero
        inds = np.isnan(data_file)
        data_file[inds] = 0
        
        if index == 0:
            field = data_file
        else:
            field = np.concatenate((field,data_file),axis=0)
            
    # populate SuperDARN-sampled field
    filtered_field = np.zeros((np.shape(field)[0],10))
    longitude_field = np.empty(10)

    for t in range(np.shape(field)[0]): 
        for i in range(10):
            filtered_field[t,i] = field[t,int(lat_inds[i]),int(lon_inds[i])]
            longitude_field[i] = nav_lons[int(lon_inds[i])]
            
    return filtered_field, longitude_field


def NAVGEM_field_SD_full_sampled(height_i,
                                 var='v',
                                 year='2014'):
    
    """
    Function loads interpolated NAVGEM-HA fields sampled along longitude circle
    at 60N. Index 9 in NAVGEM-HA files corresponds to 60 degrees North.
        
    top_i = 00 at 95km, -1.25 km intervals down to 80 km.
    
    function revised to include 2015: 10-12-2019.
    
    Works with NAVGEM-HA interpolated to an index [141:157] latitude band.

    """
    
    # import meta data
    folder = '/home/wim/Desktop/NTNU/NAVGEM_HA_Data/' + year + '_interpolated/' + year + '_interpolated_' + var + '_month_'
    months = '01 02 03 04 05 06 07 08 09 10 11 12'.split()
    
    # print meta data of sampled fields
    print('NAVGEM-HA along 60N. \n')
    print('Height: ',np.arange(95000,79999,-1250, dtype='float')[height_i], ' [m]')
    print('Variable: ',var,' [m/s]')
    print('Year: ',year)
       
    # NAVGEM-HA field from interpolated numpy arrays
    for index, month in enumerate(months):                  
        data_file = np.load(folder + month + '.npy')[:,height_i,:,:]

        # set extremely rare interpolation error to zero
        inds = np.isnan(data_file)
        data_file[inds] = 0
        
        if index == 0:
            field = data_file
        else:
            field = np.concatenate((field,data_file),axis=0)
            
    # populate SuperDARN-sampled field
    filtered_field = field[:,9,:] # all longitudes and times at 60N
    longitude_field = np.linspace(0,2 * np.pi,360,endpoint=False)

    return filtered_field, longitude_field


def NAVGEM_field_SD_full_sampled_varlat(height_i,
                                        var='v',
                                        year='2014',
                                        latitude=9,
                                        lon_step=1):
    
    """
    Function loads interpolated NAVGEM-HA fields sampled along longitude circle
    at 60N. Index 9 in NAVGEM-HA files corresponds to 60 degrees North. 
        
    top_i = 00 at 95km, -1.25 km intervals down to 80 km.
    
    function revised to include 2015: 10-12-2019.
    
    Works with NAVGEM-HA interpolated to an index [141:157] latitude band,
    16 latitude indices total.

    """
    
    # import meta data
    folder = '/home/wim/Desktop/NTNU/NAVGEM_HA_Data/' + year + '_interpolated/' + year + '_interpolated_extended_' + var + '_month_'
    months = '01 02 03 04 05 06 07 08 09 10 11 12'.split()
    
    # print meta data of sampled fields
    print('NAVGEM-HA along 60N. \n')
    print('Height: ',np.arange(95000,79999,-1250, dtype='float')[height_i], ' [m]')
    print('Variable: ',var,' [m/s]')
    print('Year: ',year)
    print('Latitude index: ', latitude)  
    print('Longitude step: ', lon_step)     
    
    # NAVGEM-HA field from interpolated numpy arrays
    for index, month in enumerate(months):                  
        data_file = np.load(folder + month + '.npy')[:,height_i,:,:]

        # set extremely rare interpolation error to zero
        inds = np.isnan(data_file)
        data_file[inds] = 0
        
        if index == 0:
            field = data_file
        else:
            field = np.concatenate((field,data_file),axis=0)
            
    # populate SuperDARN-sampled field
    filtered_field = field[:,latitude,::lon_step] 
    longitude_field = np.linspace(0,2 * np.pi,
                                  np.shape(filtered_field)[1],
                                  endpoint=False)

    return filtered_field, longitude_field


def superdarn_lon_nan_map():
    
    """
    returns longitudinal nan-map of SD locs in NAVGEM (checked)
    
    """
    
    stations_lon = np.array([-158.5,-150.1,-123.2,-105.2,-91.9,-83.3,
                        -60.3,-26.9,-18.0,25.2 - 360]) 

    lons = (360 + stations_lon) / 180 * np.pi
    navgem_lons = np.linspace(0,2 * np.pi, 360, endpoint=False)
    nan_map = np.zeros(360)
    
    for index, item in enumerate(lons):
        ind = np.argmin(np.abs(navgem_lons - item))
        nan_map[ind] = np.nan
        
    return nan_map


def NAVGEM_full_nan_map():
    
    """
    all 360 longitude points to sample
    
    """
    
    nan_map = np.empty(360) * np.nan
        
    return nan_map


def NAVGEM_half_nan_map():
    
    """
    all 360 longitude points to sample
    
    """
    
    nan_map = np.zeros(360) * np.nan
    nan_map[25:201] = 0
        
    return nan_map


def dt(x,
       A,
       phi):
    
    """
    x: [0,2 * np.pi] radians
    A: Amplitude (m/s)
    phi: phase (rad)
    
    """
    
    return A * np.sin(1 * x + phi)


def sdt(x,
        A,
        phi):
    
    """
    x: [0,2 * np.pi] radians
    A: Amplitude (m/s)
    phi: phase (rad)
    
    """
    
    return A * np.sin(2 * x + phi)


def tdt(x,
        A,
        phi):
    
    """
    x: [0,2 * np.pi] radians
    A: Amplitude (m/s)
    phi: phase (rad)
    
    """
    
    return A * np.sin(3 * x + phi)


def dt_sdt_tdt(x,
              A,
              B,
              C,
              phi_A,
              phi_B,
              phi_C):
    
    """
    x: [0,2 * np.pi] radians
    
    """
    
    return (A * np.sin(1 * x + phi_A) + 
            B * np.sin(2 * x + phi_B) +
            C * np.sin(3 * x + phi_C))


def qdt(x,
        A,
        phi):
    
    """
    x: [0,2 * np.pi] radians
    A: Amplitude (m/s)
    phi: phase (rad)
    
    """
    
    return A * np.sin(4 * x + phi)


def q2dw_50hr(x,
              A,
              phi):
    
    """
    x: [0,2 * np.pi] radians
    A: Amplitude (m/s)
    phi: phase (rad)
    
    """
    
    return A * np.sin((24 / 50) * x + phi)


def navgem_dt_sdt_tdt(x,
                      A,
                      B,
                      C,
                      phi_A,
                      phi_B,
                      phi_C):
    
    """
    function checked 23-05-2019
    
    """
    
    return (  A * np.sin(1 * ( x + (x // (4 * np.pi) % 8) * 2 * np.pi / 8 ) + phi_A)
            + B * np.sin(2 * ( x + (x // (4 * np.pi) % 8) * 2 * np.pi / 8 ) + phi_B)
            + C * np.sin(3 * ( x + (x // (4 * np.pi) % 8) * 2 * np.pi / 8 ) + phi_C))


def navgem_dt_sdt_tdt_mean(x,
                           A,
                           B,
                           C,
                           phi_A,
                           phi_B,
                           phi_C,
                           mean):
    
    """
    function checked 23-05-2019
    
    """
    
    return (  A * np.sin(1 * ( x + (x // (4 * np.pi) % 8) * 2 * np.pi / 8 ) + phi_A)
            + B * np.sin(2 * ( x + (x // (4 * np.pi) % 8) * 2 * np.pi / 8 ) + phi_B)
            + C * np.sin(3 * ( x + (x // (4 * np.pi) % 8) * 2 * np.pi / 8 ) + phi_C)
            + mean)

def navgem_dt_sdt_tdt_mean_sw1_sw3(x,
                                   A,
                                   B,
                                   C,
                                   phi_A,
                                   phi_B,
                                   phi_C,
                                   D,
                                   E,
                                   phi_D,
                                   phi_E,
                                   mean):
    
    """
    function checked 23-05-2019
    
    """
    
    return (  A * np.sin(1 * ( x + (x // (4 * np.pi) % 8) * 2 * np.pi / 8 ) + phi_A)
            + B * np.sin(2 * ( x + (x // (4 * np.pi) % 8) * 2 * np.pi / 8 ) + phi_B)
            + C * np.sin(3 * ( x + (x // (4 * np.pi) % 8) * 2 * np.pi / 8 ) + phi_C)
            + D * np.sin(1 * ( x + (x // (4 * np.pi) % 8) * 2 * np.pi / 4 ) + phi_D)
            + E * np.sin(3 * ( x + (x // (4 * np.pi) % 8) * 2 * np.pi / 12 ) + phi_E)
            + mean)


def navgem_dt_sdt_tdt_mean_sw1_se2(x,
                                   A,
                                   B,
                                   C,
                                   phi_A,
                                   phi_B,
                                   phi_C,
                                   D,
                                   E,
                                   phi_D,
                                   phi_E,
                                   mean):
    
    """
    function checked 23-05-2019
    
    """
    
    return (  A * np.sin(1 * ( x + (x // (4 * np.pi) % 8) * 2 * np.pi / 8 ) + phi_A)
            + B * np.sin(2 * ( x + (x // (4 * np.pi) % 8) * 2 * np.pi / 8 ) + phi_B)
            + C * np.sin(3 * ( x + (x // (4 * np.pi) % 8) * 2 * np.pi / 8 ) + phi_C)
            + D * np.sin(1 * ( x + (x // (4 * np.pi) % 8) * 2 * np.pi / 4 ) + phi_D)
            + E * np.sin(2 * ( x - (x // (4 * np.pi) % 8) * 2 * np.pi / 8 ) + phi_E)
            + mean)


def navgem_dt_sdt_tdt_mean_sw1_se1(x,
                                   A,
                                   B,
                                   C,
                                   phi_A,
                                   phi_B,
                                   phi_C,
                                   D,
                                   E,
                                   phi_D,
                                   phi_E,
                                   mean):
    
    """
    function checked 23-05-2019
    
    """
    
    return (  A * np.sin(1 * ( x + (x // (4 * np.pi) % 8) * 2 * np.pi / 8 ) + phi_A)
            + B * np.sin(2 * ( x + (x // (4 * np.pi) % 8) * 2 * np.pi / 8 ) + phi_B)
            + C * np.sin(3 * ( x + (x // (4 * np.pi) % 8) * 2 * np.pi / 8 ) + phi_C)
            + D * np.sin(1 * ( x + (x // (4 * np.pi) % 8) * 2 * np.pi / 4 ) + phi_D)
            + E * np.sin(1 * ( x - (x // (4 * np.pi) % 8) * 2 * np.pi / 4 ) + phi_E)
            + mean)


def navgem_dt_sdt_tdt_mean_sw1(x,
                               A,
                               B,
                               C,
                               phi_A,
                               phi_B,
                               phi_C,
                               D,
                               phi_D,
                               mean):
    
    """
    function checked 23-05-2019
    
    """
    
    return (  A * np.sin(1 * ( x + (x // (4 * np.pi) % 8) * 2 * np.pi / 8 ) + phi_A)
            + B * np.sin(2 * ( x + (x // (4 * np.pi) % 8) * 2 * np.pi / 8 ) + phi_B)
            + C * np.sin(3 * ( x + (x // (4 * np.pi) % 8) * 2 * np.pi / 8 ) + phi_C)
            + D * np.sin(1 * ( x + (x // (4 * np.pi) % 8) * 2 * np.pi / 4 ) + phi_D)
            + mean)
    


def navgem_nonmigrating(x,
                        A,
                        B,
                        C,
                        phi_A,
                        phi_B,
                        phi_C,
                        mean,
                        F,
                        phi_F):
    
    """
    function checked 23-05-2019
    
    """
    
    return (  A * np.sin(1 * ( x + (x // (4 * np.pi) % 8) * 2 * np.pi / 8 ) + phi_A)
            + B * np.sin(2 * ( x + (x // (4 * np.pi) % 8) * 2 * np.pi / 8 ) + phi_B)
            + C * np.sin(3 * ( x + (x // (4 * np.pi) % 8) * 2 * np.pi / 8 ) + phi_C)
            + mean
            + F * np.sin(1 * x + 2 * (x // (4 * np.pi) % 8) * 2 * np.pi / 8  + phi_F))


def superdarn_dt_sdt_tdt(x,
                         A,
                         B,
                         C,
                         phi_A,
                         phi_B,
                         phi_C):
    
    return (  A * np.sin(1 * ( x + (x // (4 * np.pi) % 24) * 2 * np.pi / 24 ) + phi_A)
            + B * np.sin(2 * ( x + (x // (4 * np.pi) % 24) * 2 * np.pi / 24 ) + phi_B)
            + C * np.sin(3 * ( x + (x // (4 * np.pi) % 24) * 2 * np.pi / 24 ) + phi_C))

def superdarn_dt_sdt_tdt_mean(x,
                              A,
                              B,
                              C,
                              phi_A,
                              phi_B,
                              phi_C,
                              M):
    
    return (  A * np.sin(1 * ( x + (x // (4 * np.pi) % 24) * 2 * np.pi / 24 ) + phi_A)
            + B * np.sin(2 * ( x + (x // (4 * np.pi) % 24) * 2 * np.pi / 24 ) + phi_B)
            + C * np.sin(3 * ( x + (x // (4 * np.pi) % 24) * 2 * np.pi / 24 ) + phi_C)
            + M)

def superdarn_dt_m2_tdt_mean(x,
                             A,
                             B,
                             C,
                             phi_A,
                             phi_B,
                             phi_C,
                             M):
    
    return (  A * np.sin(1 * ( x + (x // (4 * np.pi) % 24) * 2 * np.pi / 24 ) + phi_A)
            + B * np.sin(2 * ( x + (x // (4 * np.pi) * 2 * np.pi / 24.84 )) + phi_B)
            + C * np.sin(3 * ( x + (x // (4 * np.pi) % 24) * 2 * np.pi / 24 ) + phi_C)
            + M)

def superdarn_dt_sdt_tdt_mean_sw1(x,
                                  A,
                                  B,
                                  C,
                                  phi_A,
                                  phi_B,
                                  phi_C,
                                  D,
                                  phi_D,
                                  M):
    
    return (  A * np.sin(1 * ( x + (x // (4 * np.pi) % 24) * 2 * np.pi / 24 ) + phi_A)
            + B * np.sin(2 * ( x + (x // (4 * np.pi) % 24) * 2 * np.pi / 24 ) + phi_B)
            + C * np.sin(3 * ( x + (x // (4 * np.pi) % 24) * 2 * np.pi / 24 ) + phi_C)
            + D * np.sin(1 * ( x + (x // (4 * np.pi) % 24) * 2 * np.pi / 12 ) + phi_D)
            + M)

def superdarn_dt_sdt_tdt_mean_sw1_sw3(x,
                                      A,
                                      B,
                                      C,
                                      phi_A,
                                      phi_B,
                                      phi_C,
                                      D,
                                      E,
                                      phi_D,
                                      phi_E,
                                      M):
    
    return (  A * np.sin(1 * ( x + (x // (4 * np.pi) % 24) * 2 * np.pi / 24 ) + phi_A)
            + B * np.sin(2 * ( x + (x // (4 * np.pi) % 24) * 2 * np.pi / 24 ) + phi_B)
            + C * np.sin(3 * ( x + (x // (4 * np.pi) % 24) * 2 * np.pi / 24 ) + phi_C)
            + D * np.sin(1 * ( x + (x // (4 * np.pi) % 24) * 2 * np.pi / 12 ) + phi_D)
            + E * np.sin(3 * ( x + (x // (4 * np.pi) % 24) * 2 * np.pi / 36 ) + phi_E)
            + M)

def superdarn_dt_sdt_tdt_mean_sw1_se2(x,
                                      A,
                                      B,
                                      C,
                                      phi_A,
                                      phi_B,
                                      phi_C,
                                      D,
                                      E,
                                      phi_D,
                                      phi_E,
                                      M):
    
    return (  A * np.sin(1 * ( x + (x // (4 * np.pi) % 24) * 2 * np.pi / 24 ) + phi_A)
            + B * np.sin(2 * ( x + (x // (4 * np.pi) % 24) * 2 * np.pi / 24 ) + phi_B)
            + C * np.sin(3 * ( x + (x // (4 * np.pi) % 24) * 2 * np.pi / 24 ) + phi_C)
            + D * np.sin(1 * ( x + (x // (4 * np.pi) % 24) * 2 * np.pi / 12 ) + phi_D)
            + E * np.sin(2 * ( x - (x // (4 * np.pi) % 24) * 2 * np.pi / 24 ) + phi_E)
            + M)

def superdarn_dt_sdt_tdt_mean_sw1_se1(x,
                                      A,
                                      B,
                                      C,
                                      phi_A,
                                      phi_B,
                                      phi_C,
                                      D,
                                      E,
                                      phi_D,
                                      phi_E,
                                      M):
    
    return (  A * np.sin(1 * ( x + (x // (4 * np.pi) % 24) * 2 * np.pi / 24 ) + phi_A)
            + B * np.sin(2 * ( x + (x // (4 * np.pi) % 24) * 2 * np.pi / 24 ) + phi_B)
            + C * np.sin(3 * ( x + (x // (4 * np.pi) % 24) * 2 * np.pi / 24 ) + phi_C)
            + D * np.sin(1 * ( x + (x // (4 * np.pi) % 24) * 2 * np.pi / 12 ) + phi_D)
            + E * np.sin(1 * ( x - (x // (4 * np.pi) % 24) * 2 * np.pi / 12 ) + phi_E)
            + M)

    
def superdarn_dt_sdt_tdt_flip(x,
                              A,
                              B,
                              C,
                              phi_A,
                              phi_B,
                              phi_C):
    
    return (  A * np.sin(1 * ( x + (x // (4 * np.pi) % 24) * 2 * np.pi / 24 ) + phi_A)
            + B * np.sin(2 * ( x - (x // (4 * np.pi) % 24) * 2 * np.pi / 24 ) + phi_B)
            + C * np.sin(3 * ( x + (x // (4 * np.pi) % 24) * 2 * np.pi / 24 ) + phi_C))
    

def superdarn_SL_SDT(x,
                     A,
                     B,
                     phi_A,
                     phi_B):
    
    return (  A * np.sin(2 * ( x + x // (4 * np.pi) * 2 * np.pi * 50 / 1200 ) + phi_A)
            + B * np.sin(2 * ( x + x // (4 * np.pi) * 2 * np.pi * 50 / 1242 ) + phi_B)) 
    
    
def navgem_SL_SDT(x,
                  A,
                  B,
                  phi_A,
                  phi_B):
    
    return (  A * np.sin(2 * ( x + x // (4 * np.pi) * 2 * np.pi * 50 / 1200 * 3) + phi_A)
            + B * np.sin(2 * ( x + x // (4 * np.pi) * 2 * np.pi * 50 / 1242 * 3) + phi_B)) 
    
    
def navgem_SL_m2n2(x,
                   A,
                   B,
                   phi_A,
                   phi_B,
                   D,
                   phi_D):
    
    return (  A * np.sin(2 * ( x + x // (4 * np.pi) * 2 * np.pi * 50 / 1200 * 3) + phi_A)
            + B * np.sin(2 * ( x + x // (4 * np.pi) * 2 * np.pi * 50 / 1242 * 3) + phi_B)
            + D * np.sin(2 * ( x + x // (4 * np.pi) * 2 * np.pi * 50 / 1266 * 3) + phi_D)) 

def superdarn_SL_SDT_m2n2(x,
                          A,
                          B,
                          phi_A,
                          phi_B,
                          C,
                          phi_C):
    
    return (  A * np.sin(2 * ( x + x // (4 * np.pi) * 2 * np.pi * 50 / 1200 ) + phi_A)
            + B * np.sin(2 * ( x + x // (4 * np.pi) * 2 * np.pi * 50 / 1242 ) + phi_B)
            + C * np.sin(2 * ( x + x // (4 * np.pi) * 2 * np.pi * 50 / 1266 ) + phi_C))
    
    
def superdarn_SL_SDT_m2n2_eastwest(x,
                                   A,
                                   B,
                                   phi_A,
                                   phi_B,
                                   C,
                                   D,
                                   E,
                                   F,
                                   phi_C,
                                   phi_D,
                                   phi_E,
                                   phi_F):
    
    return (  A * np.sin(2 * ( x + x // (4 * np.pi) * 2 * np.pi * 50 / 1200 ) + phi_A)
            + B * np.sin(2 * ( x + x // (4 * np.pi) * 2 * np.pi * 50 / 1242 ) + phi_B)
            + C * np.sin(2 * ( x + x // (4 * np.pi) * 2 * np.pi * 50 / 1266 ) + phi_C)
            + D * np.sin(2 * ( x - x // (4 * np.pi) * 2 * np.pi * 50 / 1200 ) + phi_D)
            + E * np.sin(2 * ( x - x // (4 * np.pi) * 2 * np.pi * 50 / 1242 ) + phi_E)
            + F * np.sin(2 * ( x - x // (4 * np.pi) * 2 * np.pi * 50 / 1266 ) + phi_F))


def superdarn_SL_SDT_full(x,
                          A,
                          B,
                          phi_A,
                          phi_B,
                          C,
                          D,
                          phi_C,
                          phi_D):
    
    """
    Fit S2 and M2 with DT and TDT included in the fits
    
    """
    
    return (  A * np.sin(2 * ( x + x // (4 * np.pi) * 2 * np.pi * 50 / 1200 ) + phi_A)
            + B * np.sin(2 * ( x + x // (4 * np.pi) * 2 * np.pi * 50 / 1242 ) + phi_B)
            + C * np.sin(1 * ( x + x // (4 * np.pi) * 2 * np.pi / 24 ) + phi_C)
            + D * np.sin(3 * ( x + x // (4 * np.pi) * 2 * np.pi / 24 ) + phi_D))
       

def dt_fitter(data,
              rads):
    
    """
    data: fit multiple of 12 hr periods
    rads: 2 * np.pi radians per 24 hrs in data interval
    
    """
    
    popt, pcov = curve_fit(dt,
                           rads,
                           data,
                           maxfev=20000)
    
    # return amplitude and phase, respectively
    return popt[0], popt[1]


def sdt_fitter(data,
               rads):
    
    """
    data: fit multiple of 12 hr periods
    rads: 2 * np.pi radians per 24 hrs in data interval
    
    """
    
    popt, pcov = curve_fit(sdt, 
                           rads, 
                           data,
                           maxfev=20000)
    
    # return amplitude and phase, respectively
    return popt[0], popt[1]


def tdt_fitter(data,
               rads):
    
    """
    data: fit multiple of 12 hr periods
    rads: 2 * np.pi radians per 24 hrs in data interval
    
    """
    
    popt, pcov = curve_fit(tdt,
                           rads,
                           data,
                           maxfev=20000)
    
    # return amplitude and phase, respectively
    return popt[0], popt[1]


def dt_sdt_tdt_fitter(data,
                      rads):
    
    """
    data: fit multiple of 12 hr periods
    rads: 2 * np.pi radians per 24 hrs in data interval
    
    """
    
    popt, pcov = curve_fit(dt_sdt_tdt, 
                           rads, 
                           data,
                           maxfev=20000)
    
    # return amplitude and phase, respectively
    return popt


def qdt_fitter(data,
               rads):
    
    """
    data: fit multiple of 12 hr periods
    rads: 2 * np.pi radians per 24 hrs in data interval
    
    """
    
    popt, pcov = curve_fit(qdt,
                           rads,
                           data,
                           maxfev=20000)
    
    # return amplitude and phase, respectively
    return popt[0], popt[1]


def q2dw_50hr_fitter(data,
                     rads):
    
    """
    data: fit multiple of 12 hr periods
    rads: 2 * np.pi radians per 24 hrs in data interval
    
    """
    
    popt, pcov = curve_fit(q2dw_50hr, 
                           rads, 
                           data,
                           bounds=([0,-np.inf],
                                   [np.inf,np.inf]),
                           maxfev=10000)
    
    # return amplitude and phase, respectively
    return popt[0], popt[1]


def navgem_dt_sdt_tdt_fitter(measurements,
                             longitudes,
                             x_dim,
                             window_w,
                             tstep):
    
    """ 
    longitudes in radians
    
    function checked 23-05-2019
    
    """
    
    # number of measurements (iterations)
    iters = int(np.shape(measurements)[0] - window_w)
    it_d = int(window_w / 2)
    
    # array to store the fit variables in 
    fit_vars = np.zeros((iters,6))
    
    # array to store wave fit in
    wave_0 = np.zeros((iters,x_dim))
    wave_1 = np.zeros((iters,x_dim)) 
    wave_2 = np.zeros((iters,x_dim)) 
    
    # x-axis on which fit is projected
    lon_fit = np.linspace(0,2*np.pi,x_dim,endpoint=False)
    
    t_steps = np.arange(it_d, iters + it_d)[::tstep]
    
    for t in t_steps:
        # select data and filter np.nans from longitude and data array
        v = measurements[(t - it_d):(t + it_d),:].flatten()    
        
        # fill longitude array with (arbitrary) 4 pi spacing between t-steps
        lons = np.empty((window_w,np.shape(longitudes)[0]))
        for k in range(window_w):
            lons[k,:] = longitudes + k * 4 * np.pi

        lons.flatten()
        
        inds = np.where(np.isnan(v))
        winds = np.delete(v,inds)
        lons = np.delete(lons,inds)
        
        # unbounded fit for better / more reliable convergence
        popt, pcov = curve_fit(navgem_dt_sdt_tdt,
                               lons, 
                               winds,
                               maxfev=20000)
                   
        # store the resolved parameters of the fit (amplitude, phase)       
        fit_vars[t - it_d,:] = popt[:]
        
        wave_0[t - it_d,:] = dt(lon_fit,popt[0],popt[3])
        wave_1[t - it_d,:] = sdt(lon_fit,popt[1],popt[4])
        wave_2[t - it_d,:] = tdt(lon_fit,popt[2],popt[5])

    return fit_vars, wave_0, wave_1, wave_2


def navgem_dt_sdt_tdt_mean_fitter(measurements,
                                  longitudes,
                                  x_dim,
                                  window_w,
                                  tstep):
    
    """ 
    longitudes in radians
    
    function checked 23-05-2019
    
    """
    
    # number of measurements (iterations)
    iters = int(np.shape(measurements)[0] - window_w)
    it_d = int(window_w / 2)
    
    # array to store the fit variables in 
    fit_vars = np.zeros((iters,6))
    
    # array to store wave fit in
    wave_0 = np.zeros((iters,x_dim))
    wave_1 = np.zeros((iters,x_dim)) 
    wave_2 = np.zeros((iters,x_dim)) 
    
    # x-axis on which fit is projected
    lon_fit = np.linspace(0,2*np.pi,x_dim,endpoint=False)
    
    t_steps = np.arange(it_d, iters + it_d)[::tstep]
    
    for t in t_steps:
        # select data and filter np.nans from longitude and data array
        v = measurements[(t - it_d):(t + it_d),:].flatten()    
        
        # fill longitude array with (arbitrary) 4 pi spacing between t-steps
        lons = np.empty((window_w,np.shape(longitudes)[0]))
        for k in range(window_w):
            lons[k,:] = longitudes + k * 4 * np.pi

        lons.flatten()
        
        inds = np.where(np.isnan(v))
        winds = np.delete(v,inds)
        lons = np.delete(lons,inds)
        
        # unbounded fit for better / more reliable convergence
        popt, pcov = curve_fit(navgem_dt_sdt_tdt_mean,
                               lons, 
                               winds,
                               maxfev=20000)
                   
        # store the resolved parameters of the fit (amplitude, phase)       
        fit_vars[t - it_d,:] = popt[:6]
        
        wave_0[t - it_d,:] = dt(lon_fit,popt[0],popt[3])
        wave_1[t - it_d,:] = sdt(lon_fit,popt[1],popt[4])
        wave_2[t - it_d,:] = tdt(lon_fit,popt[2],popt[5])

    return fit_vars, wave_0, wave_1, wave_2

def navgem_dt_sdt_tdt_mean_sw1_sw3_fitter(measurements,
                                          longitudes,
                                          x_dim,
                                          window_w,
                                          tstep):
    
    """ 
    longitudes in radians
    
    function checked 23-05-2019
    
    """
    
    # number of measurements (iterations)
    iters = int(np.shape(measurements)[0] - window_w)
    it_d = int(window_w / 2)
    
    # array to store the fit variables in 
    fit_vars = np.zeros((iters,10))
    
    # array to store wave fit in
    wave_0 = np.zeros((iters,x_dim))
    wave_1 = np.zeros((iters,x_dim)) 
    wave_2 = np.zeros((iters,x_dim)) 
    
    # x-axis on which fit is projected
    lon_fit = np.linspace(0,2*np.pi,x_dim,endpoint=False)
    
    t_steps = np.arange(it_d, iters + it_d)[::tstep]
    
    for t in t_steps:
        # select data and filter np.nans from longitude and data array
        v = measurements[(t - it_d):(t + it_d),:].flatten()    
        
        # fill longitude array with (arbitrary) 4 pi spacing between t-steps
        lons = np.empty((window_w,np.shape(longitudes)[0]))
        for k in range(window_w):
            lons[k,:] = longitudes + k * 4 * np.pi

        lons.flatten()
        
        inds = np.where(np.isnan(v))
        winds = np.delete(v,inds)
        lons = np.delete(lons,inds)
        
        # unbounded fit for better / more reliable convergence
        popt, pcov = curve_fit(navgem_dt_sdt_tdt_mean_sw1_sw3,
                               lons, 
                               winds,
                               maxfev=20000)
                   
        # store the resolved parameters of the fit (amplitude, phase)       
        fit_vars[t - it_d,:] = popt[:10]
        
        wave_0[t - it_d,:] = dt(lon_fit,popt[0],popt[3])
        wave_1[t - it_d,:] = sdt(lon_fit,popt[1],popt[4])
        wave_2[t - it_d,:] = tdt(lon_fit,popt[2],popt[5])

    return fit_vars, wave_0, wave_1, wave_2


def navgem_dt_sdt_tdt_mean_sw1_se2_fitter(measurements,
                                          longitudes,
                                          x_dim,
                                          window_w,
                                          tstep):
    
    """ 
    longitudes in radians
    
    function checked 23-05-2019
    
    """
    
    # number of measurements (iterations)
    iters = int(np.shape(measurements)[0] - window_w)
    it_d = int(window_w / 2)
    
    # array to store the fit variables in 
    fit_vars = np.zeros((iters,10))
    
    # array to store wave fit in
    wave_0 = np.zeros((iters,x_dim))
    wave_1 = np.zeros((iters,x_dim)) 
    wave_2 = np.zeros((iters,x_dim)) 
    
    # x-axis on which fit is projected
    lon_fit = np.linspace(0,2*np.pi,x_dim,endpoint=False)
    
    t_steps = np.arange(it_d, iters + it_d)[::tstep]
    
    for t in t_steps:
        # select data and filter np.nans from longitude and data array
        v = measurements[(t - it_d):(t + it_d),:].flatten()    
        
        # fill longitude array with (arbitrary) 4 pi spacing between t-steps
        lons = np.empty((window_w,np.shape(longitudes)[0]))
        for k in range(window_w):
            lons[k,:] = longitudes + k * 4 * np.pi

        lons.flatten()
        
        inds = np.where(np.isnan(v))
        winds = np.delete(v,inds)
        lons = np.delete(lons,inds)
        
        # unbounded fit for better / more reliable convergence
        popt, pcov = curve_fit(navgem_dt_sdt_tdt_mean_sw1_se2,
                               lons, 
                               winds,
                               maxfev=20000)
                   
        # store the resolved parameters of the fit (amplitude, phase)       
        fit_vars[t - it_d,:] = popt[:10]
        
        wave_0[t - it_d,:] = dt(lon_fit,popt[0],popt[3])
        wave_1[t - it_d,:] = sdt(lon_fit,popt[1],popt[4])
        wave_2[t - it_d,:] = tdt(lon_fit,popt[2],popt[5])

    return fit_vars, wave_0, wave_1, wave_2

def navgem_dt_sdt_tdt_mean_sw1_se1_fitter(measurements,
                                          longitudes,
                                          x_dim,
                                          window_w,
                                          tstep):
    
    """ 
    longitudes in radians
    
    function checked 23-05-2019
    
    """
    
    # number of measurements (iterations)
    iters = int(np.shape(measurements)[0] - window_w)
    it_d = int(window_w / 2)
    
    # array to store the fit variables in 
    fit_vars = np.zeros((iters,10))
    
    # array to store wave fit in
    wave_0 = np.zeros((iters,x_dim))
    wave_1 = np.zeros((iters,x_dim)) 
    wave_2 = np.zeros((iters,x_dim)) 
    
    # x-axis on which fit is projected
    lon_fit = np.linspace(0,2*np.pi,x_dim,endpoint=False)
    
    t_steps = np.arange(it_d, iters + it_d)[::tstep]
    
    for t in t_steps:
        # select data and filter np.nans from longitude and data array
        v = measurements[(t - it_d):(t + it_d),:].flatten()    
        
        # fill longitude array with (arbitrary) 4 pi spacing between t-steps
        lons = np.empty((window_w,np.shape(longitudes)[0]))
        for k in range(window_w):
            lons[k,:] = longitudes + k * 4 * np.pi

        lons.flatten()
        
        inds = np.where(np.isnan(v))
        winds = np.delete(v,inds)
        lons = np.delete(lons,inds)
        
        # unbounded fit for better / more reliable convergence
        popt, pcov = curve_fit(navgem_dt_sdt_tdt_mean_sw1_se1,
                               lons, 
                               winds,
                               maxfev=20000)
                   
        # store the resolved parameters of the fit (amplitude, phase)       
        fit_vars[t - it_d,:] = popt[:10]
        
        wave_0[t - it_d,:] = dt(lon_fit,popt[0],popt[3])
        wave_1[t - it_d,:] = sdt(lon_fit,popt[1],popt[4])
        wave_2[t - it_d,:] = tdt(lon_fit,popt[2],popt[5])

    return fit_vars, wave_0, wave_1, wave_2


def navgem_dt_sdt_tdt_mean_sw1_fitter(measurements,
                                      longitudes,
                                      x_dim,
                                      window_w,
                                      tstep):
    
    """ 
    longitudes in radians
    
    function checked 23-05-2019
    
    """
    
    # number of measurements (iterations)
    iters = int(np.shape(measurements)[0] - window_w)
    it_d = int(window_w / 2)
    
    # array to store the fit variables in 
    fit_vars = np.zeros((iters,8))
    
    # array to store wave fit in
    wave_0 = np.zeros((iters,x_dim))
    wave_1 = np.zeros((iters,x_dim)) 
    wave_2 = np.zeros((iters,x_dim)) 
    
    # x-axis on which fit is projected
    lon_fit = np.linspace(0,2*np.pi,x_dim,endpoint=False)
    
    t_steps = np.arange(it_d, iters + it_d)[::tstep]
    
    for t in t_steps:
        # select data and filter np.nans from longitude and data array
        v = measurements[(t - it_d):(t + it_d),:].flatten()    
        
        # fill longitude array with (arbitrary) 4 pi spacing between t-steps
        lons = np.empty((window_w,np.shape(longitudes)[0]))
        for k in range(window_w):
            lons[k,:] = longitudes + k * 4 * np.pi

        lons.flatten()
        
        inds = np.where(np.isnan(v))
        winds = np.delete(v,inds)
        lons = np.delete(lons,inds)
        
        # unbounded fit for better / more reliable convergence
        popt, pcov = curve_fit(navgem_dt_sdt_tdt_mean_sw1,
                               lons, 
                               winds,
                               maxfev=20000)
                   
        # store the resolved parameters of the fit (amplitude, phase)       
        fit_vars[t - it_d,:] = popt[:8]
        
        wave_0[t - it_d,:] = dt(lon_fit,popt[0],popt[3])
        wave_1[t - it_d,:] = sdt(lon_fit,popt[1],popt[4])
        wave_2[t - it_d,:] = tdt(lon_fit,popt[2],popt[5])

    return fit_vars, wave_0, wave_1, wave_2


def navgem_nonmigrating_fitter(measurements,
                               longitudes,
                               x_dim,
                               window_w,
                               tstep):
    
    """ 
    longitudes in radians
    
    function checked 23-05-2019
    
    """
    
    # number of measurements (iterations)
    iters = int(np.shape(measurements)[0] - window_w)
    it_d = int(window_w / 2)
    
    # array to store the fit variables in 
    fit_vars = np.zeros((iters,6))
    
    # array to store wave fit in
    wave_0 = np.zeros((iters,x_dim))
    wave_1 = np.zeros((iters,x_dim)) 
    wave_2 = np.zeros((iters,x_dim)) 
    
    # x-axis on which fit is projected
    lon_fit = np.linspace(0,2*np.pi,x_dim,endpoint=False)
    
    t_steps = np.arange(it_d, iters + it_d)[::tstep]
    
    for t in t_steps:
        # select data and filter np.nans from longitude and data array
        v = measurements[(t - it_d):(t + it_d),:].flatten()    
        
        # fill longitude array with (arbitrary) 4 pi spacing between t-steps
        lons = np.empty((window_w,np.shape(longitudes)[0]))
        for k in range(window_w):
            lons[k,:] = longitudes + k * 4 * np.pi

        lons.flatten()
        
        inds = np.where(np.isnan(v))
        winds = np.delete(v,inds)
        lons = np.delete(lons,inds)
        
        # unbounded fit for better / more reliable convergence
        popt, pcov = curve_fit(navgem_nonmigrating,
                               lons, 
                               winds,
                               maxfev=20000)
                   
        # store the resolved parameters of the fit (amplitude, phase)       
        fit_vars[t - it_d,:] = popt[:6]
        
        wave_0[t - it_d,:] = dt(lon_fit,popt[0],popt[3])
        wave_1[t - it_d,:] = sdt(lon_fit,popt[1],popt[4])
        wave_2[t - it_d,:] = tdt(lon_fit,popt[2],popt[5])

    return fit_vars, wave_0, wave_1, wave_2



def superdarn_dt_sdt_tdt_fitter(measurements,
                                longitudes,
                                x_dim,
                                window_w,
                                tstep,
                                data_threshold):
    
    """ 
    longitudes in radians
    
    """
    
    # number of measurements (iterations)
    iters = int(np.shape(measurements)[0] - window_w)
    it_d = int(window_w / 2)
    
    # array to store the fit variables in 
    fit_vars = np.zeros((iters,6))
    fit_data = np.empty(iters)
    
    # array to store wave fit in
    wave_0 = np.zeros((iters,x_dim))
    wave_1 = np.zeros((iters,x_dim)) 
    wave_2 = np.zeros((iters,x_dim)) 
    
    # x-axis on which fit is projected
    lon_fit = np.linspace(0,2 * np.pi,x_dim,endpoint=False)
    
    t_steps = np.arange(it_d, iters + it_d)[::tstep]
    
    for t in t_steps:
        # select data and filter np.nans from longitude and data array
        v = measurements[(t - it_d):(t + it_d),:].flatten()
        
        # number of data points in fit
        fit_data[t - it_d] = int(window_w * 10 - np.isnan(v).sum())
        
        # fill longitude array with (arbitrary) 4pi spacing between t-steps
        lons = np.empty((window_w,np.shape(longitudes)[0]))
        for k in range(window_w):
            lons[k,:] = longitudes + k * 4 * np.pi

        lons.flatten()
        
        # if more nan values than threshold, reject fit
        threshold = int(window_w * 10 - data_threshold)
        
        if int(np.isnan(v).sum()) > threshold:
            # store the resolved parameters of the fit (amplitude, phase)
            fit_vars[t - it_d,:] = np.nan
            
            wave_0[t - it_d,:] = np.nan
            wave_1[t - it_d,:] = np.nan
            wave_2[t - it_d,:] = np.nan
            
            # do not take non-fits into account by setting to nan
            fit_data[t - it_d] = np.nan 
        
        else:
            inds = np.where(np.isnan(v))
            winds = np.delete(v,inds)
            lons = np.delete(lons,inds)
                    
            # bounds force the ampltitude to be positive
            popt, pcov = curve_fit(superdarn_dt_sdt_tdt,
                                   lons, 
                                   winds,
                                   maxfev=20000)
                       
            # store the resolved parameters of the fit (amplitude, phase)
            fit_vars[t - it_d,:] = popt[:]
            
            wave_0[t - it_d,:] = dt(lon_fit,popt[0],popt[3])
            wave_1[t - it_d,:] = sdt(lon_fit,popt[1],popt[4])
            wave_2[t - it_d,:] = tdt(lon_fit,popt[2],popt[5])
            
    return fit_vars, wave_0, wave_1, wave_2, fit_data


def superdarn_dt_sdt_tdt_mean_fitter(measurements,
                                     longitudes,
                                     x_dim,
                                     window_w,
                                     tstep,
                                     data_threshold):
    
    """ 
    longitudes in radians
    
    """
    
    # number of measurements (iterations)
    iters = int(np.shape(measurements)[0] - window_w)
    it_d = int(window_w / 2)
    
    # array to store the fit variables in 
    fit_vars = np.zeros((iters,6))
    fit_data = np.empty(iters)
    
    # array to store wave fit in
    wave_0 = np.zeros((iters,x_dim))
    wave_1 = np.zeros((iters,x_dim)) 
    wave_2 = np.zeros((iters,x_dim)) 
    
    # x-axis on which fit is projected
    lon_fit = np.linspace(0,2 * np.pi,x_dim,endpoint=False)
    
    t_steps = np.arange(it_d, iters + it_d)[::tstep]
    
    print('whaddup nurd')
    
    for t in t_steps:
        # select data and filter np.nans from longitude and data array
        v = measurements[(t - it_d):(t + it_d),:].flatten()
        
        # number of data points in fit
        fit_data[t - it_d] = int(window_w * 10 - np.isnan(v).sum())
        
        # fill longitude array with (arbitrary) 4pi spacing between t-steps
        lons = np.empty((window_w,np.shape(longitudes)[0]))
        for k in range(window_w):
            lons[k,:] = longitudes + k * 4 * np.pi

        lons.flatten()
        
        print(np.shape(lons))
        
        # if more nan values than threshold, reject fit
        threshold = int(window_w * 10 - data_threshold)
        
        if int(np.isnan(v).sum()) > threshold:
            # store the resolved parameters of the fit (amplitude, phase)
            fit_vars[t - it_d,:] = np.nan
            
            wave_0[t - it_d,:] = np.nan
            wave_1[t - it_d,:] = np.nan
            wave_2[t - it_d,:] = np.nan
            
            # do not take non-fits into account by setting to nan
            fit_data[t - it_d] = np.nan 
        
        else:
            inds = np.where(np.isnan(v))
            winds = np.delete(v,inds)
            lons = np.delete(lons,inds)
                    
            # bounds force the ampltitude to be positive
            popt, pcov = curve_fit(superdarn_dt_sdt_tdt_mean,
                                   lons, 
                                   winds,
                                   maxfev=20000)
                       
            # store the resolved parameters of the fit (amplitude, phase)
            fit_vars[t - it_d,:] = popt[:6]
            
            wave_0[t - it_d,:] = dt(lon_fit,popt[0],popt[3])
            wave_1[t - it_d,:] = sdt(lon_fit,popt[1],popt[4])
            wave_2[t - it_d,:] = tdt(lon_fit,popt[2],popt[5])
            
    return fit_vars, wave_0, wave_1, wave_2, fit_data


def superdarn_dt_m2_tdt_mean_fitter(measurements,
                                    longitudes,
                                    x_dim,
                                    window_w,
                                    tstep,
                                    data_threshold):
    
    """ 
    longitudes in radians
    
    """
    
    # number of measurements (iterations)
    iters = int(np.shape(measurements)[0] - window_w)
    it_d = int(window_w / 2)
    
    # array to store the fit variables in 
    fit_vars = np.zeros((iters,6))
    fit_data = np.empty(iters)
    
    # array to store wave fit in
    wave_0 = np.zeros((iters,x_dim))
    wave_1 = np.zeros((iters,x_dim)) 
    wave_2 = np.zeros((iters,x_dim)) 
    
    # x-axis on which fit is projected
    lon_fit = np.linspace(0,2 * np.pi,x_dim,endpoint=False)
    
    t_steps = np.arange(it_d, iters + it_d)[::tstep]
    
    for t in t_steps:
        # select data and filter np.nans from longitude and data array
        v = measurements[(t - it_d):(t + it_d),:].flatten()
        
        # number of data points in fit
        fit_data[t - it_d] = int(window_w * 10 - np.isnan(v).sum())
        
        # fill longitude array with (arbitrary) 4pi spacing between t-steps
        lons = np.empty((window_w,np.shape(longitudes)[0]))
        for k in range(window_w):
            lons[k,:] = longitudes + k * 4 * np.pi

        lons.flatten()
        
        # if more nan values than threshold, reject fit
        threshold = int(window_w * 10 - data_threshold)
        
        if int(np.isnan(v).sum()) > threshold:
            # store the resolved parameters of the fit (amplitude, phase)
            fit_vars[t - it_d,:] = np.nan
            
            wave_0[t - it_d,:] = np.nan
            wave_1[t - it_d,:] = np.nan
            wave_2[t - it_d,:] = np.nan
            
            # do not take non-fits into account by setting to nan
            fit_data[t - it_d] = np.nan 
        
        else:
            inds = np.where(np.isnan(v))
            winds = np.delete(v,inds)
            lons = np.delete(lons,inds)
                    
            # bounds force the ampltitude to be positive
            popt, pcov = curve_fit(superdarn_dt_m2_tdt_mean,
                                   lons, 
                                   winds,
                                   maxfev=20000)
                       
            # store the resolved parameters of the fit (amplitude, phase)
            fit_vars[t - it_d,:] = popt[:6]
            
            wave_0[t - it_d,:] = dt(lon_fit,popt[0],popt[3])
            wave_1[t - it_d,:] = sdt(lon_fit,popt[1],popt[4])
            wave_2[t - it_d,:] = tdt(lon_fit,popt[2],popt[5])
            
    return fit_vars, wave_0, wave_1, wave_2, fit_data



def superdarn_dt_sdt_tdt_mean_sw1_fitter(measurements,
                                         longitudes,
                                         x_dim,
                                         window_w,
                                         tstep,
                                         data_threshold):
    
    """ 
    longitudes in radians
    
    """
    
    # number of measurements (iterations)
    iters = int(np.shape(measurements)[0] - window_w)
    it_d = int(window_w / 2)
    
    # array to store the fit variables in 
    fit_vars = np.zeros((iters,8))
    fit_data = np.empty(iters)
    
    # array to store wave fit in
    wave_0 = np.zeros((iters,x_dim))
    wave_1 = np.zeros((iters,x_dim)) 
    wave_2 = np.zeros((iters,x_dim)) 
    
    # x-axis on which fit is projected
    lon_fit = np.linspace(0,2 * np.pi,x_dim,endpoint=False)
    
    t_steps = np.arange(it_d, iters + it_d)[::tstep]
    
    for t in t_steps:
        # select data and filter np.nans from longitude and data array
        v = measurements[(t - it_d):(t + it_d),:].flatten()
        
        # number of data points in fit
        fit_data[t - it_d] = int(window_w * 10 - np.isnan(v).sum())
        
        # fill longitude array with (arbitrary) 4pi spacing between t-steps
        lons = np.empty((window_w,np.shape(longitudes)[0]))
        for k in range(window_w):
            lons[k,:] = longitudes + k * 4 * np.pi

        lons.flatten()
        
        # if more nan values than threshold, reject fit
        threshold = int(window_w * 10 - data_threshold)
        
        if int(np.isnan(v).sum()) > threshold:
            # store the resolved parameters of the fit (amplitude, phase)
            fit_vars[t - it_d,:] = np.nan
            
            wave_0[t - it_d,:] = np.nan
            wave_1[t - it_d,:] = np.nan
            wave_2[t - it_d,:] = np.nan
            
            # do not take non-fits into account by setting to nan
            fit_data[t - it_d] = np.nan 
        
        else:
            inds = np.where(np.isnan(v))
            winds = np.delete(v,inds)
            lons = np.delete(lons,inds)
                    
            # bounds force the ampltitude to be positive
            popt, pcov = curve_fit(superdarn_dt_sdt_tdt_mean_sw1,
                                   lons, 
                                   winds,
                                   maxfev=20000)
                       
            # store the resolved parameters of the fit (amplitude, phase)
            fit_vars[t - it_d,:] = popt[:8]
            
            wave_0[t - it_d,:] = dt(lon_fit,popt[0],popt[3])
            wave_1[t - it_d,:] = sdt(lon_fit,popt[1],popt[4])
            wave_2[t - it_d,:] = tdt(lon_fit,popt[2],popt[5])
            
    return fit_vars, wave_0, wave_1, wave_2, fit_data

def superdarn_dt_sdt_tdt_mean_sw1_sw3_fitter(measurements,
                                             longitudes,
                                             x_dim,
                                             window_w,
                                             tstep,
                                             data_threshold):
    
    """ 
    longitudes in radians
    
    """
    
    # number of measurements (iterations)
    iters = int(np.shape(measurements)[0] - window_w)
    it_d = int(window_w / 2)
    
    # array to store the fit variables in 
    fit_vars = np.zeros((iters,10))
    fit_data = np.empty(iters)
    
    # array to store wave fit in
    wave_0 = np.zeros((iters,x_dim))
    wave_1 = np.zeros((iters,x_dim)) 
    wave_2 = np.zeros((iters,x_dim)) 
    
    # x-axis on which fit is projected
    lon_fit = np.linspace(0,2 * np.pi,x_dim,endpoint=False)
    
    t_steps = np.arange(it_d, iters + it_d)[::tstep]
    
    for t in t_steps:
        # select data and filter np.nans from longitude and data array
        v = measurements[(t - it_d):(t + it_d),:].flatten()
        
        # number of data points in fit
        fit_data[t - it_d] = int(window_w * 10 - np.isnan(v).sum())
        
        # fill longitude array with (arbitrary) 4pi spacing between t-steps
        lons = np.empty((window_w,np.shape(longitudes)[0]))
        for k in range(window_w):
            lons[k,:] = longitudes + k * 4 * np.pi

        lons.flatten()
        
        # if more nan values than threshold, reject fit
        threshold = int(window_w * 10 - data_threshold)
        
        if int(np.isnan(v).sum()) > threshold:
            # store the resolved parameters of the fit (amplitude, phase)
            fit_vars[t - it_d,:] = np.nan
            
            wave_0[t - it_d,:] = np.nan
            wave_1[t - it_d,:] = np.nan
            wave_2[t - it_d,:] = np.nan
            
            # do not take non-fits into account by setting to nan
            fit_data[t - it_d] = np.nan 
        
        else:
            inds = np.where(np.isnan(v))
            winds = np.delete(v,inds)
            lons = np.delete(lons,inds)
                    
            # bounds force the ampltitude to be positive
            popt, pcov = curve_fit(superdarn_dt_sdt_tdt_mean_sw1_sw3,
                                   lons, 
                                   winds,
                                   maxfev=20000)
                       
            # store the resolved parameters of the fit (amplitude, phase)
            fit_vars[t - it_d,:] = popt[:10]
            
            wave_0[t - it_d,:] = dt(lon_fit,popt[0],popt[3])
            wave_1[t - it_d,:] = sdt(lon_fit,popt[1],popt[4])
            wave_2[t - it_d,:] = tdt(lon_fit,popt[2],popt[5])
            
    return fit_vars, wave_0, wave_1, wave_2, fit_data

def superdarn_dt_sdt_tdt_mean_sw1_se2_fitter(measurements,
                                             longitudes,
                                             x_dim,
                                             window_w,
                                             tstep,
                                             data_threshold):
    
    """ 
    longitudes in radians
    
    """
    
    # number of measurements (iterations)
    iters = int(np.shape(measurements)[0] - window_w)
    it_d = int(window_w / 2)
    
    # array to store the fit variables in 
    fit_vars = np.zeros((iters,10))
    fit_data = np.empty(iters)
    
    # array to store wave fit in
    wave_0 = np.zeros((iters,x_dim))
    wave_1 = np.zeros((iters,x_dim)) 
    wave_2 = np.zeros((iters,x_dim)) 
    
    # x-axis on which fit is projected
    lon_fit = np.linspace(0,2 * np.pi,x_dim,endpoint=False)
    
    t_steps = np.arange(it_d, iters + it_d)[::tstep]
    
    for t in t_steps:
        # select data and filter np.nans from longitude and data array
        v = measurements[(t - it_d):(t + it_d),:].flatten()
        
        # number of data points in fit
        fit_data[t - it_d] = int(window_w * 10 - np.isnan(v).sum())
        
        # fill longitude array with (arbitrary) 4pi spacing between t-steps
        lons = np.empty((window_w,np.shape(longitudes)[0]))
        for k in range(window_w):
            lons[k,:] = longitudes + k * 4 * np.pi

        lons.flatten()
        
        # if more nan values than threshold, reject fit
        threshold = int(window_w * 10 - data_threshold)
        
        if int(np.isnan(v).sum()) > threshold:
            # store the resolved parameters of the fit (amplitude, phase)
            fit_vars[t - it_d,:] = np.nan
            
            wave_0[t - it_d,:] = np.nan
            wave_1[t - it_d,:] = np.nan
            wave_2[t - it_d,:] = np.nan
            
            # do not take non-fits into account by setting to nan
            fit_data[t - it_d] = np.nan 
        
        else:
            inds = np.where(np.isnan(v))
            winds = np.delete(v,inds)
            lons = np.delete(lons,inds)
                    
            # bounds force the ampltitude to be positive
            popt, pcov = curve_fit(superdarn_dt_sdt_tdt_mean_sw1_se2,
                                   lons, 
                                   winds,
                                   maxfev=20000)
                       
            # store the resolved parameters of the fit (amplitude, phase)
            fit_vars[t - it_d,:] = popt[:10]
            
            wave_0[t - it_d,:] = dt(lon_fit,popt[0],popt[3])
            wave_1[t - it_d,:] = sdt(lon_fit,popt[1],popt[4])
            wave_2[t - it_d,:] = tdt(lon_fit,popt[2],popt[5])
            
    return fit_vars, wave_0, wave_1, wave_2, fit_data


def superdarn_dt_sdt_tdt_mean_sw1_se1_fitter(measurements,
                                             longitudes,
                                             x_dim,
                                             window_w,
                                             tstep,
                                             data_threshold):
    
    """ 
    longitudes in radians
    
    """
    
    # number of measurements (iterations)
    iters = int(np.shape(measurements)[0] - window_w)
    it_d = int(window_w / 2)
    
    # array to store the fit variables in 
    fit_vars = np.zeros((iters,10))
    fit_data = np.empty(iters)
    
    # array to store wave fit in
    wave_0 = np.zeros((iters,x_dim))
    wave_1 = np.zeros((iters,x_dim)) 
    wave_2 = np.zeros((iters,x_dim)) 
    
    # x-axis on which fit is projected
    lon_fit = np.linspace(0,2 * np.pi,x_dim,endpoint=False)
    
    t_steps = np.arange(it_d, iters + it_d)[::tstep]
    
    for t in t_steps:
        # select data and filter np.nans from longitude and data array
        v = measurements[(t - it_d):(t + it_d),:].flatten()
        
        # number of data points in fit
        fit_data[t - it_d] = int(window_w * 10 - np.isnan(v).sum())
        
        # fill longitude array with (arbitrary) 4pi spacing between t-steps
        lons = np.empty((window_w,np.shape(longitudes)[0]))
        for k in range(window_w):
            lons[k,:] = longitudes + k * 4 * np.pi

        lons.flatten()
        
        # if more nan values than threshold, reject fit
        threshold = int(window_w * 10 - data_threshold)
        
        if int(np.isnan(v).sum()) > threshold:
            # store the resolved parameters of the fit (amplitude, phase)
            fit_vars[t - it_d,:] = np.nan
            
            wave_0[t - it_d,:] = np.nan
            wave_1[t - it_d,:] = np.nan
            wave_2[t - it_d,:] = np.nan
            
            # do not take non-fits into account by setting to nan
            fit_data[t - it_d] = np.nan 
        
        else:
            inds = np.where(np.isnan(v))
            winds = np.delete(v,inds)
            lons = np.delete(lons,inds)
                    
            # bounds force the ampltitude to be positive
            popt, pcov = curve_fit(superdarn_dt_sdt_tdt_mean_sw1_se1,
                                   lons, 
                                   winds,
                                   maxfev=20000)
                       
            # store the resolved parameters of the fit (amplitude, phase)
            fit_vars[t - it_d,:] = popt[:10]
            
            wave_0[t - it_d,:] = dt(lon_fit,popt[0],popt[3])
            wave_1[t - it_d,:] = sdt(lon_fit,popt[1],popt[4])
            wave_2[t - it_d,:] = tdt(lon_fit,popt[2],popt[5])
            
    return fit_vars, wave_0, wave_1, wave_2, fit_data

def superdarn_dt_sdt_tdt_mean_fitter(measurements,
                                     longitudes,
                                     x_dim,
                                     window_w,
                                     tstep,
                                     data_threshold):
    
    """ 
    longitudes in radians
    
    """
    
    # number of measurements (iterations)
    iters = int(np.shape(measurements)[0] - window_w)
    it_d = int(window_w / 2)
    
    # array to store the fit variables in 
    fit_vars = np.zeros((iters,6))
    fit_data = np.empty(iters)
    
    # array to store wave fit in
    wave_0 = np.zeros((iters,x_dim))
    wave_1 = np.zeros((iters,x_dim)) 
    wave_2 = np.zeros((iters,x_dim)) 
    
    # x-axis on which fit is projected
    lon_fit = np.linspace(0,2 * np.pi,x_dim,endpoint=False)
    
    t_steps = np.arange(it_d, iters + it_d)[::tstep]
    
    for t in t_steps:
        # select data and filter np.nans from longitude and data array
        v = measurements[(t - it_d):(t + it_d),:].flatten()
        
        # number of data points in fit
        fit_data[t - it_d] = int(window_w * 10 - np.isnan(v).sum())
        
        # fill longitude array with (arbitrary) 4pi spacing between t-steps
        lons = np.empty((window_w,np.shape(longitudes)[0]))
        for k in range(window_w):
            lons[k,:] = longitudes + k * 4 * np.pi

        lons.flatten()
            
        # if more nan values than threshold, reject fit
        threshold = int(window_w * 10 - data_threshold)
        
        if int(np.isnan(v).sum()) > threshold:
            # store the resolved parameters of the fit (amplitude, phase)
            fit_vars[t - it_d,:] = np.nan
            
            wave_0[t - it_d,:] = np.nan
            wave_1[t - it_d,:] = np.nan
            wave_2[t - it_d,:] = np.nan
            
            # do not take non-fits into account by setting to nan
            fit_data[t - it_d] = np.nan 
        
        else:
            inds = np.where(np.isnan(v))
            winds = np.delete(v,inds)
            lons = np.delete(lons,inds)
                    
            # bounds force the ampltitude to be positive
            popt, pcov = curve_fit(superdarn_dt_sdt_tdt_mean,
                                   lons, 
                                   winds,
                                   maxfev=20000)
                       
            # store the resolved parameters of the fit (amplitude, phase)
            fit_vars[t - it_d,:] = popt[:6]
            
            wave_0[t - it_d,:] = dt(lon_fit,popt[0],popt[3])
            wave_1[t - it_d,:] = sdt(lon_fit,popt[1],popt[4])
            wave_2[t - it_d,:] = tdt(lon_fit,popt[2],popt[5])
            
    return fit_vars, wave_0, wave_1, wave_2, fit_data


def superdarn_dt_sdt_tdt_fitter_flip(measurements,
                                     longitudes,
                                     x_dim,
                                     window_w,
                                     tstep,
                                     data_threshold):
    
    """ 
    longitudes in radians
    
    """
    
    # number of measurements (iterations)
    iters = int(np.shape(measurements)[0] - window_w)
    it_d = int(window_w / 2)
    
    # array to store the fit variables in 
    fit_vars = np.zeros((iters,6))
    fit_data = np.empty(iters)
    
    # array to store wave fit in
    wave_0 = np.zeros((iters,x_dim))
    wave_1 = np.zeros((iters,x_dim)) 
    wave_2 = np.zeros((iters,x_dim)) 
    
    # x-axis on which fit is projected
    lon_fit = np.linspace(0,2 * np.pi,x_dim,endpoint=False)
    
    t_steps = np.arange(it_d, iters + it_d)[::tstep]
    
    for t in t_steps:
        # select data and filter np.nans from longitude and data array
        v = measurements[(t - it_d):(t + it_d),:].flatten()
        
        # number of data points in fit
        fit_data[t - it_d] = int(window_w * 10 - np.isnan(v).sum())
        
        # fill longitude array with (arbitrary) 4pi spacing between t-steps
        lons = np.empty((window_w,np.shape(longitudes)[0]))
        for k in range(window_w):
            lons[k,:] = longitudes + k * 4 * np.pi

        lons.flatten()
        
        # if more nan values than threshold, reject fit
        threshold = int(window_w * 10 - data_threshold)
        
        if int(np.isnan(v).sum()) > threshold:
            # store the resolved parameters of the fit (amplitude, phase)
            fit_vars[t - it_d,:] = np.nan
            
            wave_0[t - it_d,:] = np.nan
            wave_1[t - it_d,:] = np.nan
            wave_2[t - it_d,:] = np.nan
            
            # do not take non-fits into account by setting to nan
            fit_data[t - it_d] = np.nan 
        
        else:
            inds = np.where(np.isnan(v))
            winds = np.delete(v,inds)
            lons = np.delete(lons,inds)
                    
            # bounds force the ampltitude to be positive
            popt, pcov = curve_fit(superdarn_dt_sdt_tdt_flip,
                                   lons, 
                                   winds,
                                   maxfev=20000)
                       
            # store the resolved parameters of the fit (amplitude, phase)
            fit_vars[t - it_d,:] = popt[:]
            
            wave_0[t - it_d,:] = dt(lon_fit,popt[0],popt[3])
            wave_1[t - it_d,:] = sdt(lon_fit,popt[1],popt[4])
            wave_2[t - it_d,:] = tdt(lon_fit,popt[2],popt[5])
            
    return fit_vars, wave_0, wave_1, wave_2, fit_data


def superdarn_SLT_sdt_fit(measurements,
                          longitudes,
                          x_dim,
                          window_w,
                          tstep,
                          data_threshold,
                          only_lunar=False,
                          only_solar=False):
    
    """ 
    longitudes must be in radians
    
    """
    
    # number of measurements (iterations)
    iters = int(np.shape(measurements)[0] - window_w)
    it_d = int(window_w / 2)
    
    # array to store the fit variables in 
    fit_vars = np.zeros((iters,4))
    fit_data = np.empty(iters)
    
    # array to store wave fit in
    wave_0 = np.zeros((iters,x_dim))
    wave_1 = np.zeros((iters,x_dim))  
    
    # x-axis on which fit is projected
    lon_fit = np.linspace(0,2 * np.pi,x_dim,endpoint=False)
    
    t_steps = np.arange(it_d, iters + it_d)[::tstep]
    
    for t in t_steps:
        # select data and filter np.nans from longitude and data array
        v = measurements[(t - it_d):(t + it_d),:].flatten()
        
        # number of data points in fit
        fit_data[t - it_d] = int(window_w * 10 - np.isnan(v).sum())
        
        # fill longitude array with (arbitrary) 4pi spacing between t-steps
        lons = np.empty((window_w,np.shape(longitudes)[0]))
        for k in range(window_w):
            lons[k,:] = longitudes + k * 4 * np.pi

        lons.flatten()
        
        # if more nan values than threshold, reject fit
        threshold = int(window_w * 10 - data_threshold)
        
        if int(np.isnan(v).sum()) > threshold:
            # store the resolved parameters of the fit (amplitude, phase)
            fit_vars[t - it_d,:] = np.nan
            
            wave_0[t - it_d,:] = np.nan
            wave_1[t - it_d,:] = np.nan
            
            # do not take non-fits into account by setting to nan
            fit_data[t - it_d] = np.nan 
        
        else:
            inds = np.where(np.isnan(v))
            winds = np.delete(v,inds)
            lons = np.delete(lons,inds)
                    
            if only_lunar:
                popt, pcov = curve_fit(superdarn_SL_SDT,
                                       lons,
                                       winds,
                                       maxfev=20000)
            elif only_solar:
                popt, pcov = curve_fit(superdarn_SL_SDT,
                                       lons,
                                       winds,
                                       maxfev=20000)
            else:
                popt, pcov = curve_fit(superdarn_SL_SDT,
                                       lons,
                                       winds,
                                       maxfev=20000)
                
            # store the resolved parameters of the fit (amplitude, phase)
            fit_vars[t - it_d,:] = popt[:]
            
            wave_0[t - it_d,:] = sdt(lon_fit,popt[0],popt[2])
            wave_1[t - it_d,:] = sdt(lon_fit,popt[1],popt[3])
            
    return fit_vars, wave_0, wave_1, fit_data


def navgem_SLT_sdt_fit(measurements,
                       longitudes,
                       x_dim,
                       window_w,
                       tstep,
                       data_threshold):
    
    """ 
    longitudes must be in radians
    
    """
    
    # number of measurements (iterations)
    iters = int(np.shape(measurements)[0] - window_w)
    it_d = int(window_w / 2)
    
    # array to store the fit variables in 
    fit_vars = np.zeros((iters,4))
    fit_data = np.empty(iters)
    
    # array to store wave fit in
    wave_0 = np.zeros((iters,x_dim))
    wave_1 = np.zeros((iters,x_dim)) 
    
    # x-axis on which fit is projected
    lon_fit = np.linspace(0,2 * np.pi,x_dim,endpoint=False)
    
    t_steps = np.arange(it_d, iters + it_d)[::tstep]
    
    for t in t_steps:
        # select data and filter np.nans from longitude and data array
        v = measurements[(t - it_d):(t + it_d),:].flatten()
        
        # number of data points in fit
        fit_data[t - it_d] = int(window_w * 10 - np.isnan(v).sum())
        
        # fill longitude array with (arbitrary) 4pi spacing between t-steps
        lons = np.empty((window_w,np.shape(longitudes)[0]))
        for k in range(window_w):
            lons[k,:] = longitudes + k * 4 * np.pi

        lons.flatten()
        
        # if more nan values than threshold, reject fit
        threshold = int(window_w * 10 - data_threshold)
        
        if int(np.isnan(v).sum()) > threshold:
            # store the resolved parameters of the fit (amplitude, phase)
            fit_vars[t - it_d,:] = np.nan
            
            wave_0[t - it_d,:] = np.nan
            wave_1[t - it_d,:] = np.nan
            
            # do not take non-fits into account by setting to nan
            fit_data[t - it_d] = np.nan 
        
        else:
            inds = np.where(np.isnan(v))
            winds = np.delete(v,inds)
            lons = np.delete(lons,inds)
                    
            popt, pcov = curve_fit(navgem_SL_SDT,
                                   lons,
                                   winds,
                                   maxfev=20000)
                
            # store the resolved parameters of the fit (amplitude, phase)
            fit_vars[t - it_d,:] = popt[:]
            
            wave_0[t - it_d,:] = sdt(lon_fit,popt[0],popt[2])
            wave_1[t - it_d,:] = sdt(lon_fit,popt[1],popt[3])
            
    return fit_vars, wave_0, wave_1, fit_data


def navgem_SLT_m2n2_fit(measurements,
                        longitudes,
                        x_dim,
                        window_w,
                        tstep,
                        data_threshold):
    
    """ 
    longitudes must be in radians
    
    """
    
    # number of measurements (iterations)
    iters = int(np.shape(measurements)[0] - window_w)
    it_d = int(window_w / 2)
    
    # array to store the fit variables in 
    fit_vars = np.zeros((iters,6))
    fit_data = np.empty(iters)
    
    # array to store wave fit in
    wave_0 = np.zeros((iters,x_dim))
    wave_1 = np.zeros((iters,x_dim)) 
    
    # x-axis on which fit is projected
    lon_fit = np.linspace(0,2 * np.pi,x_dim,endpoint=False)
    
    t_steps = np.arange(it_d, iters + it_d)[::tstep]
    
    for t in t_steps:
        # select data and filter np.nans from longitude and data array
        v = measurements[(t - it_d):(t + it_d),:].flatten()
        
        # number of data points in fit
        fit_data[t - it_d] = int(window_w * 10 - np.isnan(v).sum())
        
        # fill longitude array with (arbitrary) 4pi spacing between t-steps
        lons = np.empty((window_w,np.shape(longitudes)[0]))
        for k in range(window_w):
            lons[k,:] = longitudes + k * 4 * np.pi

        lons.flatten()
        
        # if more nan values than threshold, reject fit
        threshold = int(window_w * 10 - data_threshold)
        
        if int(np.isnan(v).sum()) > threshold:
            # store the resolved parameters of the fit (amplitude, phase)
            fit_vars[t - it_d,:] = np.nan
            
            wave_0[t - it_d,:] = np.nan
            wave_1[t - it_d,:] = np.nan
            
            # do not take non-fits into account by setting to nan
            fit_data[t - it_d] = np.nan 
        
        else:
            inds = np.where(np.isnan(v))
            winds = np.delete(v,inds)
            lons = np.delete(lons,inds)
                    
            popt, pcov = curve_fit(navgem_SL_m2n2,
                                   lons,
                                   winds,
                                   maxfev=20000)
                
            # store the resolved parameters of the fit (amplitude, phase)
            fit_vars[t - it_d,:] = popt[:]
            
            wave_0[t - it_d,:] = sdt(lon_fit,popt[0],popt[2])
            wave_1[t - it_d,:] = sdt(lon_fit,popt[1],popt[3])
            
    return fit_vars, wave_0, wave_1, fit_data


def superdarn_SLT_sdt_full_fit(measurements,
                               longitudes,
                               x_dim,
                               window_w,
                               tstep,
                               data_threshold,
                               only_lunar=False,
                               only_solar=False):
    
    """ 
    longitudes must be in radians
    
    """
    
    # number of measurements (iterations)
    iters = int(np.shape(measurements)[0] - window_w)
    it_d = int(window_w / 2)
    
    # array to store the fit variables in 
    fit_vars = np.zeros((iters,8))
    fit_data = np.empty(iters)
    
    # array to store wave fit in
    wave_0 = np.zeros((iters,x_dim))
    wave_1 = np.zeros((iters,x_dim)) 
    wave_2 = np.zeros((iters,x_dim)) 
    
    # x-axis on which fit is projected
    lon_fit = np.linspace(0,2 * np.pi,x_dim,endpoint=False)
    
    t_steps = np.arange(it_d, iters + it_d)[::tstep]
    
    for t in t_steps:
        # select data and filter np.nans from longitude and data array
        v = measurements[(t - it_d):(t + it_d),:].flatten()
        
        # number of data points in fit
        fit_data[t - it_d] = int(window_w * 10 - np.isnan(v).sum())
        
        # fill longitude array with (arbitrary) 4pi spacing between t-steps
        lons = np.empty((window_w,np.shape(longitudes)[0]))
        for k in range(window_w):
            lons[k,:] = longitudes + k * 4 * np.pi

        lons.flatten()
        
        # if more nan values than threshold, reject fit
        threshold = int(window_w * 10 - data_threshold)
        
        if int(np.isnan(v).sum()) > threshold:
            # store the resolved parameters of the fit (amplitude, phase)
            fit_vars[t - it_d,:] = np.nan
            
            wave_0[t - it_d,:] = np.nan
            wave_1[t - it_d,:] = np.nan
            wave_2[t - it_d,:] = np.nan
            
            # do not take non-fits into account by setting to nan
            fit_data[t - it_d] = np.nan 
        
        else:
            inds = np.where(np.isnan(v))
            winds = np.delete(v,inds)
            lons = np.delete(lons,inds)
                    
            popt, pcov = curve_fit(superdarn_SL_SDT_full,
                                   lons,
                                   winds,
                                   maxfev=20000)
                
            # store the resolved parameters of the fit (amplitude, phase)
            fit_vars[t - it_d,:] = popt[:]
            
            wave_0[t - it_d,:] = sdt(lon_fit,popt[0],popt[2])
            wave_1[t - it_d,:] = sdt(lon_fit,popt[1],popt[3])
            
    return fit_vars, wave_0, wave_1, fit_data


def superdarn_SLT_sdt_m2n2_fit(measurements,
                               longitudes,
                               x_dim,
                               window_w,
                               tstep,
                               data_threshold,
                               only_lunar=False,
                               only_solar=False):
    
    """ 
    longitudes must be in radians
    
    """
    
    # number of measurements (iterations)
    iters = int(np.shape(measurements)[0] - window_w)
    it_d = int(window_w / 2)
    
    # array to store the fit variables in 
    fit_vars = np.zeros((iters,6))
    fit_data = np.empty(iters)
    
    # array to store wave fit in
    wave_0 = np.zeros((iters,x_dim))
    wave_1 = np.zeros((iters,x_dim)) 
    wave_2 = np.zeros((iters,x_dim)) 
    
    # x-axis on which fit is projected
    lon_fit = np.linspace(0,2 * np.pi,x_dim,endpoint=False)
    
    t_steps = np.arange(it_d, iters + it_d)[::tstep]
    
    for t in t_steps:
        # select data and filter np.nans from longitude and data array
        v = measurements[(t - it_d):(t + it_d),:].flatten()
        
        # number of data points in fit
        fit_data[t - it_d] = int(window_w * 10 - np.isnan(v).sum())
        
        # fill longitude array with (arbitrary) 4pi spacing between t-steps
        lons = np.empty((window_w,np.shape(longitudes)[0]))
        for k in range(window_w):
            lons[k,:] = longitudes + k * 4 * np.pi

        lons.flatten()
        
        # if more nan values than threshold, reject fit
        threshold = int(window_w * 10 - data_threshold)
        
        if int(np.isnan(v).sum()) > threshold:
            # store the resolved parameters of the fit (amplitude, phase)
            fit_vars[t - it_d,:] = np.nan
            
            wave_0[t - it_d,:] = np.nan
            wave_1[t - it_d,:] = np.nan
            wave_2[t - it_d,:] = np.nan
            
            # do not take non-fits into account by setting to nan
            fit_data[t - it_d] = np.nan 
        
        else:
            inds = np.where(np.isnan(v))
            winds = np.delete(v,inds)
            lons = np.delete(lons,inds)
                    
            popt, pcov = curve_fit(superdarn_SL_SDT_m2n2,
                                   lons,
                                   winds,
                                   maxfev=20000)
                
            # store the resolved parameters of the fit (amplitude, phase)
            fit_vars[t - it_d,:] = popt[:]
            
            wave_0[t - it_d,:] = sdt(lon_fit,popt[0],popt[2])
            wave_1[t - it_d,:] = sdt(lon_fit,popt[1],popt[3])
            
    return fit_vars, wave_0, wave_1, fit_data


def superdarn_SLT_sdt_m2n2_fit_ex(measurements,
                                  longitudes,
                                  x_dim,
                                  window_w,
                                  tstep,
                                  data_threshold,
                                  only_lunar=False,
                                  only_solar=False):
    
    """ 
    longitudes must be in radians
    
    """
    
    # number of measurements (iterations)
    iters = int(np.shape(measurements)[0] - window_w)
    it_d = int(window_w / 2)
    
    # array to store the fit variables in 
    fit_vars = np.zeros((iters,6))
    t_steps = np.arange(it_d, iters + it_d)[::tstep]
    
    for t in t_steps:
        # select data and filter np.nans from longitude and data array
        v = measurements[(t - it_d):(t + it_d),:].flatten()
                
        # fill longitude array with (arbitrary) 4pi spacing between t-steps
        lons = np.empty((window_w,np.shape(longitudes)[0]))
        for k in range(window_w):
            lons[k,:] = longitudes + k * 4 * np.pi

        lons.flatten()
        
        # if more nan values than threshold, reject fit
        threshold = int(window_w * 10 - data_threshold)
        
        if int(np.isnan(v).sum()) > threshold:
            # store the resolved parameters of the fit (amplitude, phase)
            fit_vars[t - it_d,:] = np.nan
        
        else:
            inds = np.where(np.isnan(v))
            winds = np.delete(v,inds)
            lons = np.delete(lons,inds)
                    
            popt, pcov = curve_fit(superdarn_SL_SDT_m2n2,
                                   lons,
                                   winds,
                                   maxfev=20000)
                
            # store the resolved parameters of the fit (amplitude, phase)
            fit_vars[t - it_d,:] = popt[:]

    return fit_vars


def superdarn_SLT_sdt_m2n2_eastwest(measurements,
                                    longitudes,
                                    x_dim,
                                    window_w,
                                    tstep,
                                    data_threshold):
    
    """ 
    longitudes must be in radians
    
    """
    
    # number of measurements (iterations)
    iters = int(np.shape(measurements)[0] - window_w)
    it_d = int(window_w / 2)
    
    # array to store the fit variables in 
    fit_vars = np.zeros((iters,12))
    t_steps = np.arange(it_d, iters + it_d)[::tstep]
    
    for t in t_steps:
        # select data and filter np.nans from longitude and data array
        v = measurements[(t - it_d):(t + it_d),:].flatten()
                
        # fill longitude array with (arbitrary) 4pi spacing between t-steps
        lons = np.empty((window_w,np.shape(longitudes)[0]))
        for k in range(window_w):
            lons[k,:] = longitudes + k * 4 * np.pi

        lons.flatten()
        
        # if more nan values than threshold, reject fit
        threshold = int(window_w * 10 - data_threshold)
        
        if int(np.isnan(v).sum()) > threshold:
            # store the resolved parameters of the fit (amplitude, phase)
            fit_vars[t - it_d,:] = np.nan
        
        else:
            inds = np.where(np.isnan(v))
            winds = np.delete(v,inds)
            lons = np.delete(lons,inds)
                    
            popt, pcov = curve_fit(superdarn_SL_SDT_m2n2_eastwest,
                                   lons,
                                   winds,
                                   maxfev=20000)
                
            # store the resolved parameters of the fit (amplitude, phase)
            fit_vars[t - it_d,:] = popt[:]

    return fit_vars


def superdarn_single_fit(measurements,
                         longitudes,
                         window_w,
                         tstep,
                         data_threshold,
                         period=None):
    
    """ 
    - longitudes must be in radians
    - function accepts period. 1200 = 12.00 hr wave
    
    """
    
    # number of measurements (iterations)
    iters = int(np.shape(measurements)[0] - window_w)
    it_d = int(window_w / 2)
    
    # array to store the fit variables in 
    fit_vars = np.zeros((iters,2))
    fit_data = np.empty(iters)
        
    t_steps = np.arange(it_d, iters + it_d)[::tstep]
    
    def superdarn_single(x,
                         A,
                         phi_A):
        
        """
        nested single wave-form function
        
        """
    
        return A * np.sin(2 * ( x + x // (4 * np.pi) * 2 * np.pi * 50 / period) + phi_A)

    
    for t in t_steps:
        # select data and filter np.nans from longitude and data array
        v = measurements[(t - it_d):(t + it_d),:].flatten()
        
        # number of data points in fit
        fit_data[t - it_d] = int(window_w * 10 - np.isnan(v).sum())
        
        # fill longitude array with (arbitrary) 4pi spacing between t-steps
        lons = np.empty((window_w,np.shape(longitudes)[0]))
        for k in range(window_w):
            lons[k,:] = longitudes + k * 4 * np.pi

        lons.flatten()
        
        # if more nan values than threshold, reject fit
        threshold = int(window_w * 10 - data_threshold)
        
        if int(np.isnan(v).sum()) > threshold:
            # store the resolved parameters of the fit (amplitude, phase)
            fit_vars[t - it_d,:] = np.nan
                        
            # do not take non-fits into account by setting to nan
            fit_data[t - it_d] = np.nan 
        
        else:
            inds = np.where(np.isnan(v))
            winds = np.delete(v,inds)
            lons = np.delete(lons,inds)
                    
            popt, pcov = curve_fit(superdarn_single,
                                   lons,
                                   winds,
                                   maxfev=20000)
                
            # store the resolved parameters of the fit (amplitude, phase)
            fit_vars[t - it_d,:] = popt[:]
            
    return fit_vars, fit_data


def superdarn_single_s1(measurements,
                        longitudes,
                        window_w,
                        tstep,
                        data_threshold,
                        period=None):
    
    """ 
    - longitudes must be in radians
    - function accepts period. 1200 = 12.00 hr wave
    
    """
    
    # number of measurements (iterations)
    iters = int(np.shape(measurements)[0] - window_w)
    it_d = int(window_w / 2)
    
    # array to store the fit variables in 
    fit_vars = np.zeros((iters,2))
    fit_data = np.empty(iters)
        
    t_steps = np.arange(it_d, iters + it_d)[::tstep]
    
    def superdarn_single(x,
                         A,
                         phi_A):
        
        """
        nested single wave-form function
        
        """
    
        return A * np.sin(2 * ( x - x // (4 * np.pi) * 2 * np.pi * 1 / (period * 2)) + phi_A)

    
    for t in t_steps:
        # select data and filter np.nans from longitude and data array
        v = measurements[(t - it_d):(t + it_d),:].flatten()
        
        # number of data points in fit
        fit_data[t - it_d] = int(window_w * 10 - np.isnan(v).sum())
        
        # fill longitude array with (arbitrary) 4pi spacing between t-steps
        lons = np.empty((window_w,np.shape(longitudes)[0]))
        for k in range(window_w):
            lons[k,:] = longitudes + k * 4 * np.pi

        lons.flatten()
        
        # if more nan values than threshold, reject fit
        threshold = int(window_w * 10 - data_threshold)
        
        if int(np.isnan(v).sum()) > threshold:
            # store the resolved parameters of the fit (amplitude, phase)
            fit_vars[t - it_d,:] = np.nan
                        
            # do not take non-fits into account by setting to nan
            fit_data[t - it_d] = np.nan 
        
        else:
            inds = np.where(np.isnan(v))
            winds = np.delete(v,inds)
            lons = np.delete(lons,inds)
                    
            popt, pcov = curve_fit(superdarn_single,
                                   lons,
                                   winds,
                                   maxfev=20000)
                
            # store the resolved parameters of the fit (amplitude, phase)
            fit_vars[t - it_d,:] = popt[:]
            
    return fit_vars, fit_data


def superdarn_stationary_SDT(measurements,
                             longitudes,
                             window_w,
                             tstep,
                             data_threshold,
                             period=None):
    
    """ 
    - longitudes must be in radians
    - function accepts period. 1200 = 12.00 hr wave
    
    """
    
    # number of measurements (iterations)
    iters = int(np.shape(measurements)[0] - window_w)
    it_d = int(window_w / 2)
    
    # array to store the fit variables in 
    fit_vars = np.zeros((iters,2))
    fit_data = np.empty(iters)
        
    t_steps = np.arange(it_d, iters + it_d)[::tstep]
    
#    def superdarn_single(x,
#                         A,
#                         phi_A,
#                         B,
#                         phi_B):
#        
#        """
#        nested single wave-form function
#        
#        """
#    
#        return (  A * np.sin(2 * ( x + x // (4 * np.pi) * 2 * np.pi * 1 / (period * 2)) + phi_A) 
#                + A * np.sin(2 * ( x - x // (4 * np.pi) * 2 * np.pi * 1 / (period * 2)) + phi_A)
#                + B * np.sin(2 * ( x + x // (4 * np.pi) * 2 * np.pi * 1 / (period * 2)) + phi_B))

    def superdarn_single(x,
                         B,
                         phi_B):
        
        """
        nested single wave-form function
        
        """
    
        return (B * np.sin(2 * ( x + x // (4 * np.pi) * 2 * np.pi * 1 / (period * 2)) + phi_B))

    
    for t in t_steps:
        # select data and filter np.nans from longitude and data array
        v = measurements[(t - it_d):(t + it_d),:].flatten()
        
        # number of data points in fit
        fit_data[t - it_d] = int(window_w * 10 - np.isnan(v).sum())
        
        # fill longitude array with (arbitrary) 4pi spacing between t-steps
        lons = np.empty((window_w,np.shape(longitudes)[0]))
        for k in range(window_w):
            lons[k,:] = longitudes + k * 4 * np.pi

        lons.flatten()
        
        # if more nan values than threshold, reject fit
        threshold = int(window_w * 10 - data_threshold)
        
        if int(np.isnan(v).sum()) > threshold:
            # store the resolved parameters of the fit (amplitude, phase)
            fit_vars[t - it_d,:] = np.nan
                        
            # do not take non-fits into account by setting to nan
            fit_data[t - it_d] = np.nan 
        
        else:
            inds = np.where(np.isnan(v))
            winds = np.delete(v,inds)
            lons = np.delete(lons,inds)
                    
            popt, pcov = curve_fit(superdarn_single,
                                   lons,
                                   winds,
                                   maxfev=20000)
                
            # store the resolved parameters of the fit (amplitude, phase)
            fit_vars[t - it_d,:] = popt[:]
            
    return fit_vars, fit_data


def superdarn_double_fit(measurements,
                         longitudes,
                         window_w,
                         tstep,
                         data_threshold,
                         period_1=1200,
                         period_2=1200):
    
    """ 
    - longitudes must be in radians
    - function accepts period. 1200 = 12.00 hr wave
    
    """
    
    # number of measurements (iterations)
    iters = int(np.shape(measurements)[0] - window_w)
    it_d = int(window_w / 2)
    
    # array to store the fit variables in 
    fit_vars = np.zeros((iters,4))
    fit_data = np.empty(iters)
        
    t_steps = np.arange(it_d, iters + it_d)[::tstep]
    
    def superdarn_single(x,
                         A,
                         B,
                         phi_A,
                         phi_B):
        
        """
        nested single double wave-form. [B, phi_B] correspond to exact SD2W
        
        """
    
        return (A * np.sin(2 * ( x + x // (4 * np.pi) * 2 * np.pi * 50 / period_1) + phi_A) +
                B * np.sin(2 * ( x + x // (4 * np.pi) * 2 * np.pi * 50 / period_2) + phi_B))

    
    for t in t_steps:
        # select data and filter np.nans from longitude and data array
        v = measurements[(t - it_d):(t + it_d),:].flatten()
        
        # number of data points in fit
        fit_data[t - it_d] = int(window_w * 10 - np.isnan(v).sum())
        
        # fill longitude array with (arbitrary) 4pi spacing between t-steps
        lons = np.empty((window_w,np.shape(longitudes)[0]))
        for k in range(window_w):
            lons[k,:] = longitudes + k * 4 * np.pi

        lons.flatten()
        
        # if more nan values than threshold, reject fit
        threshold = int(window_w * 10 - data_threshold)
        
        if int(np.isnan(v).sum()) > threshold:
            # store the resolved parameters of the fit (amplitude, phase)
            fit_vars[t - it_d,:] = np.nan
                        
            # do not take non-fits into account by setting to nan
            fit_data[t - it_d] = np.nan 
        
        else:
            inds = np.where(np.isnan(v))
            winds = np.delete(v,inds)
            lons = np.delete(lons,inds)
                    
            popt, pcov = curve_fit(superdarn_single,
                                   lons,
                                   winds,
                                   maxfev=20000)
                
            # store the resolved parameters of the fit (amplitude, phase)
            fit_vars[t - it_d,:] = popt[:]
            
    return fit_vars, fit_data


def superdarn_tripple_fit(measurements,
                          longitudes,
                          window_w,
                          tstep,
                          data_threshold,
                          period_1=1200,
                          period_2=1200,
                          period_3=1200):
    
    """ 
    - longitudes must be in radians
    - function accepts period. 1200 = 12.00 hr wave
    
    """
    
    # number of measurements (iterations)
    iters = int(np.shape(measurements)[0] - window_w)
    it_d = int(window_w / 2)
    
    # array to store the fit variables in 
    fit_vars = np.zeros((iters,6))
    fit_data = np.empty(iters)
        
    t_steps = np.arange(it_d, iters + it_d)[::tstep]
    
    def superdarn_single(x,
                         A,
                         B,
                         C,
                         phi_A,
                         phi_B,
                         phi_C):
        
        """
        nested single double wave-form. [B, phi_B] correspond to exact SD2W
        
        """
    
        return (A * np.sin(2 * ( x + x // (4 * np.pi) * 2 * np.pi * 50 / period_1) + phi_A) +
                B * np.sin(2 * ( x + x // (4 * np.pi) * 2 * np.pi * 50 / period_2) + phi_B) + 
                C * np.sin(2 * ( x + x // (4 * np.pi) * 2 * np.pi * 50 / period_3) + phi_C))

    
    for t in t_steps:
        # select data and filter np.nans from longitude and data array
        v = measurements[(t - it_d):(t + it_d),:].flatten()
        
        # number of data points in fit
        fit_data[t - it_d] = int(window_w * 10 - np.isnan(v).sum())
        
        # fill longitude array with (arbitrary) 4pi spacing between t-steps
        lons = np.empty((window_w,np.shape(longitudes)[0]))
        for k in range(window_w):
            lons[k,:] = longitudes + k * 4 * np.pi

        lons.flatten()
        
        # if more nan values than threshold, reject fit
        threshold = int(window_w * 10 - data_threshold)
        
        if int(np.isnan(v).sum()) > threshold:
            # store the resolved parameters of the fit (amplitude, phase)
            fit_vars[t - it_d,:] = np.nan
                        
            # do not take non-fits into account by setting to nan
            fit_data[t - it_d] = np.nan 
        
        else:
            inds = np.where(np.isnan(v))
            winds = np.delete(v,inds)
            lons = np.delete(lons,inds)
                    
            popt, pcov = curve_fit(superdarn_single,
                                   lons,
                                   winds,
                                   maxfev=20000)
                
            # store the resolved parameters of the fit (amplitude, phase)
            fit_vars[t - it_d,:] = popt[:]
            
    return fit_vars, fit_data


def running_mean(x, 
                 N):
    
    """
    calculates running mean via convolution, does not handle nan values
    
    """
    
    array = np.copy(x) 
    
    for t in range(int(N / 2), np.shape(array)[0] - int(N / 2)):
        array[t] = np.nanmean(x[int(t - N / 2):int(t + N / 2)])
        
    return array


def running_mean_nans(x,
                      N):
    
    """
    running mean for arrays containing nan values
    
    """

    array = np.copy(x) 
    
    for t in range(int(N / 2), np.shape(array)[0] - int(N / 2)):
        array[t] = np.nanmean(x[int(t - N / 2):int(t + N / 2)])
        
    return array
    

def monthly_means(input_array,
                  n_yrs):
    
    mms = np.zeros(n_yrs * 12)
    month_indices = [30,31,30,31,30,31,31,30,30,31,30,30]
    
    f = 0 
    for k in range(n_yrs):
        m = 0
        for l in month_indices:
            mms[f] = np.nanmean(input_array[(k * 365 * 24 + m * 24):(k * 365 * 24 + (m + l) * 24)])       
            m += l
            f += 1
            
    return mms


def SD_remove_leap_days_95_16(input_array):
    
    """
    29th of February removed from 1995-2016 dataset
    
    """
    
    day_inds = [423,1884,3345,4806,6268]
    
    inds = []
    
    for day in day_inds:
        a = np.arange(day * 24, (day + 1) * 24) 
        inds = np.concatenate((inds,a))
        
    skimmed_array = np.delete(input_array,(inds),axis=0)
        
    return skimmed_array

def remove_leap_days_95_16(input_array):
    
    """
    29th of February removed from 1995-2016 dataset
    
    """
    
    day_inds = [423,1884,3345,4806,6268]
    
    inds = []
    
    for day in day_inds:
        a = np.arange(day * 24, (day + 1) * 24) 
        inds = np.concatenate((inds,a),axis=0)
    
    inds = np.array(inds,dtype='int')
    skimmed_array = np.delete(input_array,inds,axis=0)
        
    return skimmed_array


def single_station_tide(data,
                        window_w,
                        threshold,
                        tide=None):
    
    """
    -- tide indices: 0. DT, 1. SDT, 2. TDT
    -- window_w in number of hours (e.g. 10 * 24 for 10 days)
    
    """
    
    tides = ['dt','sdt','tdt'][tide]
    
    # set meta-data for fit routine
    n_days = int(window_w / 24)
    xdim = np.shape(data)[1]
    rads = np.linspace(0, 2 * np.pi * n_days, window_w, endpoint=False)
    len_fit = int(np.shape(data)[0] - window_w)
    
    # number of nans needed to exceed set threshold
    thresh = int(window_w - threshold) 
    
    # arrays to project data onto
    amps = np.zeros((len_fit,xdim))
    phis = np.copy(amps)
    fit_data = np.copy(amps)
    
    for k in range(xdim):
        for j in range(len_fit):
            temp_dat = data[j:(j + window_w),k]
            fit_data[j,k] = int(window_w - np.isnan(temp_dat).sum())
                      
            if int(np.isnan(temp_dat).sum()) > thresh:
                 amps[j,k],phis[j,k],fit_data[j,k] = (np.nan,np.nan,np.nan)
            else:
                # remove nans
                inds = np.where(np.isnan(temp_dat))
                temp_dat = np.delete(temp_dat,inds)
                radi = np.delete(rads,inds)
                
                if tides == 'dt':
                    amps[j,k],phis[j,k] = dt_fitter(temp_dat,radi)
                elif tides == 'sdt':
                    amps[j,k],phis[j,k] = sdt_fitter(temp_dat,radi)
                elif tides == 'tdt':
                    amps[j,k],phis[j,k] = tdt_fitter(temp_dat,radi)
                    
    return amps, phis, fit_data


def mae(stations,
        surface_fit,
        threshold):
    
    """
    -- calculate mean absolute error at each time step
    -- threshold sets minimum number of data points per MAE calculation
    
    checked 6-06-2019
    
    """
    
    t_steps = np.shape(surface_fit)[0]
    mae = np.empty(t_steps)
    
    for t in range(t_steps):
        migrating_A = surface_fit[t]
        local_A = stations[t,:]
        
        data_p = np.isnan(local_A).sum()
        n = 10 - data_p

        if migrating_A == np.nan:
            mae[t] = np.nan
        elif data_p > (10 - threshold):
            mae[t] = np.nan
        else:
            inds = np.where(np.isnan(local_A))
            local_A = np.delete(local_A,inds)
            
            # calculate Mean Absolute Error
            mae[t] = np.sum(np.abs(local_A - migrating_A)) / n

    return mae


def rmse(stations,
         surface_fit,
         threshold):
    """
    """
    
    t_steps = np.shape(surface_fit)[0]
    rmse = np.empty(t_steps)
    
    for t in range(t_steps):
        migrating_A = surface_fit[t]
        local_A = stations[t,:]
        
        data_p = np.isnan(local_A).sum()
        n = 10 - data_p

        if migrating_A == np.nan:
            rmse[t] = np.nan
        elif data_p > (10 - threshold):
            rmse[t] = np.nan
        else:
            inds = np.where(np.isnan(local_A))
            local_A = np.delete(local_A,inds)
            
            rmse[t] = np.sqrt(np.sum((local_A - migrating_A) ** 2) / n)

    return rmse   
    

def dimdata_to_power_of2(data,
                         N,
                         M):
    """
    Map data onto 2^N * 2^M grid where empty spaces are filled with zeros. This
    format of the data goes well with the stockwell transform algorithm.
    
    Standard N, M = 6, 13. Corresponds to 64 and 8192 grid points
    """  
    
    grid = np.zeros((2**M, 2**N))
    
    dim_dat = np.shape(data)[0]
    grid[:dim_dat,:] = data
    
    return grid


def east_west_2DFFT_split(data,
                          wave_n):
    
    """
    Input data with dimensions being of the form 2^k
    """

    y_center = int(np.shape(data)[0] / 2)
    end_dim = (np.shape(data)[1] - wave_n)
    
    transform = np.fft.fft2(data)
    
    transform_east = np.zeros(np.shape(transform),dtype=np.complex_)
    transform_west = np.zeros(np.shape(transform),dtype=np.complex_)
    
    transform_east[y_center:,wave_n] = transform[y_center:, wave_n]
    transform_east[:y_center,end_dim] = transform[:y_center, end_dim]
    
    transform_west[:y_center,wave_n] = transform[:y_center, wave_n]
    transform_west[y_center:,end_dim] = transform[y_center:, end_dim]
    
    eastward_field = np.fft.ifft2(transform_east)
    westward_field = np.fft.ifft2(transform_west)
    
    eastward_field = np.real(eastward_field)
    westward_field = np.real(westward_field)
    
    return eastward_field, westward_field


def superdarn_SDT_climatology(window,
                              var):
    
    """
    Function input: - Window: size in (integer) number of days
                    - Variable: zonal (U) or meridional (V) wind
                    
    Output: - Yearly climatology of the SDT amplitude
            - Standard deviation of the year to year variability in amplitude
            - Yearly climatology of 24-hr subsampled SDT phase
            - Standard deviation of the year to year variability in phase
    """
        
    ww = 'ww_' + str(int(window / 24)) + 'd_'
    
    save_f = '/home/wim/Desktop/Projects/SLT-beat/Fits/'
    varibs = np.load(save_f + 'SD_' + var + '_variables_1995_2016_' + ww + '.npy')
    
    # remove leap days to make each year fit in 365 days
    varibs = SD_remove_leap_days_95_16(varibs)
    
    # number of data points per year
    one_yr = 365 * 24
    
    climatology_sdt = np.empty((21,one_yr))
    phase_climatology_sdt = np.empty((21,one_yr))
    
    # adjust phase by 180 degrees when fitted amplitude is negative
    inds_sdt = np.where(varibs[:,1] < 0)
    varibs[:,4][inds_sdt] += np.pi 
      
    # after adjusting phase, absolute values of amplitudes can be used
    for k in range(21):
        climatology_sdt[k,:] = np.abs(varibs[(k * one_yr):((k + 1) * one_yr),1]) 
        phase_climatology_sdt[k,:] = varibs[(k * one_yr):((k + 1) * one_yr),4]
    
    # amplitude climatology
    sdt_amp_climatology = np.nanmean(climatology_sdt,axis=0)
    sdt_amp_std = np.nanstd(climatology_sdt,axis=0)
    
    # deal with nans and calculate climatological phase means and std
    phases_sdt = np.empty(one_yr)
    stds_sdt = np.empty(one_yr)
    
    for m in range(one_yr):   
        temp_sdt = phase_climatology_sdt[:,m]
        inds_sdt = np.where(np.isnan(temp_sdt))
        temp_sdt = np.delete(temp_sdt,inds_sdt)
        
        phases_sdt[m] = circmean(temp_sdt)
        stds_sdt[m] = circstd(temp_sdt)
        
    # sample and center phase progression around 6 hrs
    sample_r = 24
    center = np.pi
    
    phase_cen = phases_sdt[::sample_r] - circmean(phases_sdt[::sample_r]) + center
    sdt_phase_climatology = circmean(np.stack((phase_cen,phase_cen),axis=1),axis=1)
    sdt_phase_std = stds_sdt[::sample_r]
    
    return sdt_amp_climatology, sdt_amp_std, sdt_phase_climatology, sdt_phase_std


def navgem_SDT(window,
               var,
               md):
    """
    """

    md = 12              # number of days in running mean
    sample_r = 8         # phases sampled every sample_r hours
    window = int(6 * 8)  # window width of surface fit

    ww = 'ww_' + str(window / 8) + 'd'
    save_f = '/home/wim/Desktop/Projects/Tides SuperDARN/NAVGEM Experiments/Migrating tides fit/'
    vars_sd = np.load(save_f + 'NAVGEM_' + "SD_latlon" + '_78.750km_' + var + '_variables_' + ww + '.npy')
    vars_na = np.load(save_f + 'NAVGEM_' + "NAVGEM_full" + '_78.750km_' + var + '_variables_' + ww + '.npy')

    # correct phase by 180 degrees when fitted amplitude is negative
    inds_sdt_sd = np.where(vars_sd[:,1] < 0)
    vars_sd[:,4][inds_sdt_sd] += np.pi

    inds_sdt_na = np.where(vars_na[:,1] < 0)
    vars_na[:,4][inds_sdt_na] += np.pi

    center = np.pi
    phase_sd = vars_sd[::sample_r,4] - circmean(vars_sd[::sample_r,4]) + center
    phase_na = vars_na[::sample_r,4] - circmean(vars_na[::sample_r,4]) + center
    
    phase_sd = circmean(np.stack((phase_sd,phase_sd),axis=1),axis=1)
    phase_na = circmean(np.stack((phase_na,phase_na),axis=1),axis=1)
    
    amp_sd = running_mean(np.abs(vars_sd[:,1]),md)
    amp_na = running_mean(np.abs(vars_na[:,1]),md)
    
    return amp_sd, phase_sd, amp_na, phase_na


def navgem_SDT_band_fitter(window,
                           var,
                           md,
                           height_i,
                           lat_i):
    
    sample_r = 8         # phases sampled every sample_r hours
    
    data_sd, longitudes_sd = interpolated_field_SDlat(height_i,var)
    varibs, dt, sdt, tdt = navgem_dt_sdt_tdt_fitter(data_sd,
                                                    longitudes_sd,
                                                    32,
                                                    window,
                                                    1)
    
    vars_sd = np.copy(varibs)
    
    nan_map = NAVGEM_full_nan_map()
    data, longitudes = interpolated_field_fixlat(height_i,lat_i,nan_map,var)
            
    varibs, dt, sdt, tdt = navgem_dt_sdt_tdt_fitter(data,
                                                    longitudes,
                                                    32,
                                                    window,
                                                    1)
    
    vars_na = np.copy(varibs)

    # correct phase by 180 degrees when fitted amplitude is negative
    inds_sdt_sd = np.where(vars_sd[:,1] < 0)
    vars_sd[:,4][inds_sdt_sd] += np.pi

    inds_sdt_na = np.where(vars_na[:,1] < 0)
    vars_na[:,4][inds_sdt_na] += np.pi

    center = np.pi
    phase_sd = vars_sd[::sample_r,4] - circmean(vars_sd[::sample_r,4]) + center
    phase_na = vars_na[::sample_r,4] - circmean(vars_na[::sample_r,4]) + center
    
    phase_sd = circmean(np.stack((phase_sd,phase_sd),axis=1),axis=1)
    phase_na = circmean(np.stack((phase_na,phase_na),axis=1),axis=1)
    
    amp_sd = running_mean(np.abs(vars_sd[:,1]),md)
    amp_na = running_mean(np.abs(vars_na[:,1]),md)
    
    return amp_sd, phase_sd, amp_na, phase_na


def superdarn_2014(window,
                   var,
                   year):
    
    """
    Function input: - Window: size in (integer) number of days
                    - Variable: zonal (U) or meridional (V) wind
                    - year (index, 19 = 2014)
                    
    Output: - 2014 SDT amps SuperDARN
            - 2014 SDT phases centered over 6 hrs

    """
        
    ww = 'ww_' + str(int(window / 24)) + 'd_'
    
    save_f = '/home/wim/Desktop/Projects/SLT-beat/Fits/'
    varibs = np.load(save_f + 'SD_' + var + '_variables_1995_2016_' + ww + '.npy')
    
    # remove leap days to make each year fit in 365 days
    varibs = SD_remove_leap_days_95_16(varibs)
    
    # number of data points per year and year index (19 = 2014)
    one_yr = 365 * 24
    yr = year
    
    # adjust phase by 180 degrees when fitted amplitude is negative
    inds_sdt = np.where(varibs[:,1] < 0)
    varibs[:,4][inds_sdt] += np.pi 
      
    sdt_amp = np.abs(varibs[(yr * one_yr):((yr + 1) * one_yr):,1])
    phis = varibs[(yr * one_yr):((yr + 1) * one_yr):,4]
      
    # sample and center phase progression around 6 hrs
    sample_r = 24
    center = np.pi
    
    phase_cen = phis[::sample_r] - circmean(phis[::sample_r]) + center
    sdt_phase = circmean(np.stack((phase_cen,phase_cen),axis=1),axis=1)

    return sdt_amp, sdt_phase


def psd_data_sim(years=None,
                 base_period=None,
                 jump_1_month=None,
                 jump_2_month=None,
                 jump_degree=None,
                 jump_spread=None,
                 noise_std=None,
                 A_365=None,
                 A_180=None,
                 A_120=None,
                 A_90=None,
                 A_70=None,
                 A_0=None):
    
    """
    - jump_degree in radians, pi = 180 degrees. 
    - modulation periods corresponding to peak powers in climatology.
    
    """
    
    x = np.arange(years * 24 * 365,dtype='float')
    inds_base = np.arange(24 * 30 * jump_1_month, 24 * 30 * jump_2_month)
    
    for k in range(years):
        inds = inds_base + k * 24 * 365 + random.randint(-jump_spread / 2, jump_spread / 2)
        x[inds] = x[inds] + int(jump_degree / np.pi * base_period / 2)         
    
    base = np.sin(2 * np.pi / (base_period) * x)
    mod_1 = A_365 * np.sin(2 * np.pi / (365 * 24) * x)
    mod_2 = A_180 * np.sin(2 * np.pi / (180 * 24) * x)
    mod_3 = A_120 * np.sin(2 * np.pi / (120 * 24) * x)
    mod_4 = A_90 * np.sin(2 * np.pi / (90 * 24) * x)
    mod_5 = A_70 * np.sin(2 * np.pi / (60 * 24) * x)
    
    y = base * (mod_1 + mod_2 + mod_3 + mod_4 + mod_5 + A_0) + np.random.normal(0,noise_std,size=np.shape(x)[0])
    t = np.arange(years * 24 * 365,dtype='float')
    
    return y, t


def model_SLT_fits(measurements,
                   longitudes,
                   window_w,
                   tstep,
                   func_ind):
    
    """ 
    - longitudes must be in radians
    - assumes 3-hourly input data
    - func_ind (0,1,2) = (S2M2,S2,M2) surface fits
    - function checked 23-10-2019
    
    """
    
    # number of measurements (iterations)
    iters = int(np.shape(measurements)[0] - window_w)
    it_d = int(window_w / 2)
    
    # arrays to store the fit variables in
    if func_ind == 0:
        fit_vars = np.zeros((iters,4))
    else:
        fit_vars = np.zeros((iters,2))
        
    # subsample full time array to save time
    t_steps = np.arange(it_d, iters + it_d)[::tstep]
        
    def S2M2(x,
             A,
             B,
             phi_A,
             phi_B):
    
        return (  A * np.sin(2 * ( x + x // (4 * np.pi) * 2 * np.pi * 50 / 1200 * 3) + phi_A)
                + B * np.sin(2 * ( x + x // (4 * np.pi) * 2 * np.pi * 50 / 1242 * 3) + phi_B))
        
    def S2(x,
           A,
           phi_A):
    
        return (  A * np.sin(2 * ( x + x // (4 * np.pi) * 2 * np.pi * 50 / 1200 * 3) + phi_A))
    
    def M2(x,
           A,
           phi_A):
    
        return (  A * np.sin(2 * ( x + x // (4 * np.pi) * 2 * np.pi * 50 / 1242 * 3) + phi_A))
    
    fit_functions = [S2M2,S2,M2]
    
    for t in t_steps:
        # select data and filter np.nans from longitude and data array
        v = measurements[(t - it_d):(t + it_d),:].flatten()
                
        # fill longitude array with (arbitrarily) 4pi spacing between t-steps
        lons = np.empty((window_w,np.shape(longitudes)[0]))
        for k in range(window_w):
            lons[k,:] = longitudes + k * 4 * np.pi

        lons = lons.flatten()
        
        popt, pcov = curve_fit(fit_functions[func_ind],
                               lons,
                               v,
                               maxfev=20000)
            
        # store the resolved parameters of the fit (amplitude, phase)
        fit_vars[t - it_d,:] = popt[:]
                        
    return fit_vars


def geometric_to_geopotential(z):
    """
    Equation referenced in 'geopot_to_geometric.pdf' in /Support files.
    
    Function checked 16-09-2019 (dd/mm/year).
    
    """
    
    r_earth = 6371008 # mean earth radius (meters)
    
    # rounded to two decimal places
    return np.round(z * r_earth / (r_earth + z),2)


def geopotential_to_geometric(z):
    """
    Equation referenced in 'geopot_to_geometric.pdf' in /Support files.
    
    Function checked 16-09-2019 (dd/mm/year).
    
    """
    
    r_earth = 6371008 # mean earth radius (meters)
    
    # rounded to two decimal places
    return np.round(z * r_earth / (r_earth - z),2)

def trondheim_timeseries(alt_ind,
                         years_list):
    """
    u and v from hourly Trondheim meteor radar from 2013 to 2018.
    
    alt_ind: [0,1,2,3,4,5] = [82,85,88,91,94,97] km
    
    Possible years: - ['2012','2013','2014','2015','2016','2017','2018']
    
    Function checked: .....
    
    """
    
    path = '/home/wim/Desktop/NTNU/Trondheim radar/trondheim_meteor_radar_'
    months = ['1','2','3','4','5','6','7','8','9','10','11','12']

    for index, item in enumerate(years_list):
        for month in months:
            if index == 0 and month == '1':
                data = Dataset(path + item + '_' + month + '.nc', mode='r')
                u = data['U'][:]
                v = data['V'][:]
                
                altitude = data['altitude'][:]
                print('Trondheim meteor data at: ', altitude[alt_ind], ' km')
                ind_nan_u = np.where(u == -9999)
                ind_nan_v = np.where(v == -9999)
                u[ind_nan_u] = np.nan
                v[ind_nan_v] = np.nan
            else:
                data = Dataset(path + item + '_' + month + '.nc', mode='r')
                u_temp = data['U'][:]
                v_temp = data['V'][:]
                
                ind_nan_u = np.where(u_temp == -9999)
                ind_nan_v = np.where(v_temp == -9999)
                u_temp[ind_nan_u] = np.nan
                v_temp[ind_nan_v] = np.nan
                
                u = np.concatenate((u,u_temp))
                v = np.concatenate((v,v_temp))
                
    return u[:,alt_ind], v[:,alt_ind]


def trondheim_tide_fit(data,
                       window_w,
                       threshold):
    
    """
    Data: Hourly Trondheim (all single station) wind data.
    window_w: Width of the sliding window (hrs).
    threshold: Max number of missing data points per fit.
    
    Output: Fit parameters for diurnal, semidiurnal, terdiurnal tides and mean.
    
    can be used to calculate 'climatology' (i.e. corresponding output length).
    
    Function checked: .....
    
    """
    
    n_days = window_w / 24
    n_steps = int(np.shape(data)[0])
    x = np.linspace(0,2 * np.pi * n_days,window_w,endpoint=False)
    fit_vars = np.zeros((n_steps,7)) * np.nan
    
    def fit_function(x,
                     A_dt,
                     A_sdt,
                     A_tdt,
                     phi_dt,
                     phi_sdt,
                     phi_tdt,
                     mean):
        
        return (A_dt * np.sin(x + phi_dt) + A_sdt * np.sin(2 * x + phi_sdt)
               + A_tdt * np.sin(3 * x + phi_tdt) + mean)
        
    for t in range(int(window_w / 2), int(n_steps - window_w / 2)):
        temp_data = data[int(t - window_w / 2):int(t + window_w / 2)]
        
        if int(np.isnan(temp_data).sum()) > threshold:
            # nan if not enough data points are present
            fit_vars[t,:] = np.nan
                    
        else:
            inds = np.where(np.isnan(temp_data))
            winds = np.delete(temp_data,inds)
            lons = np.delete(x,inds)
                    
            popt, pcov = curve_fit(fit_function,
                                   lons,
                                   winds,
                                   maxfev=20000)
                
            # store the resolved parameters of the fit (amplitude, phase)
            fit_vars[t,:] = popt[:]
    
    return fit_vars

def trondheim_m2s2_fit(data,
                       window_w,
                       threshold):
    
    """
    Data: Hourly Trondheim (all single station) wind data.
    window_w: Width of the sliding window (hrs).
    threshold: Max number of missing data points per fit.
    
    Output: Fit parameters for diurnal, semidiurnal, terdiurnal tides and mean.
    
    can be used to calculate 'climatology' (i.e. corresponding output length).
    
    Function checked: .....
    
    """
    
    n_days = window_w / 24
    n_steps = int(np.shape(data)[0])
    x = np.linspace(0,2 * np.pi * n_days,window_w,endpoint=False)
    fit_vars = np.zeros((n_steps,9)) * np.nan
    
    def fit_function(x,
                     A_s2,
                     A_m2,
                     A_dt,
                     A_td,
                     phi_s2,
                     phi_m2,
                     phi_dt,
                     phi_td,
                     mean):
        
        return (  A_s2 * np.sin((24 / 12.00) * x + phi_s2) 
                + A_m2 * np.sin((24 / 12.42) * x + phi_m2)
                + A_dt * np.sin((24 / 24.00) * x + phi_dt)
                + A_td * np.sin((24 / 8.000) * x + phi_td)
                + mean)
        
    for t in range(int(window_w / 2), int(n_steps - window_w / 2)):
        temp_data = data[int(t - window_w / 2):int(t + window_w / 2)]
        
        if int(np.isnan(temp_data).sum()) > threshold:
            # nan if not enough data points are present
            fit_vars[t,:] = np.nan
                    
        else:
            inds = np.where(np.isnan(temp_data))
            winds = np.delete(temp_data,inds)
            lons = np.delete(x,inds)
                    
            popt, pcov = curve_fit(fit_function,
                                   lons,
                                   winds,
                                   maxfev=20000)
                
            # store the resolved parameters of the fit (amplitude, phase)
            fit_vars[t,:] = popt[:]
    
    return fit_vars


def trondheim_m2s2_fit2(data,
                        window_w,
                        threshold,
                        lunar_phase):
    
    """
    Data: Hourly Trondheim (all single station) wind data.
    window_w: Width of the sliding window (hrs).
    threshold: Max number of missing data points per fit.
    
    Output: Fit parameters for diurnal, semidiurnal, terdiurnal tides and mean.
    
    can be used to calculate 'climatology' (i.e. corresponding output length).
    
    Function checked: .....
    
    """
    
    n_days = window_w / 24
    n_steps = int(np.shape(data)[0])
    x = np.linspace(0,2 * np.pi * n_days,window_w,endpoint=False)
    fit_vars = np.zeros((n_steps,5)) * np.nan
    
    def fit_function(x,
                     A_dt,
                     A_sdt,
                     phi_dt,
                     phi_sdt,
                     mean):
        
        return (A_dt * np.sin((24 / 12) * x + phi_dt) + A_sdt * np.sin(((24 / 12.42)) * x + phi_sdt)
               + mean)
        
    for t in range(int(window_w / 2), int(n_steps - window_w / 2)):
        temp_data = data[int(t - window_w / 2):int(t + window_w / 2)]
        
        if int(np.isnan(temp_data).sum()) > threshold:
            # nan if not enough data points are present
            fit_vars[t,:] = np.nan
                    
        else:
            inds = np.where(np.isnan(temp_data))
            winds = np.delete(temp_data,inds)
            lons = np.delete(x,inds)
                    
            popt, pcov = curve_fit(fit_function,
                                   lons,
                                   winds,
                                   maxfev=20000)
                
            # store the resolved parameters of the fit (amplitude, phase)
            fit_vars[t,:] = popt[:]
    
    return fit_vars

def trondheim_m2_fit(data,
                     window_w,
                     threshold):
    
    """
    Data: Hourly Trondheim (all single station) wind data.
    window_w: Width of the sliding window (hrs).
    threshold: Max number of missing data points per fit.
    
    Output: Fit parameters for diurnal, semidiurnal, terdiurnal tides and mean.
    
    can be used to calculate 'climatology' (i.e. corresponding output length).
    
    Function checked: .....
    
    """
    
    n_days = window_w / 24
    n_steps = int(np.shape(data)[0])
    x = np.linspace(0,2 * np.pi * n_days,window_w,endpoint=False)
    fit_vars = np.zeros((n_steps,3)) * np.nan
    
    def fit_function(x,
                     A_sdt,
                     phi_sdt,
                     mean):
        
        return (A_sdt * np.sin(((24 / 12.42)) * x + phi_sdt)
               + mean)
        
    for t in range(int(window_w / 2), int(n_steps - window_w / 2)):
        temp_data = data[int(t - window_w / 2):int(t + window_w / 2)]
        
        if int(np.isnan(temp_data).sum()) > threshold:
            # nan if not enough data points are present
            fit_vars[t,:] = np.nan
                    
        else:
            inds = np.where(np.isnan(temp_data))
            winds = np.delete(temp_data,inds)
            lons = np.delete(x,inds)
                    
            popt, pcov = curve_fit(fit_function,
                                   lons,
                                   winds,
                                   maxfev=20000)
                
            # store the resolved parameters of the fit (amplitude, phase)
            fit_vars[t,:] = popt[:]
    
    return fit_vars


def gaussian_weighted_average_2(altitude,
                              parameter,
                              mean,
                              FWHM):
    
    """
    Returns vertical Gaussian weighted average of input parameter.
    
    - altitude input: Can be an array of altitudes (geometric, meters)
    - parameter input: Array of same shape as altitude
    - mean: center of Gaussian (meters)
    - FWHM: full width half-max of the Gaussian (meters)
    
    Function checked 17-09-2019 (dd/mm/year).
    
    """
    
    # translate FWHM to equivalent standard deviation
    std = FWHM / 2.355
       
    # calculate weighted mean
    weights = scipy.stats.norm(mean,std).pdf(altitude)
    mean_param = np.sum(weights * parameter) / np.sum(weights)
    
    return mean_param


def gaussian_weighted_average(altitude,
                              parameter,
                              mean,
                              FWHM):
    
    """
    Returns vertical Gaussian weighted average of input parameter.
    
    - altitude input: Can be an array of altitudes (geometric, meters)
    - parameter input: Array of same shape as altitude
    - mean: center of Gaussian (meters)
    - FWHM: full width half-max of the Gaussian (meters)
    
    Function checked 17-09-2019 (dd/mm/year).
    
    """
    
    # translate FWHM to equivalent standard deviation
    std = FWHM / 2.355
    
    # calculate weighted mean
    weights = scipy.stats.norm(mean,std).pdf(altitude)
    mean_param = np.sum(weights * parameter) / np.sum(weights)
    
    return mean_param, weights


# def superdarn_dt_sdt_tdt_mean_fitter_std(measurements,
#                                          stds,
#                                          longitudes,
#                                          x_dim,
#                                          window_w,
#                                          tstep,
#                                          data_threshold,
#                                          err):
    
#     """ 
#     longitudes in radians
    
#     """
    
#     # number of measurements (iterations)
#     iters = int(np.shape(measurements)[0] - window_w)
#     it_d = int(window_w / 2)
    
#     # array to store the fit variables in 
#     fit_vars = np.zeros((iters,6))
#     errors = np.zeros((iters,6))
#     fit_data = np.empty(iters)
    
#     # array to store wave fit in
#     wave_0 = np.zeros((iters,x_dim))
#     wave_1 = np.zeros((iters,x_dim)) 
#     wave_2 = np.zeros((iters,x_dim)) 
    
#     # x-axis on which fit is projected
#     lon_fit = np.linspace(0,2 * np.pi,x_dim,endpoint=False)
    
#     t_steps = np.arange(it_d, iters + it_d)[::tstep]
    
#     for t in t_steps:
#         # select data and filter np.nans from longitude and data array
#         v = measurements[(t - it_d):(t + it_d),:].flatten()
        
#         stand_err = stds[(t - it_d):(t + it_d),:].flatten()
        
#         # number of data points in fit
#         fit_data[t - it_d] = int(window_w * 10 - np.isnan(v).sum())
        
#         # fill longitude array with (arbitrary) 4pi spacing between t-steps
#         lons = np.empty((window_w,np.shape(longitudes)[0]))
#         for k in range(window_w):
#             lons[k,:] = longitudes + k * 4 * np.pi

#         lons.flatten()
        
#         # if more nan values than threshold, reject fit
#         threshold = int(window_w * 10 - data_threshold)
        
#         if int(np.isnan(v).sum()) > threshold:
#             # store the resolved parameters of the fit (amplitude, phase)
#             fit_vars[t - it_d,:] = np.nan
            
#             wave_0[t - it_d,:] = np.nan
#             wave_1[t - it_d,:] = np.nan
#             wave_2[t - it_d,:] = np.nan
            
#             # do not take non-fits into account by setting to nan
#             fit_data[t - it_d] = np.nan 
        
#         else:
#             inds = np.where(np.isnan(v))
#             winds = np.delete(v,inds)
#             lons = np.delete(lons,inds)
            
#             std_errors = np.delete(stand_err,inds)
            
#             # fix missing STDs by filling with above-average 'error'; this does nothing at the moment
            
# #            if np.isnan(std_errors).sum() > 0:
# #                print('you fool!!!', np.isnan(std_errors).sum())
                
#             # this makes no difference...
#             err_inds = np.where(np.isnan(std_errors))
#             std_errors[err_inds] = err
                                
#             # perform lsqt fit
#             popt, pcov = curve_fit(superdarn_dt_sdt_tdt_mean,
#                                    lons, 
#                                    winds,
#                                    sigma=std_errors,  
#                                    absolute_sigma=True,
#                                    maxfev=20000)
                       
#             # store the resolved parameters of the fit (amplitude, phase)
#             fit_vars[t - it_d,:] = popt[:6]
#             errors[t - it_d,:] = np.sqrt(np.diag(pcov))[:6]
                        
#             wave_0[t - it_d,:] = dt(lon_fit,popt[0],popt[3])
#             wave_1[t - it_d,:] = sdt(lon_fit,popt[1],popt[4])
#             wave_2[t - it_d,:] = tdt(lon_fit,popt[2],popt[5])
            
#     return fit_vars, wave_0, wave_1, wave_2, errors


def superdarn_dt_sdt_tdt_mean_fitter_std(measurements,
                                         stds,
                                         longitudes,
                                         x_dim,
                                         window_w,
                                         tstep,
                                         data_threshold,
                                         err):
    
    """ 
    longitudes in radians
    
    """
    
    # number of measurements (iterations)
    iters = int(np.shape(measurements)[0] - window_w)
    it_d = int(window_w / 2)
    
    # array to store the fit variables in 
    fit_vars = np.zeros((iters,6))
    errors = np.zeros((iters,6))
    fit_data = np.empty(iters)
    
    # array to store wave fit in
    wave_0 = np.zeros((iters,x_dim))
    wave_1 = np.zeros((iters,x_dim)) 
    wave_2 = np.zeros((iters,x_dim)) 
    
    # x-axis on which fit is projected
    lon_fit = np.linspace(0,2 * np.pi,x_dim,endpoint=False)
    
    t_steps = np.arange(it_d, iters + it_d)[::tstep]
    
    for t in t_steps:
        # select data and filter np.nans from longitude and data array
        v = measurements[(t - it_d):(t + it_d),:].flatten()
        
        stand_err = stds[(t - it_d):(t + it_d),:].flatten()
        
        # number of data points in fit
        fit_data[t - it_d] = int(window_w * 10 - np.isnan(v).sum())
        
        # fill longitude array with (arbitrary) 4pi spacing between t-steps
        lons = np.empty((window_w,np.shape(longitudes)[0]))
        for k in range(window_w):
            lons[k,:] = longitudes + k * 4 * np.pi

        lons.flatten()
        
        # if more nan values than threshold, reject fit
        threshold = int(window_w * 10 - data_threshold)
        
        if int(np.isnan(v).sum()) > threshold:
            # store the resolved parameters of the fit (amplitude, phase)
            fit_vars[t - it_d,:] = np.nan
            
            wave_0[t - it_d,:] = np.nan
            wave_1[t - it_d,:] = np.nan
            wave_2[t - it_d,:] = np.nan
            
            # do not take non-fits into account by setting to nan
            fit_data[t - it_d] = np.nan 
        
        else:
            inds = np.where(np.isnan(v))
            winds = np.delete(v,inds)
            lons = np.delete(lons,inds)
            
            std_errors = np.delete(stand_err,inds)
            
            # fix missing STDs by filling with above-average 'error'; this does nothing at the moment
            
#            if np.isnan(std_errors).sum() > 0:
#                print('you fool!!!', np.isnan(std_errors).sum())
                
            # this makes no difference...
            err_inds = np.where(np.isnan(std_errors))
            std_errors[err_inds] = err
                        
            std_errors = np.array(std_errors,dtype=np.float64)
            
#            std_errors = np.random.normal(0,1,size=np.shape(std_errors))
                                
#            std_errors = np.ones(np.shape(std_errors)) * 100
                                
            # bounds force the ampltitude to be positive
            popt, pcov = curve_fit(superdarn_dt_sdt_tdt_mean,
                                   lons, 
                                   winds,
                                   sigma=(std_errors**1),
                                   absolute_sigma=True,
                                   maxfev=20000)
            
            # bounds=([0,0,0,-np.inf,-np.inf,-np.inf,-np.inf],
            #                               [30,30,30,np.inf,np.inf,np.inf,np.inf])
                       
            # store the resolved parameters of the fit (amplitude, phase)
            fit_vars[t - it_d,:] = popt[:6]
            errors[t - it_d,:] = np.sqrt(np.diag(pcov))[:6]
                        
            wave_0[t - it_d,:] = dt(lon_fit,popt[0],popt[3])
            wave_1[t - it_d,:] = sdt(lon_fit,popt[1],popt[4])
            wave_2[t - it_d,:] = tdt(lon_fit,popt[2],popt[5])
            
    return fit_vars, wave_0, wave_1, wave_2, errors



def superdarn_dt_sdt_tdt_mean_sw1_sw3_fitter_std(measurements,
                                                 stds,
                                                 longitudes,
                                                 x_dim,
                                                 window_w,
                                                 tstep,
                                                 data_threshold,
                                                 err):
    
    """ 
    longitudes in radians
    
    """
    
    # number of measurements (iterations)
    iters = int(np.shape(measurements)[0] - window_w)
    it_d = int(window_w / 2)
    
    # array to store the fit variables in 
    fit_vars = np.zeros((iters,10))
    errors = np.zeros((iters,10))
    fit_data = np.empty(iters)
    
    # array to store wave fit in
    wave_0 = np.zeros((iters,x_dim))
    wave_1 = np.zeros((iters,x_dim)) 
    wave_2 = np.zeros((iters,x_dim)) 
    
    # x-axis on which fit is projected
    lon_fit = np.linspace(0,2 * np.pi,x_dim,endpoint=False)
    
    t_steps = np.arange(it_d, iters + it_d)[::tstep]
    
    for t in t_steps:
        # select data and filter np.nans from longitude and data array
        v = measurements[(t - it_d):(t + it_d),:].flatten()
        
        stand_err = stds[(t - it_d):(t + it_d),:].flatten()
        
        # number of data points in fit
        fit_data[t - it_d] = int(window_w * 10 - np.isnan(v).sum())
        
        # fill longitude array with (arbitrary) 4pi spacing between t-steps
        lons = np.empty((window_w,np.shape(longitudes)[0]))
        for k in range(window_w):
            lons[k,:] = longitudes + k * 4 * np.pi

        lons.flatten()
        
        # if more nan values than threshold, reject fit
        threshold = int(window_w * 10 - data_threshold)
        
        if int(np.isnan(v).sum()) > threshold:
            # store the resolved parameters of the fit (amplitude, phase)
            fit_vars[t - it_d,:] = np.nan
            
            wave_0[t - it_d,:] = np.nan
            wave_1[t - it_d,:] = np.nan
            wave_2[t - it_d,:] = np.nan
            
            # do not take non-fits into account by setting to nan
            fit_data[t - it_d] = np.nan 
        
        else:
            inds = np.where(np.isnan(v))
            winds = np.delete(v,inds)
            lons = np.delete(lons,inds)
            
            std_errors = np.delete(stand_err,inds)
            
            # fix missing STDs by filling with above-average 'error'; this does nothing at the moment
            
#            if np.isnan(std_errors).sum() > 0:
#                print('you fool!!!', np.isnan(std_errors).sum())
                
            # this makes no difference...
            err_inds = np.where(np.isnan(std_errors))
            std_errors[err_inds] = err
                        
            std_errors = np.array(std_errors,dtype=np.float64)
            
#            std_errors = np.random.normal(0,1,size=np.shape(std_errors))
                                
#            std_errors = np.ones(np.shape(std_errors)) * 100
                                
            # bounds force the ampltitude to be positive
            popt, pcov = curve_fit(superdarn_dt_sdt_tdt_mean_sw1_sw3,
                                   lons, 
                                   winds,
                                   sigma=(std_errors**1),
                                   absolute_sigma=True,
                                   maxfev=20000)
            
            # bounds=([0,0,0,-np.inf,-np.inf,-np.inf,-np.inf],
            #                               [30,30,30,np.inf,np.inf,np.inf,np.inf])
                       
            # store the resolved parameters of the fit (amplitude, phase)
            fit_vars[t - it_d,:] = popt[:10]
            errors[t - it_d,:] = np.sqrt(np.diag(pcov))[:10]
                        
            wave_0[t - it_d,:] = dt(lon_fit,popt[0],popt[3])
            wave_1[t - it_d,:] = sdt(lon_fit,popt[1],popt[4])
            wave_2[t - it_d,:] = tdt(lon_fit,popt[2],popt[5])
            
    return fit_vars, wave_0, wave_1, wave_2, errors

#%% LEGACY FUNCTION

def legacy_updated_superdarn_mc_skim(start_date,
                                     end_date,
                                     var='v',
                                     optimized_spread=False,
                                     skim_n='25'):    
    
    """
    -- start and end date: datetime objects
    -- var:  "v" or "u"
    -- mc_skim --> skimmed based on meteor echo count rate,
       set count skim rate in input file (25,50,75,100)
    
    Sign of U has been corrected. 
    
    The (updated) pre-processing of the datafiles used in this function can be
    traced back to the Python script in the datafile's folder.
    
    Configured 11-12-2019.
    
    """

    path = '/home/wim/Desktop/NTNU/SuperDARN data/'
    stations = ['Salmon','Kodiak','PrinceGeorge','Saskatoon',
                'RankinInlet','Kapuskasing','GooseBay','Stokkseyri',
                'Pykkvibaer','Hankasalmi']
    
    # number of amplitude measurements (t_steps)
    t_steps = int((end_date - start_date) / timedelta(hours=1))

    # array to hold hourly interval measurements
    measurements = np.zeros((t_steps,np.shape(stations)[0]))
    
    # longitudes (rad) scaled to NAVGEM; corresponding to data measurements
    longitudes = np.array([-158.5,-150.1,-123.2,-105.2,-91.9,-83.3,
                           -60.3,-26.9,-18.0,25.2 - 360])
    longitudes = (longitudes + 360) / 180 * np.pi
    
    if var == 'v':
        var = 7
    elif var == 'u':
        var = 8
        
    # fill measurement field using station data from start to end date
    for index,item in enumerate(stations):
        file_mc = 'SD_wind_correct_signs_mc_' + skim_n
        print(file_mc)
        data = np.load(path + file_mc + item + '.npy')
        
        not_found = True
        c = 0
                
        while not_found:
            if datetime(int(data[c,0]),int(data[c,1]),
                        int(data[c,2]),int(data[c,3])) == start_date:
                measurements[:,index] = data[c:(c + t_steps),var]
                not_found = False
            else:
                c += 1
                
    # sample stations such as to improve equidistant longitudinal spread
    if optimized_spread:
            for t in range(t_steps):
                # choose KSR over Kod
                if np.isfinite(measurements[t,0]) and np.isfinite(measurements[t,1]):
                    measurements[t,1] = np.nan
                # choose Kap over Rkn
                if np.isfinite(measurements[t,4]) and np.isfinite(measurements[t,5]):
                    measurements[t,4] = np.nan
                # choose Sto over Pyk
                if np.isfinite(measurements[t,7]) and np.isfinite(measurements[t,8]):
                    measurements[t,8] = np.nan
                
    return measurements, longitudes
