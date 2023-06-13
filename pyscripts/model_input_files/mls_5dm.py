#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 19:03:01 2020

@author: wim

"""

import h5py
import numpy as np
from matplotlib import pyplot as plt
from netCDF4 import Dataset
from wrf import interplevel

folder = '/media/wim/datadisk/ERA5 Netcdf/UVZT_2014/'
data = Dataset(folder + '201401.nc',mode='r')

# flip latitude to match NAVGEM-HA and HWM+MSISE
era_lats = np.array(data['latitude'][::-1])

#%%

f = h5py.File('/home/wim/Desktop/MLS-Aura_L2GP-Temperature_v03-33-c01_2013d001.he5', 'r')

swaths = f['HDFEOS']['SWATHS']
temps = swaths['Temperature']
latitudes = swaths['Temperature']['Geolocation Fields']['Latitude']

# pressure levels of temperature measurements
pressure = swaths['Temperature']['Geolocation Fields']['Pressure'][:]

# MLS documentation: impact of clouds neglibile above 100 hPa
range_bot = 150   # hPa, ~ 13 km altitude; stratosphere
range_top = 0.001 # hPa top to where thermal wind balance is a decent approx

ind_bot = np.argmin(np.abs(pressure - range_bot))
ind_top = np.argmin(np.abs(pressure - range_top)) + 1

# pressure range of meaningful data and number of pressure levels
pressure = pressure[ind_bot:ind_top]
n_levs = np.shape(pressure)[0]

# number of measurements
n_steps = np.shape(temps['nTimes'])[0]

status = swaths['Temperature']['Data Fields']['Status']
temperature = swaths['Temperature']['Data Fields']['Temperature']

precision = swaths['Temperature']['Data Fields']['TemperaturePrecision']
quality = swaths['Temperature']['Data Fields']['Quality']
convergence = swaths['Temperature']['Data Fields']['Convergence']

# populate measurement field
dat_arr = np.empty((n_steps,n_levs))
lat_arr = np.empty(n_steps)

# do not include bad retrievals by documentation selection criterea
for m in range(n_steps):
    if status[m] == 0 and quality[m] > 0.65 and convergence[m] < 1.2:
        
        temp_swath = temperature[m,ind_bot:ind_top]
        prec_swath = precision[m,ind_bot:ind_top]
        
        # only data where precision is a positive number
        inds = np.where(prec_swath < 0)
        temp_swath[inds] = np.nan
        
        dat_arr[m,:] = temp_swath
        lat_arr[m] = latitudes[m]
        
    else:
        dat_arr[m,:] = np.nan
        lat_arr[m] = np.nan

nan_inds = np.where(np.isnan(lat_arr))
dat_arr = np.delete(dat_arr,nan_inds,axis=0)
lat_arr = np.delete(lat_arr,nan_inds)
        
inds = np.argsort(lat_arr)
sorted_tems = dat_arr[inds,:]
sorted_lats = lat_arr[inds]

# daily mean measurements over latitude bins
lat_bins = np.arange(-81,82,2.0) # means centered on lat (degrees)
mean_t = np.empty((np.shape(sorted_tems)[1],np.shape(lat_bins)[0]))

k = 0
m = 0

for index,item in enumerate(lat_bins):
    while sorted_lats[k] < lat_bins[index]:
        k = k + 1     
    lat_bin = sorted_tems[m:k,:]
    mean_t[:,index] = np.nanmean(lat_bin,axis=0)
    m = np.copy(k)

# glue on the data past 81 degrees latitude as a final lat / bin step
mean_t = np.column_stack((mean_t,np.nanmean(sorted_tems[m:,:],axis=0)))
lat_bins = np.concatenate((lat_bins,[83]))

# extend data to the poles and center bin coordinates
mls_lat_coords = lat_bins - 1.0
mls_lat_coords = np.concatenate(([-91],mls_lat_coords,[91]))
mean_t = np.column_stack((mean_t[:,0],mean_t,mean_t[:,-1]))

#interpolate to era5 latitude grid
interp_field_1 = np.zeros((np.shape(mean_t)[0],np.shape(era_lats)[0]))

for i in range(np.shape(interp_field_1)[0]):
    interp_field_1[i,:] = np.interp(era_lats,mls_lat_coords,mean_t[i,:])
    

f = h5py.File('/home/wim/Desktop/MLS-Aura_L2GP-Temperature_v03-33-c01_2013d004.he5', 'r')

swaths = f['HDFEOS']['SWATHS']
temps = swaths['Temperature']
latitudes = swaths['Temperature']['Geolocation Fields']['Latitude']

# pressure levels of temperature measurements
pressure = swaths['Temperature']['Geolocation Fields']['Pressure'][:]

# MLS documentation: impact of clouds neglibile above 100 hPa
range_bot = 150   # hPa, ~ 13 km altitude; stratosphere
range_top = 0.001 # hPa top to where thermal wind balance is a decent approx

ind_bot = np.argmin(np.abs(pressure - range_bot))
ind_top = np.argmin(np.abs(pressure - range_top)) + 1

# pressure range of meaningful data and number of pressure levels
pressure = pressure[ind_bot:ind_top]
n_levs = np.shape(pressure)[0]

# number of measurements
n_steps = np.shape(temps['nTimes'])[0]

status = swaths['Temperature']['Data Fields']['Status']
temperature = swaths['Temperature']['Data Fields']['Temperature']

precision = swaths['Temperature']['Data Fields']['TemperaturePrecision']
quality = swaths['Temperature']['Data Fields']['Quality']
convergence = swaths['Temperature']['Data Fields']['Convergence']

# populate measurement field
dat_arr = np.empty((n_steps,n_levs))
lat_arr = np.empty(n_steps)

# do not include bad retrievals by documentation selection criterea
for m in range(n_steps):
    if status[m] == 0 and quality[m] > 0.65 and convergence[m] < 1.2:
        
        temp_swath = temperature[m,ind_bot:ind_top]
        prec_swath = precision[m,ind_bot:ind_top]
        
        # only data where precision is a positive number
        inds = np.where(prec_swath < 0)
        temp_swath[inds] = np.nan
        
        dat_arr[m,:] = temp_swath
        lat_arr[m] = latitudes[m]
        
    else:
        dat_arr[m,:] = np.nan
        lat_arr[m] = np.nan

nan_inds = np.where(np.isnan(lat_arr))
dat_arr = np.delete(dat_arr,nan_inds,axis=0)
lat_arr = np.delete(lat_arr,nan_inds)
        
inds = np.argsort(lat_arr)
sorted_tems = dat_arr[inds,:]
sorted_lats = lat_arr[inds]

# daily mean measurements over latitude bins
lat_bins = np.arange(-81,82,2.0) # means centered on lat (degrees)
mean_t = np.empty((np.shape(sorted_tems)[1],np.shape(lat_bins)[0]))

k = 0
m = 0

for index,item in enumerate(lat_bins):
    while sorted_lats[k] < lat_bins[index]:
        k = k + 1     
    lat_bin = sorted_tems[m:k,:]
    mean_t[:,index] = np.nanmean(lat_bin,axis=0)
    m = np.copy(k)

# glue on the data past 81 degrees latitude as a final lat / bin step
mean_t = np.column_stack((mean_t,np.nanmean(sorted_tems[m:,:],axis=0)))
lat_bins = np.concatenate((lat_bins,[83]))

# extend data to the poles and center bin coordinates
mls_lat_coords = lat_bins - 1.0
mls_lat_coords = np.concatenate(([-91],mls_lat_coords,[91]))
mean_t = np.column_stack((mean_t[:,0],mean_t,mean_t[:,-1]))

#interpolate to era5 latitude grid
interp_field_2 = np.zeros((np.shape(mean_t)[0],np.shape(era_lats)[0]))

for i in range(np.shape(interp_field_1)[0]):
    interp_field_2[i,:] = np.interp(era_lats,mls_lat_coords,mean_t[i,:])
    
f = h5py.File('/home/wim/Desktop/MLS-Aura_L2GP-Temperature_v03-33-c01_2013d008.he5', 'r')

swaths = f['HDFEOS']['SWATHS']
temps = swaths['Temperature']
latitudes = swaths['Temperature']['Geolocation Fields']['Latitude']

# pressure levels of temperature measurements
pressure = swaths['Temperature']['Geolocation Fields']['Pressure'][:]

# MLS documentation: impact of clouds neglibile above 100 hPa
range_bot = 150   # hPa, ~ 13 km altitude; stratosphere
range_top = 0.001 # hPa top to where thermal wind balance is a decent approx

ind_bot = np.argmin(np.abs(pressure - range_bot))
ind_top = np.argmin(np.abs(pressure - range_top)) + 1

# pressure range of meaningful data and number of pressure levels
pressure = pressure[ind_bot:ind_top]
n_levs = np.shape(pressure)[0]

# number of measurements
n_steps = np.shape(temps['nTimes'])[0]

status = swaths['Temperature']['Data Fields']['Status']
temperature = swaths['Temperature']['Data Fields']['Temperature']

precision = swaths['Temperature']['Data Fields']['TemperaturePrecision']
quality = swaths['Temperature']['Data Fields']['Quality']
convergence = swaths['Temperature']['Data Fields']['Convergence']

# populate measurement field
dat_arr = np.empty((n_steps,n_levs))
lat_arr = np.empty(n_steps)

# do not include bad retrievals by documentation selection criterea
for m in range(n_steps):
    if status[m] == 0 and quality[m] > 0.65 and convergence[m] < 1.2:
        
        temp_swath = temperature[m,ind_bot:ind_top]
        prec_swath = precision[m,ind_bot:ind_top]
        
        # only data where precision is a positive number
        inds = np.where(prec_swath < 0)
        temp_swath[inds] = np.nan
        
        dat_arr[m,:] = temp_swath
        lat_arr[m] = latitudes[m]
        
    else:
        dat_arr[m,:] = np.nan
        lat_arr[m] = np.nan

nan_inds = np.where(np.isnan(lat_arr))
dat_arr = np.delete(dat_arr,nan_inds,axis=0)
lat_arr = np.delete(lat_arr,nan_inds)
        
inds = np.argsort(lat_arr)
sorted_tems = dat_arr[inds,:]
sorted_lats = lat_arr[inds]

# daily mean measurements over latitude bins
lat_bins = np.arange(-81,82,2.0) # means centered on lat (degrees)
mean_t = np.empty((np.shape(sorted_tems)[1],np.shape(lat_bins)[0]))

k = 0
m = 0

for index,item in enumerate(lat_bins):
    while sorted_lats[k] < lat_bins[index]:
        k = k + 1     
    lat_bin = sorted_tems[m:k,:]
    mean_t[:,index] = np.nanmean(lat_bin,axis=0)
    m = np.copy(k)

# glue on the data past 81 degrees latitude as a final lat / bin step
mean_t = np.column_stack((mean_t,np.nanmean(sorted_tems[m:,:],axis=0)))
lat_bins = np.concatenate((lat_bins,[83]))

# extend data to the poles and center bin coordinates
mls_lat_coords = lat_bins - 1.0
mls_lat_coords = np.concatenate(([-91],mls_lat_coords,[91]))
mean_t = np.column_stack((mean_t[:,0],mean_t,mean_t[:,-1]))

#interpolate to era5 latitude grid
interp_field_3 = np.zeros((np.shape(mean_t)[0],np.shape(era_lats)[0]))

for i in range(np.shape(interp_field_1)[0]):
    interp_field_3[i,:] = np.interp(era_lats,mls_lat_coords,mean_t[i,:])
    
f = h5py.File('/home/wim/Desktop/MLS-Aura_L2GP-Temperature_v03-33-c01_2013d012.he5', 'r')

swaths = f['HDFEOS']['SWATHS']
temps = swaths['Temperature']
latitudes = swaths['Temperature']['Geolocation Fields']['Latitude']

# pressure levels of temperature measurements
pressure = swaths['Temperature']['Geolocation Fields']['Pressure'][:]

# MLS documentation: impact of clouds neglibile above 100 hPa
range_bot = 150   # hPa, ~ 13 km altitude; stratosphere
range_top = 0.001 # hPa top to where thermal wind balance is a decent approx

ind_bot = np.argmin(np.abs(pressure - range_bot))
ind_top = np.argmin(np.abs(pressure - range_top)) + 1

# pressure range of meaningful data and number of pressure levels
pressure = pressure[ind_bot:ind_top]
n_levs = np.shape(pressure)[0]

# number of measurements
n_steps = np.shape(temps['nTimes'])[0]

status = swaths['Temperature']['Data Fields']['Status']
temperature = swaths['Temperature']['Data Fields']['Temperature']

precision = swaths['Temperature']['Data Fields']['TemperaturePrecision']
quality = swaths['Temperature']['Data Fields']['Quality']
convergence = swaths['Temperature']['Data Fields']['Convergence']

# populate measurement field
dat_arr = np.empty((n_steps,n_levs))
lat_arr = np.empty(n_steps)

# do not include bad retrievals by documentation selection criterea
for m in range(n_steps):
    if status[m] == 0 and quality[m] > 0.65 and convergence[m] < 1.2:
        
        temp_swath = temperature[m,ind_bot:ind_top]
        prec_swath = precision[m,ind_bot:ind_top]
        
        # only data where precision is a positive number
        inds = np.where(prec_swath < 0)
        temp_swath[inds] = np.nan
        
        dat_arr[m,:] = temp_swath
        lat_arr[m] = latitudes[m]
        
    else:
        dat_arr[m,:] = np.nan
        lat_arr[m] = np.nan

nan_inds = np.where(np.isnan(lat_arr))
dat_arr = np.delete(dat_arr,nan_inds,axis=0)
lat_arr = np.delete(lat_arr,nan_inds)
        
inds = np.argsort(lat_arr)
sorted_tems = dat_arr[inds,:]
sorted_lats = lat_arr[inds]

# daily mean measurements over latitude bins
lat_bins = np.arange(-81,82,2.0) # means centered on lat (degrees)
mean_t = np.empty((np.shape(sorted_tems)[1],np.shape(lat_bins)[0]))

k = 0
m = 0

for index,item in enumerate(lat_bins):
    while sorted_lats[k] < lat_bins[index]:
        k = k + 1     
    lat_bin = sorted_tems[m:k,:]
    mean_t[:,index] = np.nanmean(lat_bin,axis=0)
    m = np.copy(k)

# glue on the data past 81 degrees latitude as a final lat / bin step
mean_t = np.column_stack((mean_t,np.nanmean(sorted_tems[m:,:],axis=0)))
lat_bins = np.concatenate((lat_bins,[83]))

# extend data to the poles and center bin coordinates
mls_lat_coords = lat_bins - 1.0
mls_lat_coords = np.concatenate(([-91],mls_lat_coords,[91]))
mean_t = np.column_stack((mean_t[:,0],mean_t,mean_t[:,-1]))

#interpolate to era5 latitude grid
interp_field_4 = np.zeros((np.shape(mean_t)[0],np.shape(era_lats)[0]))

for i in range(np.shape(interp_field_1)[0]):
    interp_field_4[i,:] = np.interp(era_lats,mls_lat_coords,mean_t[i,:])

f = h5py.File('/home/wim/Desktop/MLS-Aura_L2GP-Temperature_v03-33-c01_2013d016.he5', 'r')

swaths = f['HDFEOS']['SWATHS']
temps = swaths['Temperature']
latitudes = swaths['Temperature']['Geolocation Fields']['Latitude']

# pressure levels of temperature measurements
pressure = swaths['Temperature']['Geolocation Fields']['Pressure'][:]

# MLS documentation: impact of clouds neglibile above 100 hPa
range_bot = 150   # hPa, ~ 13 km altitude; stratosphere
range_top = 0.001 # hPa top to where thermal wind balance is a decent approx

ind_bot = np.argmin(np.abs(pressure - range_bot))
ind_top = np.argmin(np.abs(pressure - range_top)) + 1

# pressure range of meaningful data and number of pressure levels
pressure = pressure[ind_bot:ind_top]
n_levs = np.shape(pressure)[0]

# number of measurements
n_steps = np.shape(temps['nTimes'])[0]

status = swaths['Temperature']['Data Fields']['Status']
temperature = swaths['Temperature']['Data Fields']['Temperature']

precision = swaths['Temperature']['Data Fields']['TemperaturePrecision']
quality = swaths['Temperature']['Data Fields']['Quality']
convergence = swaths['Temperature']['Data Fields']['Convergence']

# populate measurement field
dat_arr = np.empty((n_steps,n_levs))
lat_arr = np.empty(n_steps)

# do not include bad retrievals by documentation selection criterea
for m in range(n_steps):
    if status[m] == 0 and quality[m] > 0.65 and convergence[m] < 1.2:
        
        temp_swath = temperature[m,ind_bot:ind_top]
        prec_swath = precision[m,ind_bot:ind_top]
        
        # only data where precision is a positive number
        inds = np.where(prec_swath < 0)
        temp_swath[inds] = np.nan
        
        dat_arr[m,:] = temp_swath
        lat_arr[m] = latitudes[m]
        
    else:
        dat_arr[m,:] = np.nan
        lat_arr[m] = np.nan

nan_inds = np.where(np.isnan(lat_arr))
dat_arr = np.delete(dat_arr,nan_inds,axis=0)
lat_arr = np.delete(lat_arr,nan_inds)
        
inds = np.argsort(lat_arr)
sorted_tems = dat_arr[inds,:]
sorted_lats = lat_arr[inds]

# daily mean measurements over latitude bins
lat_bins = np.arange(-81,82,2.0) # means centered on lat (degrees)
mean_t = np.empty((np.shape(sorted_tems)[1],np.shape(lat_bins)[0]))

k = 0
m = 0

for index,item in enumerate(lat_bins):
    while sorted_lats[k] < lat_bins[index]:
        k = k + 1     
    lat_bin = sorted_tems[m:k,:]
    mean_t[:,index] = np.nanmean(lat_bin,axis=0)
    m = np.copy(k)

# glue on the data past 81 degrees latitude as a final lat / bin step
mean_t = np.column_stack((mean_t,np.nanmean(sorted_tems[m:,:],axis=0)))
lat_bins = np.concatenate((lat_bins,[83]))

# extend data to the poles and center bin coordinates
mls_lat_coords = lat_bins - 1.0
mls_lat_coords = np.concatenate(([-91],mls_lat_coords,[91]))
mean_t = np.column_stack((mean_t[:,0],mean_t,mean_t[:,-1]))

#interpolate to era5 latitude grid
interp_field_5 = np.zeros((np.shape(mean_t)[0],np.shape(era_lats)[0]))

for i in range(np.shape(interp_field_1)[0]):
    interp_field_5[i,:] = np.interp(era_lats,mls_lat_coords,mean_t[i,:])

#%%

interp_field = (interp_field_1 + interp_field_2 + interp_field_3 +
                interp_field_4 + interp_field_5) / 5
    
data = Dataset(folder + '2014' + '01' + '.nc',mode='r')

era_pres = data['level'][::-1] # hPa, like MLS
era_mls_ind = np.argmin(np.abs(era_pres - range_bot)) # 150 hPa index

t = 0

# daily mean zonal mean era5 fields
u = np.mean(np.mean(data['u'][int(t * 8):int(t * 8 + 8),::-1,:,:],axis=0),axis=2)
t = np.mean(np.mean(data['t'][int(t * 8):int(t * 8 + 8),::-1,:,:],axis=0),axis=2)

# these era5 fields are glued on beneath the MLS data
u_stack_bot = u[:era_mls_ind,:]
t_stack_bot = t[:era_mls_ind,:]

#%%

u_1 = u[era_mls_ind,:] # at 150 hpa
u_0 = u[era_mls_ind - 1,:] # at 175 hpa

#%%

def coriolis(lat):
    lat = lat / 180 * np.pi
    omega = 7.2921159e-5 # angular velocity earth
    return 2 * omega * np.sin(lat)

def beta(lat):
    lat = lat / 180 * np.pi
    omega = 7.2921159e-5 # angular velocity earth
    a = 6.3781e6 # radius earth [m]
    return 2 * omega * np.cos(lat) / a

def cart_y(lat):
    lat = lat / 180 * np.pi
    a = 6.3781e6 # radius earth [m]
    return a * np.sin(lat)
    
cor = coriolis(era_lats)  # coriolis parameter
bet = beta(era_lats)      # beta-plane parameter
y_cart = cart_y(era_lats) # cartesian y-coordinate
gas_c = 287.058 # dry gas constant

#%%

beta_plane_lat = 2

pres = np.array(np.concatenate((era_pres[:era_mls_ind],pressure)))*1e2 # Pa

t_field = np.concatenate((t_stack_bot,interp_field))
u_field = np.zeros(np.shape(t_field))
u_field[:era_mls_ind,:] = u_stack_bot

sp_lim = 5
np_lim = 5

# central du / dp 
for l in range(era_mls_ind,np.shape(t_field)[0]):
    for n in range(sp_lim,np.shape(t_field)[1] - np_lim):
        # forward temperature at SP
        if n == sp_lim:
            dt = t_field[l - 1,n + 1] - t_field[l - 1,n]        
            dy = y_cart[n + 1] - y_cart[n]
            dp = pres[l] - pres[l - 2]
            u_field[l,n] = gas_c * (cor[n] * pres[l - 1])**-1 * dt / dy * dp + u_field[l - 2,n]
        # backward temperature at NP
        elif n == (np.shape(t_field)[1] - np_lim):
            dt = t_field[l - 1,n] - t_field[l - 1,n - 1]        
            dy = y_cart[n] - y_cart[n - 1]
            dp = pres[l] - pres[l - 2]
            u_field[l,n] = gas_c * (cor[n] * pres[l - 1])**-1 * dt / dy * dp + u_field[l - 2,n]
        # l'hopital at equator
        elif era_lats[n] == 0:
            dt = t_field[l - 1,n + 1] - 2 * t_field[l - 1,n] + t_field[l - 1,n - 1]        
            dy = (y_cart[n + 1] - y_cart[n]) ** 4
            dp = pres[l] - pres[l - 2]
            u_field[l,n] = gas_c * (bet[n] * pres[l - 1])**-1 * dt / dy * dp + u_field[l - 2,n]          
        # equatorial beta plane within equatorial bounds
        elif np.abs(era_lats[n]) < beta_plane_lat:
            None
        # central difference 
        else:
            dt = t_field[l - 1,n + 1] - t_field[l - 1,n - 1]        
            dy = y_cart[n + 1] - y_cart[n - 1]
            dp = pres[l] - pres[l - 2]
            u_field[l,n] = gas_c * (cor[n] * pres[l - 1])**-1 * dt / dy * dp + u_field[l - 2,n]

plt.pcolormesh(u_field[:,:])
plt.colorbar()

#%%

beta_plane_lat = 5

pres = np.array(np.concatenate((era_pres[:era_mls_ind],pressure))) * 1e2 # Pa

t_field = np.concatenate((t_stack_bot,interp_field))
u_field = np.zeros(np.shape(t_field))
u_field[:era_mls_ind,:] = u_stack_bot

sp_lim = 10
np_lim = 10


# central du / dp 
for l in range(era_mls_ind,np.shape(t_field)[0]):
    for n in range(sp_lim,np.shape(t_field)[1] - np_lim):
        # forward temperature at SP
        if n == sp_lim:
            dt = t_field[l - 1,n + 2] - t_field[l - 1,n]        
            dy = y_cart[n + 2] - y_cart[n]
            dp = pres[l] - pres[l - 1]
            u_field[l,n] = gas_c * (cor[n] * pres[l - 1])**-1 * dt / dy * dp + u_field[l - 1,n]
        # backward temperature at NP
        elif n == (np.shape(t_field)[1] - np_lim):
            dt = t_field[l - 1,n] - t_field[l - 1,n - 2]        
            dy = y_cart[n] - y_cart[n - 2]
            dp = pres[l] - pres[l - 1]
            u_field[l,n] = gas_c * (cor[n] * pres[l - 1])**-1 * dt / dy * dp + u_field[l - 1,n]
        # l'hopital at equator
        elif era_lats[n] == 0:
            dt = t_field[l - 1,n + 5] - 2 * t_field[l - 1,n] + t_field[l - 1,n - 5]
            dy = (y_cart[n + 5] - y_cart[n]) ** 2
            dp = pres[l] - pres[l - 1]
            u_field[l,n] = gas_c * (bet[n] * pres[l - 1])**-1 * dt / dy * dp + u_field[l - 1,n]          
        # equatorial beta plane within equatorial bounds
        elif np.abs(era_lats[n]) < beta_plane_lat:
            dt = t_field[l - 1,n + 2] - t_field[l - 1,n - 2]        
            dy = y_cart[n + 10] - y_cart[n - 10]
            dp = pres[l] - pres[l - 1]
            u_field[l,n] = gas_c * (bet[91] * y_cart[n] * pres[l - 1] * 2)**-1 * dt / dy * dp + u_field[l - 1,n]
        # central difference 
        else:
            dt = t_field[l - 1,n + 5] - t_field[l - 1,n - 5]        
            dy = y_cart[n + 5] - y_cart[n - 5]
            dp = pres[l] - pres[l - 1]
            u_field[l,n] = gas_c * (cor[n] * pres[l - 1])**-1 * dt / dy * dp + u_field[l - 1,n]

plt.pcolormesh(u_field[:,:])
plt.colorbar()

#%%


beta_plane_lat = 5

pres = np.array(np.concatenate((era_pres[:era_mls_ind],pressure)))*1e2 # Pa
pres = np.log(pres)
t_field = np.concatenate((t_stack_bot,interp_field))
u_field = np.zeros(np.shape(t_field))
u_field[:era_mls_ind,:] = u_stack_bot

sp_lim = 10
np_lim = 10

# central du / dp 
for l in range(era_mls_ind,np.shape(t_field)[0]):
    for n in range(sp_lim,np.shape(t_field)[1] - np_lim):
        # forward temperature at SP
        if n == sp_lim:
            dt = t_field[l - 1,n + 2] - t_field[l - 1,n]        
            dy = y_cart[n + 2] - y_cart[n]
            dp = pres[l] - pres[l - 1]
            u_field[l,n] = gas_c * (cor[n])**-1 * dt / dy * dp + u_field[l - 1,n]
        # backward temperature at NP
        elif n == (np.shape(t_field)[1] - np_lim):
            dt = t_field[l - 1,n] - t_field[l - 1,n - 2]        
            dy = y_cart[n] - y_cart[n - 2]
            dp = pres[l] - pres[l - 1]
            u_field[l,n] = gas_c * (cor[n])**-1 * dt / dy * dp + u_field[l - 1,n]
        # l'hopital at equator
        elif era_lats[n] == 0:
            dt = t_field[l - 1,n + 5] - 2 * t_field[l - 1,n] + t_field[l - 1,n - 5]
            dy = (y_cart[n + 5] - y_cart[n]) ** 2
            dp = pres[l] - pres[l - 1]
            u_field[l,n] = gas_c * (bet[n])**-1 * dt / dy * dp + u_field[l - 1,n]          
        # equatorial beta plane within equatorial bounds
        elif np.abs(era_lats[n]) < beta_plane_lat:
            dt = t_field[l - 1,n + 5] - t_field[l - 1,n - 5]        
            dy = y_cart[n + 5] - y_cart[n - 5]
            dp = pres[l] - pres[l - 1]
            u_field[l,n] = gas_c * (bet[n] * y_cart[n])**-1 * dt / dy * dp + u_field[l - 1,n]
        # central difference 
        else:
            dt = t_field[l - 1,n + 5] - t_field[l - 1,n - 5]        
            dy = y_cart[n + 5] - y_cart[n - 5]
            dp = pres[l] - pres[l - 1]
            u_field[l,n] = gas_c * (cor[n])**-1 * dt / dy * dp + u_field[l - 1,n]

plt.pcolormesh(u_field[:,:])
plt.colorbar()

#%%

beta_plane_lat = 5

pres = np.array(np.concatenate((era_pres[:era_mls_ind],pressure)))*1e2 # Pa
pres = np.log(pres)
t_field = np.concatenate((t_stack_bot,interp_field))
u_field = np.zeros(np.shape(t_field))
u_field[:era_mls_ind,:] = u_stack_bot

sp_lim = 5
np_lim = 5

# central du / dp 
for l in range(era_mls_ind,np.shape(t_field)[0]):
    for n in range(sp_lim,np.shape(t_field)[1] - np_lim):
        # forward temperature at SP
        if n == sp_lim:
            dt = t_field[l - 1,n + 2] - t_field[l - 1,n]        
            dy = y_cart[n + 2] - y_cart[n]
            dp = pres[l] - pres[l - 1]
            u_field[l,n] = gas_c * (cor[n])**-1 * dt / dy * dp + u_field[l - 1,n]
        # backward temperature at NP
        elif n == (np.shape(t_field)[1] - np_lim):
            dt = t_field[l - 1,n] - t_field[l - 1,n - 2]        
            dy = y_cart[n] - y_cart[n - 2]
            dp = pres[l] - pres[l - 1]
            u_field[l,n] = gas_c * (cor[n])**-1 * dt / dy * dp + u_field[l - 1,n]
        # l'hopital at equator
        elif era_lats[n] == 0:
            dt = t_field[l - 1,n + 5] - 2 * t_field[l - 1,n] + t_field[l - 1,n - 5]
            dy = (y_cart[n + 5] - y_cart[n]) ** 2
            dp = pres[l] - pres[l - 1]
            u_field[l,n] = gas_c * (bet[n] * 100)**-1 * dt / dy * dp + u_field[l - 1,n]          
        # equatorial beta plane within equatorial bounds
        elif np.abs(era_lats[n]) < beta_plane_lat:
            dt = t_field[l - 1,n + 5] - t_field[l - 1,n - 5]        
            dy = y_cart[n + 5] - y_cart[n - 5]
            dp = pres[l] - pres[l - 1]
            u_field[l,n] = gas_c * (bet[n] * y_cart[n] * 50)**-1 * dt / dy * dp + u_field[l - 1,n]
        # central difference 
        else:
            dt = t_field[l - 1,n + 5] - t_field[l - 1,n - 5]        
            dy = y_cart[n + 5] - y_cart[n - 5]
            dp = pres[l] - pres[l - 1]
            u_field[l,n] = gas_c * (cor[n])**-1 * dt / dy * dp + u_field[l - 1,n]

plt.pcolormesh(u_field[:,:])
plt.colorbar()

#%%

beta_plane_lat = 0

pres = np.array(np.concatenate((era_pres[:era_mls_ind],pressure)))*1e2 # Pa
pres = np.log(pres)

t_field = np.concatenate((t_stack_bot,interp_field))
u_field = np.zeros(np.shape(t_field))
u_field[:era_mls_ind,:] = u_stack_bot

sp_lim = 10
np_lim = 10

# central du / dp 
for l in range(era_mls_ind - 1,np.shape(t_field)[0] - 1):
    for n in range(sp_lim,np.shape(t_field)[1] - np_lim):
        # forward temperature at SP
        if n == sp_lim:
            dt = t_field[l - 1,n + 2] - t_field[l - 1,n]        
            dy = y_cart[n + 2] - y_cart[n]
            dp = pres[l] - pres[l - 2]
            u_field[l,n] = gas_c * cor[n]**-1 * dt / dy * dp + u_field[l - 2,n]
        # backward temperature at NP
        elif n == (np.shape(t_field)[1] - np_lim):
            dt = t_field[l - 1,n] - t_field[l - 1,n - 2]        
            dy = y_cart[n] - y_cart[n - 2]
            dp = pres[l] - pres[l - 2]
            u_field[l,n] = gas_c * cor[n]**-1 * dt / dy * dp + u_field[l - 2,n]
        # l'hopital at equator
        elif era_lats[n] == 0:
            dt = t_field[l - 1,n + 5] - 2 * t_field[l - 1,n] + t_field[l - 1,n - 5]
            dy = (y_cart[n + 5] - y_cart[n]) ** 2
            dp = pres[l] - pres[l - 2]
            u_field[l,n] = gas_c * bet[n]**-1 * dt / dy * dp + u_field[l - 2,n]          
        # equatorial beta plane within equatorial bounds
        elif np.abs(era_lats[n]) < beta_plane_lat:
            dt = t_field[l - 1,n + 5] - t_field[l - 1,n - 5]        
            dy = y_cart[n + 5] - y_cart[n - 5]
            dp = pres[l] - pres[l - 1]
            u_field[l,n] = gas_c * (bet[n] * y_cart[n])**-1 * dt / dy * dp + u_field[l - 1,n]
        # central difference 
        else:
            dt = t_field[l - 1,n + 5] - t_field[l - 1,n - 5]        
            dy = y_cart[n + 5] - y_cart[n - 5]
            dp = pres[l + 1] - pres[l]
            u_field[l + 1,n] = gas_c * cor[n]**-1 * dt / dy * dp + u_field[l,n]

plt.pcolormesh(u_field[:,:])
plt.colorbar()

