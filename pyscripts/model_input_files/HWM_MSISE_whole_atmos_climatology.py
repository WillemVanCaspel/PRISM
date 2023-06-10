#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 19:43:09 2020

@author: wim
"""

from datetime import datetime, timedelta
import numpy as np
import model_funcs as mf
from matplotlib import pyplot as plt
from wrf import interplevel

#%%

lats = np.linspace(-89,89,90,endpoint=True) # index 0 = SP, as in DOPE
lons = np.arange(0,360)[::30] 

folder = '/home/wim/Desktop/Projects/Model/HWM14_MSISE00/indices/'

altitudes = [75,77.5,80,82.5,85,87.5,90,92.5,95,97.5,100,102.5,105,107.5,110,115,120,125,130,135,140,145,150,
             160,170,180,190,200,
             225,250,275,300,325,350,375,400][::-1]

# top to bottom thermosphere to troposphere
altitudes = np.concatenate((altitudes,np.arange(72.5,-0.1,-1.25)))

# 365 days of zonal mean fields
output = np.empty((4,365 * 8,np.shape(altitudes)[0],np.shape(lats)[0]))

# temporarily store output arrays in
temp_output = np.empty((4,8,np.shape(altitudes)[0],np.shape(lats)[0]))

# start day to base climatology on
date_temp = datetime(2010,1,1,1)

day_s = 0
n_days = 365

for day_step in range(day_s,n_days):
    for t_step in range(8):
        for index,item in enumerate(altitudes):

            u,v,t,p = mf.hwm_msise_uvtp_climatology(date_temp,
                                                    item,
                                                    lats,
                                                    lons,
                                                    folder)
            
            # zonal means
            temp_output[0,t_step,index,:] = np.mean(u,axis=1)
            temp_output[1,t_step,index,:] = np.mean(v,axis=1)
            temp_output[2,t_step,index,:] = np.mean(t,axis=1)
            temp_output[3,t_step,index,:] = np.mean(p,axis=1)
            
        date_temp = date_temp + timedelta(hours=3)
                
    # daily means
    output[0,day_step,:,:] = np.mean(temp_output[0,:,:,:],axis=0)
    output[1,day_step,:,:] = np.mean(temp_output[1,:,:,:],axis=0)
    output[2,day_step,:,:] = np.mean(temp_output[2,:,:,:],axis=0)
    output[3,day_step,:,:] = np.mean(temp_output[3,:,:,:],axis=0)
    
    print('Finished DOY: ', day_step)
    
# np.save('/media/wim/datadisk/HWM_MSISE/climatology.npy',output)
    
#%%

u = np.load('/media/wim/datadisk/HWM_MSISE/climatology.npy')[0,:,:,:]
v = np.load('/media/wim/datadisk/HWM_MSISE/climatology.npy')[1,:,:,:]
t = np.load('/media/wim/datadisk/HWM_MSISE/climatology.npy')[2,:,:,:]
p = np.load('/media/wim/datadisk/HWM_MSISE/climatology.npy')[3,:,:,:]

# set fixed surface pressure so that interpolation is always within range
p[:,-1,:] = 100005

#%%

plt.pcolormesh(t[0,::-1,:])
    
#%%

# interpolate to sigprim model levels 
dope_lvls = mf.sigprim_plvls() # Pa

# 31 day padding
full_u = np.empty((365 + 31,np.shape(dope_lvls)[0],90,90))
full_v = np.empty((365 + 31,np.shape(dope_lvls)[0],90,90))
full_t = np.empty((365 + 31,np.shape(dope_lvls)[0],90,90))

for d in range(31):
    u_temp = np.empty((95,90,90))
    v_temp = np.empty((95,90,90))
    t_temp = np.empty((95,90,90))
    p_temp = np.empty((95,90,90))
    
    # populate longitudes
    for l in range(90):
        u_temp[:,:,l] = u[0,:,:]
        v_temp[:,:,l] = v[0,:,:]
        t_temp[:,:,l] = t[0,:,:]
        p_temp[:,:,l] = p[0,:,:]
    
    interp_u = interplevel(u_temp,p_temp,dope_lvls)
    interp_v = interplevel(v_temp,p_temp,dope_lvls)
    interp_t = interplevel(t_temp,p_temp,dope_lvls)
    
    full_u[d,:,:,:] = interp_u
    full_v[d,:,:,:] = interp_v
    full_t[d,:,:,:] = interp_t
    
    print('Constructing padded day: ', d)

for d in range(31, 365 + 31):
    u_temp = np.empty((95,90,90))
    v_temp = np.empty((95,90,90))
    t_temp = np.empty((95,90,90))
    p_temp = np.empty((95,90,90))
    
    # populate longitudes
    for l in range(90):
        u_temp[:,:,l] = u[d - 31,:,:]
        v_temp[:,:,l] = v[d - 31,:,:]
        t_temp[:,:,l] = t[d - 31,:,:]
        p_temp[:,:,l] = p[d - 31,:,:]
    
    interp_u = interplevel(u_temp,p_temp,dope_lvls)
    interp_v = interplevel(v_temp,p_temp,dope_lvls)
    interp_t = interplevel(t_temp,p_temp,dope_lvls)
    
    full_u[d,:,:,:] = interp_u
    full_v[d,:,:,:] = interp_v
    full_t[d,:,:,:] = interp_t
    
    print('Constructing DOY: ', d - 31)
    
#%%

temp_v = np.copy(full_v[:,15:,:,:])
plt.contourf(temp_v[0 + 31,::-1,:,0])
plt.colorbar()

#%%

# write to file
z_in = dope_lvls / 1000 # dummy array for backwards compatibility
time_grid = np.arange(-31,365,1,dtype=np.float32)

#%%

dphi = (2) / 180 * np.pi
dlam = (4) / 180 * np.pi

# save to file
folder = '/home/wim/Fortran/Sigpy/input/'

mf.uvt_to_sigprim_zm_div(full_u,
                         full_v,
                         full_t,
                         dphi,
                         dlam,
                         z_in,
                         time_grid,
                         folder,
                         'zm_climatology_div2.dat')

#%%

