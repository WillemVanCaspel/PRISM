#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 19:43:09 2020

@author: wim
"""

from datetime import datetime, timedelta
import numpy as np
import model_funcs as mf
import timeit

start = timeit.default_timer()

lats = np.linspace(-89,89,90,endpoint=True)
lons = np.linspace(0,360,4,endpoint=False)

folder = '/home/wim/Desktop/Projects/Model/HWM14_MSISE00/indices/'

altitudes = [75,77.5,80,82.5,85,87.5,90,92.5,95,97.5,100,102.5,105,107.5,110,115,120,125,130,135,140,145,150,
             160,170,180,190,200,
             225,250,275,300,325,350,375,400,425,450,475,500,525,550,575,600,625,650,675,700,725,750][::-1]

date_temp = datetime(2014,1,1,1)

day_s = 0
n_days = 365

day_e = day_s + n_days
date_temp = date_temp + timedelta(days=day_s)

time_steps = 2

output = np.empty((4,n_days,np.shape(altitudes)[0],np.shape(lats)[0]))

for day_step in range(day_s,day_e):
    output_t = np.empty((4,time_steps,np.shape(altitudes)[0],np.shape(lats)[0],np.shape(lons)[0]))
    for t_step in range(time_steps):
        for index,item in enumerate(altitudes):
            u,v,t,p = mf.hwm_msise_uvtp_climatology(date_temp,
                                                    item,
                                                    lats,
                                                    lons,
                                                    folder)
            
            output_t[0,t_step,index,:,:] = u
            output_t[1,t_step,index,:,:] = v
            output_t[2,t_step,index,:,:] = t
            output_t[3,t_step,index,:,:] = p
            
        date_temp = date_temp + timedelta(hours=12)
        
    date_save = date_temp + timedelta(days=-1)
    
    # daily means
    output_t = np.mean(output_t,axis=1)
    
    # zonal means 
    output_t = np.mean(output_t,axis=3)
    
    # store in climatology array
    output[0,day_step,:,:] = output_t[0,:,:]
    output[1,day_step,:,:] = output_t[1,:,:]
    output[2,day_step,:,:] = output_t[2,:,:]
    output[3,day_step,:,:] = output_t[3,:,:]
   
    print('One day written')

np.save('/home/wim/Desktop/Projects/Model/HWM14_MSISE00/HWM_uuhr/' + 'climatology' + '_uuhr' + '.npy',
            output)
    
stop = timeit.default_timer()
print('Time: ', stop - start) # takes about 12 min per day

#%%
