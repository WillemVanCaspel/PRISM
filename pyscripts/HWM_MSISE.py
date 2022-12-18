#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 19:43:09 2020

@author: wim
"""

from datetime import datetime, timedelta
import numpy as np
import timeit

# these are the HWM14 and MSISE00 Python libraries; install first
from pyhwm2014 import HWM142D
from nrlmsise00 import msise_flat

#%% define function that returns u, v, T, and p fields

def hwm_msise_uvtp(date,
                   alt,
                   lats,
                   lons,
                   folder,
                   use_clim=False):
    
    """
    Function reads f107a, f107 and ap from file based on DOY. Calculates
    pressure from MSISE number densities, temperatures, and the ideal gas law. 
    
    Input:
        - date (datetime object; year, month, day)
        - altitude (km)
        - latitude array (degrees)  (e.g. np.arange(-90,91))
        - longitude array (degrees) (e.g. np.arange(0,360))
        - folder; path to .txt files containing geophysical indices 
    
    Output:
        - lat/lon NRLMSISE-00 temperature and pressure 
        - lat/lon HWM14 zonal and meridional winds
        
    f107a: float
		The observed f107a (81-day running mean of f107) centred at date.
	f107: float
		The observed f107 value on the previous day.
	ap: float
		The ap value at date.
        
    NAVGEM-HA top = 0.1435314 Pa (~ 90 km)
    DOPE top      = 3.57821e-6 Pa
    
    n_levels DOPE between NAVGEM-HA and DOPE top: 32
    """
    
    doy = (date - datetime(2008,1,1,0)).days
    
    # read in f10.7 from file based on input date
    if alt <= 80:   # the models use climatological values below this altitude
        f107a, f107, ap = (150,150,4)
    elif use_clim: 
        f107a, f107, ap = (150,150,4)
    else:
        # f10.7 and (derived) f10.7A
        date_f = float(date.strftime('%Y%m%d'))
        file = np.loadtxt(folder + 'f107_08_16.txt')
        dates, vals = file[:,0], file[:,1]
                
        ind = 0
        found = True
        while found:
            if dates[ind] == date_f:
                found = False
            else:
                ind += 1
                
        f107 = vals[int(ind - 1)] # day before today / doy
        f107a = np.mean(vals[int(ind - 40):int(ind + 41)]) # 81 day mean
                
        # ap index
        file_ap = open(folder + 'ap_08_16.txt','r')

        for m in range(doy):
            line = file_ap.readline()
            ap = float(line[53:55])

        file_ap.close()
            
    # index 5 = density, index -1 = temperature
    msise = msise_flat(date,
                       alt,
                       lats[:,None],
                       lons[None,:],
                       f107a,
                       f107,
                       ap)
        
    r = 8.314               # J / K / mol
    avg = 6.022140e23       # avocado's number
    
    # number density m**-3
    he = 1e6 * msise[:,:,0] 
    o  = 1e6 * msise[:,:,1]  
    n2 = 1e6 * msise[:,:,2] 
    o2 = 1e6 * msise[:,:,3] 
    ar = 1e6 * msise[:,:,4] 
    h  = 1e6 * msise[:,:,6]  
    n  = 1e6 * msise[:,:,7] 
    ox = 1e6 * msise[:,:,8] 
    
    n_density = he + o + n2 + o2 + ar + h + n + ox
    
    temp = msise[:,:,-1]    # temperature in [K]
    
    p = n_density / avg * r * temp # pressure using ideal gas law
        
    # horizontal winds fro HWM14. Does not use f107 and f107a
    winds = HWM142D(alt=alt, glatlim=[lats[0], lats[-1]],
                    glatstp=(lats[1] - lats[0]), glonlim=[lons[0], lons[-1]], 
                    glonstp=(lons[1] - lons[0]), option=6,verbose=False,
                    day=int(date.strftime("%j")), 
                    ut=int(date.strftime('%H')),
                    year=int(date.strftime('%Y')),
                    ap=[-1,ap])
    
    u = winds.Uwind
    v = winds.Vwind
    
    return u, v, temp, p

#%% construct atmosphere and save

# get runtime statistics of this script
start = timeit.default_timer()

lats = np.linspace(-89,89,90,endpoint=True)
lons = np.linspace(0,360,90,endpoint=False)

# output date label
date = '2013' 

# date range of constructed file
date_start = datetime(2012,12,1,0)
date_end = datetime(2013,4,1,0)

length = (date_end - date_start).days

#%

# path to geomagnetic indices used by HWM14/MSISE00
folder = '/home/willemvc/Desktop/projects/PRISM/PRISM/pyscripts/geomagnetic_indices/'

# desired altitudes (km) in output files
altitudes = [75,77.5,80,82.5,85,87.5,90,92.5,95,97.5,100,102.5,105,107.5,110,115,120,125,130,135,140,145,150,
             160,170,180,190,200,
             225,250,275,300,325,350,375,400,425,450,475,500,525,550,575,600,625,650,675,700,725,750][::-1]

# number of timesteps in a day (4 = 6 hrly)
time_steps = 4
date_temp = date_start

output = np.empty((4,
                   length,np.shape(altitudes)[0],
                   np.shape(lats)[0],
                   np.shape(lons)[0]))

for day_step in range(length):
    
    # temporary output
    output_t = np.empty((4,
                         time_steps,
                         np.shape(altitudes)[0],
                         np.shape(lats)[0],
                         np.shape(lons)[0]))
    
    for t_step in range(time_steps):
        for index,item in enumerate(altitudes):
            u,v,t,p = hwm_msise_uvtp(date_temp,
                                     item,
                                     lats,
                                     lons,
                                     folder)
            
            output_t[0,t_step,index,:,:] = u
            output_t[1,t_step,index,:,:] = v
            output_t[2,t_step,index,:,:] = t
            output_t[3,t_step,index,:,:] = p
            
        date_temp = date_temp + timedelta(hours=6)
            
    # daily means
    output_t = np.mean(output_t,axis=1)
        
    # store in climatology array
    output[0,day_step,:,:,:] = output_t[0,:,:,:]
    output[1,day_step,:,:,:] = output_t[1,:,:,:]
    output[2,day_step,:,:,:] = output_t[2,:,:,:]
    output[3,day_step,:,:,:] = output_t[3,:,:,:]
   
    print('One day written')

# save to file for later use
np.save('/home/wim/Desktop/Projects/SSW/data/hwm_msise/' + 'climatology' + '_' + date + '.npy',
            output)
    
stop = timeit.default_timer()
print('Time: ', stop - start) # takes about 12 min per day

#%%
