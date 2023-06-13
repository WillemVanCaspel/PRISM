#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  6 15:02:04 2020

@author: wim
"""

import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import circmean
import pylunar as pl
from datetime import datetime, timedelta

#%

file_n = 'lunar_hr.dat'
folder = '/home/wim/Desktop/Projects/iSSW/data/sigprim_fits/'
wind_field = np.load(folder + file_n + '_u.npy')
lats = np.load(folder + file_n + '_lat.npy')
lons = np.load(folder + file_n + '_lon.npy')
heights = np.load(folder + file_n + '_height_array.npy')

#%%

def migrating_m2_fft_lat_alt(data,
                             n_days):
    
    ww = n_days * 8
    dt, dz, dy, dx = np.shape(data)
        
    amps = np.empty((dt - ww,dz,dy))
    phas = np.empty((dt - ww,dz,dy))

    for z in range(dz):
        for y in range(dy):
            winds = data[:,z,y,:]
                
            for t in range(dt - ww):
                wind_slice = winds[t:(t + ww),:]
                transform = np.fft.fft2(wind_slice)
               
                wave = transform[n_days * 2,2]
                amps[t,z,y] = np.abs(wave) / (ww) / np.shape(winds)[1] * 2
                phas[t,z,y] = circmean(np.angle(wave))
                
    return amps, phas


def m2_fft_lat_lon(data,
                   n_days):
    
    ww = n_days * 8
    dt, dy, dx = np.shape(data)
        
    amps = np.empty((dt - ww,dy,dx))
    phas = np.empty((dt - ww,dy,dx))
    
    for x in range(dx):
        for y in range(dy):
            winds = data[:,y,x]
                
            for t in range(dt - ww):
                wind_slice = winds[t:(t + ww)]
                transform = np.fft.fft(wind_slice)
               
                wave = transform[n_days * 2]
                amps[t,y,x] = np.abs(wave) / (ww) * 2
                phas[t,y,x] = circmean(np.angle(wave))
                
    return amps, phas


def migrating_m2_fft_lat_time(data,
                              n_days):
    
    ww = n_days * 8
    dt, dy, dx = np.shape(data)
        
    amps = np.empty((dt - ww,dy))
    phas = np.empty((dt - ww,dy))

    for y in range(dy):
        winds = data[:,y,:]
                
        for t in range(dt - ww):
            wind_slice = winds[t:(t + ww),:]
            transform = np.fft.fft2(wind_slice)
           
            wave = transform[n_days * 2,2]
            amps[t,y] = np.abs(wave) / (ww) / np.shape(winds)[1] * 2
            phas[t,y] = circmean(np.angle(wave))
                
    return amps, phas
                    
#%%

n_day_fit = 1
day_s = 0
day_e = 31
a, p = migrating_m2_fft_lat_alt(wind_field[(day_s * 8):(day_e * 8),:,:,:],n_day_fit)

lat_alt_amp = np.mean(a,axis=0)

n_day_fit = 1
day_s = 180
day_e = 211
a, p = migrating_m2_fft_lat_alt(wind_field[(day_s * 8):(day_e * 8),:,:,:],n_day_fit)

lat_lon_amp = np.mean(a,axis=0)

#%


day_s = 0
day_e = 365

alt = 97000
alt_ind = np.argmin(np.abs(heights - alt))

n_day_fit = 5
a, p = migrating_m2_fft_lat_time(wind_field[(day_s * 8):(day_e * 8),alt_ind,:,:],n_day_fit)

p = p[::8,:] # 00:00 solar time every day
p_lunar = np.zeros((np.shape(p)))

dt, dy = np.shape(p_lunar)

day_0 = datetime(2015,1,1)

for t in range(dt):
    year = day_0.year
    month = day_0.month
    day_s = day_0.day
    
    phase = pl.MoonInfo((0,0,0),(0,0,0))
    phase.update((year,month,day_s))
    
    lun_phase = phase.age() / 29.53059 * 2 * np.pi
    
    for y in range(dy):
        p_lunar[t,y] = lun_phase * 2
        
    day_0 = day_0 + timedelta(days=1)
    
p_adjust = np.zeros(np.shape(p))
dt, dy = np.shape(p)
for t in range(dt):
    for y in range(dy):
        p_adjust[t,y] = circmean(p[t,y] + p_lunar[t,y])

lat_time_amp = np.copy(a)
lat_time_pha = np.copy(p_adjust)

#%%

smin = 0
smax = 35

smin_2 = 0
smax_2 = 20

nstep = 3
nstep_2 = 1.

plot_max = 50
plot_max_2 = 15

linefs = 1.0
l_fs = 22.5
contour_step = 1

lat_ticks = np.arange(-75,76,25)
lon_ticks = np.arange(30,351,60)
time_ticks = np.arange(0,361,30)
ts = 25
fs = 28
title_fs = 28

h_min = 80

dy, dt = np.shape(lat_time_amp.T)
amp_time = np.arange(dt) / 8

dy, dt = np.shape(lat_time_pha.T)
pha_time = np.arange(dt)

fig = plt.figure(figsize=(10.,10.))
ax1 = fig.add_axes([0.00, 0.65, 0.75, 0.5]) # [xmin, ymin, dx, dy]
ax2 = fig.add_axes([0.00, 0.05, 1.70, 0.45])
ax3 = fig.add_axes([0.95, 0.65, 0.75, 0.5])

ax1.contourf(lats,heights / 1000,lat_alt_amp[::1,:],np.arange(smin,smax,nstep),vmin=0,vmax=plot_max,cmap='Greys')
cs = ax1.contour(lats,heights / 1000,lat_alt_amp[::1,:], np.arange(smin,smax,nstep),linewidths=linefs,colors='black')
ax1.clabel(cs, np.arange(smin,smax,nstep)[::contour_step], inline=1, fontsize=l_fs,fmt='%1.0f')
ax1.set_xticks(lat_ticks)
ax1.set_xticklabels(lat_ticks,size=ts)
ax1.set_ylim(h_min)
ax1.tick_params(axis='y', which='major', labelsize=ts)
ax1.set_ylabel('Altitude (km)',fontsize=fs)
ax1.set_xlabel(r'Latitude ($^{\circ})$',fontsize=fs)
ax1.set_title(r'Jan migrating LSDT Amp (ms$^{-1}$) in U',fontsize=title_fs)

ax2.contourf(amp_time,lats,lat_time_amp.T,np.arange(smin_2,smax_2,nstep_2),vmin=0,vmax=plot_max_2,cmap='Greys')
cs = ax2.contour(amp_time,lats,lat_time_amp.T, np.arange(smin_2,smax_2,nstep_2),linewidths=linefs,colors='black')
ax2.clabel(cs, np.arange(smin_2,smax_2,nstep_2)[::contour_step], inline=1, fontsize=l_fs,fmt='%1.0f')
ax2.set_xticks(time_ticks)
ax2.set_xticklabels(time_ticks,size=ts)
ax2.set_yticks(lat_ticks)
ax2.set_yticklabels(lat_ticks,size=ts)
ax2.set_ylabel(r'Latitude ($^{\circ})$',fontsize=fs)
ax2.set_xlabel(r'Day of year',fontsize=fs)
ax2.set_title(r'Migrating LSDT Amp (ms$^{-1}$) at 97.5 km in U',fontsize=title_fs)

ax3.contourf(lats,heights / 1000,lat_lon_amp[::1,:],np.arange(smin,smax,nstep),vmin=0,vmax=plot_max,cmap='Greys')
cs = ax3.contour(lats,heights / 1000,lat_lon_amp[::1,:], np.arange(smin,smax,nstep),linewidths=linefs,colors='black')
ax3.clabel(cs, np.arange(smin,smax,nstep)[::contour_step], inline=1, fontsize=l_fs,fmt='%1.0f')
ax3.set_xticks(lat_ticks)
ax3.set_xticklabels(lat_ticks,size=ts)
ax3.set_ylim(h_min)
ax3.tick_params(axis='y', which='major', labelsize=ts)
ax3.set_ylabel('Altitude (km)',fontsize=fs)
ax3.set_xlabel(r'Latitude ($^{\circ})$',fontsize=fs)
ax3.set_title(r'Jul migrating LSDT Amp (ms$^{-1}$) in U',fontsize=title_fs)

fw = 'bold'
ax1.text(-82.5,82.5,'(a)',fontsize=fs + 5,fontweight=fw)
ax3.text(-82.5,82.5,'(b)',fontsize=fs + 5,fontweight=fw)
ax2.text(0,95,'(c)',fontsize=fs + 5,fontweight=fw)


save_f = '/home/wim/Desktop/Projects/iSSW/paper/figures/'
# plt.savefig(save_f + 'lunar_validation.png',dpi=300,bbox_inches="tight")

#%%

x = 1
y = -1
print(2 * np.arctan(y/(np.sqrt(x**2 + y **2) + x)))

#%%

a = 20
phi_1 = 4/5 * np.pi
phi_2 = 3/4 * np.pi

print(2 * a * np.cos(phi_1 / 2))
print(2 * a * np.cos(phi_2 / 2))










