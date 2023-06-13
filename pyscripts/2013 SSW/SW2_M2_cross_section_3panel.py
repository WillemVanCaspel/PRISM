#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  6 15:02:04 2020

@author: wim
"""

import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import circmean
# import pylunar as pl
from datetime import datetime, timedelta
import scipy.stats

#%

file_n_1 = 'test26.dat' # only lunar
file_n_2 = 'test25.dat' # only thermal
file_n_3 = 'test58.dat' # PRISM
folder = '/home/willemvc/Desktop/NTNU files/dope transfer MET/iSSW/data/sigprim_fits/'

wind_field_1 = np.load(folder + file_n_1 + '_u.npy')
wind_field_2 = np.load(folder + file_n_2 + '_u.npy')
wind_field_3 = np.load(folder + file_n_3 + '_u.npy')
# wind_field_3 = wind_field_1 + wind_field_2
lats = np.load(folder + file_n_1 + '_lat.npy')
lons = np.load(folder + file_n_1 + '_lon.npy')
heights = np.load(folder + file_n_1 + '_height_array.npy')

def migrating_m2_fft_lat_alt(data,
                             n_days):
    
    ww = n_days * 24
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
    
    ww = n_days * 24
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


def migrating_m2_fft_alt_tim(data,
                             n_days,
                             latitudes,
                             lat_loc):
    
    lat_ind = np.argmin(np.abs(latitudes - lat_loc))
    print(lat_ind)
    
    ww = n_days * 24
    dt, dz, dy, dx = np.shape(data)
        
    amps = np.empty((dt - ww,dz))
    phas = np.empty((dt - ww,dz))

    for z in range(dz):
        winds = data[:,z,lat_ind,:]
            
        for t in range(dt - ww):
            wind_slice = winds[t:(t + ww),:]
            transform = np.fft.fft2(wind_slice)
           
            wave = transform[n_days * 2,2]
            amps[t,z] = np.abs(wave) / (ww) / np.shape(winds)[1] * 2
            phas[t,z] = circmean(np.angle(wave))
                
    return amps, phas


def migrating_m2_fft_lat_time(data,
                              n_days):
    
    ww = n_days * 24
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
                    
#%

n_day_fit = 4
day_s = 0
day_e = 80

latitude = 60

a, p = migrating_m2_fft_alt_tim(wind_field_1[(day_s * 24):(day_e * 24),:,:,:],
                                n_day_fit,
                                lats,
                                latitude)

p = p[::24,:] # 00:00 solar time every day
    

alt_time_amp_1 = np.copy(a)
alt_time_pha_1 = np.copy(p)

a, p = migrating_m2_fft_alt_tim(wind_field_2[(day_s * 24):(day_e * 24),:,:,:],
                                n_day_fit,
                                lats,
                                latitude)

p = p[::24,:] # 00:00 solar time every day
    

alt_time_amp_2 = np.copy(a)
alt_time_pha_2 = np.copy(p)

a, p = migrating_m2_fft_alt_tim(wind_field_3[(day_s * 24):(day_e * 24),:,:,:],
                                n_day_fit,
                                lats,
                                latitude)

p = p[::24,:] # 00:00 solar time every day
    

alt_time_amp_3 = np.copy(a)
alt_time_pha_3 = np.copy(p)

def gaussian_weighted_average(altitude,
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
    
    return weights

heights_sd = np.arange(130001,69999,-500)
sd_obs = np.zeros(np.shape(heights_sd))

for index, item in enumerate(heights_sd):
    sd_obs[index] = gaussian_weighted_average(item,
                                              100000,
                                              30000)
    
#%%

smin = 0
smax = 86
nstep = 10

smax2 = smax # 25 normally
nstep2 = 2 # 2 normally

plot_max = 90 # 76 norally
plot_max2 = plot_max # 23 normaly

linefs = 0.3
l_fs = 14.
contour_step = 1

tmin = 20
tmax = 65

zmin = 70
zmax = 130

ts = 17.2
fs = 21.
lw_ssw = 3.0

dt, dz = np.shape(alt_time_amp_2)
t = np.arange(dt) / 24 + int(n_day_fit / 2)

cm = 'Greys'

fig = plt.figure(figsize=(6.5,6.5))

ax3 = fig.add_axes([0.0, 1.30, 1.0, 0.5]) # [xmin, ymin, dx, dy]
ax1 = fig.add_axes([0.0, 0.65, 1.0, 0.5]) 
ax2 = fig.add_axes([0.0, 0.00, 1.0, 0.5])

im2 = ax1.contourf(t,heights / 1000,alt_time_amp_2[:,:].T,np.arange(smin,smax,nstep),vmin=0,vmax=plot_max,cmap=cm)
cs = ax1.contour(t,heights / 1000,alt_time_amp_2[:,:].T, np.arange(smin,smax,nstep),linewidths=linefs,colors='black')
ax1.clabel(cs, np.arange(smin,smax,nstep)[::contour_step], inline=1, fontsize=l_fs,fmt='%1.0f')
ax1.set_ylim(zmin,zmax)
ax1.set_xlim(tmin,tmax)
ax1.tick_params(axis='x', which='major', labelsize=ts)
ax1.tick_params(axis='y', which='major', labelsize=ts)
ax1.set_ylabel(r'Altitude (km)',fontsize=fs)
ax1.set_title(r'OnlySolar SW2 at 60$^{\circ}$N in U',fontsize=fs)
# ax1.set_title(r'Only Strat SW2 at 60$^{\circ}$N in U',fontsize=fs)

im = ax2.contourf(t,heights / 1000,alt_time_amp_1[:,:].T,np.arange(smin,smax2,nstep2),vmin=0,vmax=plot_max2,cmap=cm)
cs = ax2.contour(t,heights / 1000,alt_time_amp_1[:,:].T, np.arange(smin,smax2,nstep2),linewidths=linefs,colors='black')
ax2.clabel(cs, np.arange(smin,smax2,nstep2)[::contour_step], inline=1, fontsize=l_fs,fmt='%1.0f')
ax2.set_ylim(zmin,zmax)
ax2.set_xlim(tmin,tmax)
ax2.tick_params(axis='x', which='major', labelsize=ts)
ax2.tick_params(axis='y', which='major', labelsize=ts)
ax2.set_ylabel(r'Altitude (km)',fontsize=fs)
ax2.set_title(r'OnlyLunar SW2 at 60$^{\circ}$N in U',fontsize=fs)
# ax2.set_title(r'Only Trop SW2 at 60$^{\circ}$N in U',fontsize=fs)
ax2.set_xlabel('Days since Dec 1st 2012',fontsize=fs)

im3 = ax3.contourf(t,heights / 1000,alt_time_amp_3[:,:].T,np.arange(smin,smax,nstep),vmin=0,vmax=plot_max,cmap=cm)
cs = ax3.contour(t,heights / 1000,alt_time_amp_3[:,:].T, np.arange(smin,smax,nstep),linewidths=linefs,colors='black')
ax3.clabel(cs, np.arange(smin,smax,nstep)[::contour_step], inline=1, fontsize=l_fs,fmt='%1.0f')
ax3.set_ylim(zmin,zmax)
ax3.set_xlim(tmin,tmax)
ax3.tick_params(axis='x', which='major', labelsize=ts)
ax3.tick_params(axis='y', which='major', labelsize=ts)
ax3.set_ylabel(r'Altitude (km)',fontsize=fs)
ax3.set_title(r'PRISM SW2 at 60$^{\circ}$N in U',fontsize=fs)
# ax3.set_title(r'Trop + Strat thermal SW2 at 60$^{\circ}$N in U',fontsize=fs)

cbtick = np.arange(smin,smax,10)
# cbtick = np.arange(0,130,10)
cbar_ax = fig.add_axes([1.05, 1.30, 0.04, 0.5])
cb = fig.colorbar(im3, cax=cbar_ax) 
cb.set_ticks(cbtick)
cb.set_ticklabels(cbtick)
cb.set_label(label=r'Amplitude (ms$^{-1}$)', size=fs)
cb.ax.tick_params(labelsize=ts + 2)

cbar_ax = fig.add_axes([1.05, 0.65, 0.04, 0.5])
cb = fig.colorbar(im2, cax=cbar_ax) 
cb.set_ticks(cbtick)
cb.set_ticklabels(cbtick)
cb.set_label(label=r'Amplitude (ms$^{-1}$)', size=fs)
cb.ax.tick_params(labelsize=ts + 2)

cbar_ax = fig.add_axes([1.05, 0.00, 0.04, 0.5])
cb = fig.colorbar(im, cax=cbar_ax) 
cb.set_ticks(cbtick)
cb.set_ticklabels(cbtick)
cb.set_label(label=r'Amplitude (ms$^{-1}$)', size=fs)
cb.ax.tick_params(labelsize=ts + 2)

# ms = 5
# ax1.plot(tmin + sd_obs*2.5e5,heights_sd / 1000,color='black',ls='-')
# ax1.plot(tmin + sd_obs*2.5e5,heights_sd / 1000,'o',color='black',ms=ms,markevery=5) 
# ax2.plot(tmin + sd_obs*2.5e5,heights_sd / 1000,color='black',ls='-')
# ax2.plot(tmin + sd_obs*2.5e5,heights_sd / 1000,'o',color='black',ms=ms,markevery=5) 

ssw_c = 'black'
ssw_l = 41
ssw_l2 = 34.5
ssw_l3 = 53
ax1.vlines(ssw_l,0,150,color=ssw_c,ls='--',lw=lw_ssw)  # x, ymin, ymax
ax2.vlines(ssw_l,0,150,color=ssw_c,ls='--',lw=lw_ssw)
ax3.vlines(ssw_l,0,150,color=ssw_c,ls='--',lw=lw_ssw)
ax1.vlines(ssw_l2,0,150,color=ssw_c,ls='--',lw=lw_ssw) # x, ymin, ymax
ax2.vlines(ssw_l2,0,150,color=ssw_c,ls='--',lw=lw_ssw)
ax3.vlines(ssw_l2,0,150,color=ssw_c,ls='--',lw=lw_ssw)
ax1.vlines(ssw_l3,0,150,color=ssw_c,ls='--',lw=lw_ssw) # x, ymin, ymax
ax2.vlines(ssw_l3,0,150,color=ssw_c,ls='--',lw=lw_ssw)
ax3.vlines(ssw_l3,0,150,color=ssw_c,ls='--',lw=lw_ssw)

fw = 'bold' # text font weight
ax1.text(17.5,133.5,'(b)',fontsize=fs + 4,fontweight=fw)
ax2.text(17.5,133.5,'(c)',fontsize=fs + 4,fontweight=fw)
ax3.text(17.5,133.5,'(a)',fontsize=fs + 4,fontweight=fw)


save_f = '/home/willemvc/Desktop/NTNU files/dope transfer MET/iSSW/paper/MET figures/'
# save_f = '/home/wim/Desktop/'
# plt.savefig(save_f + 'sw2_m2_cross_section_3p_HR.png',dpi=300, bbox_inches="tight")

#%%

perc = alt_time_amp_1[:,:].T / alt_time_amp_2[:,:].T

cutoff = 0.7
inds = np.where(perc > cutoff)
perc[inds] = np.nan

plt.pcolormesh(t,heights / 1000,perc)
plt.colorbar()
# plt.ylim(0,20)

#%%


sw2_lvls = np.arange(0,91,10)
#%

plt.contourf(t,heights,alt_time_amp_2[:,:].T,sw2_lvls,cmap='Greys')
plt.colorbar()
plt.xlim(20,65)

#%%

plt.figure()
plt.contourf(t,heights,alt_time_amp_1[:,:].T,cmap='Greys')
plt.colorbar()
plt.xlim(20,65)

#%%

n_day_fit = 1
day_s = 0
day_e = 31
a, p = migrating_m2_fft_lat_alt(wind_field[(day_s * 24):(day_e * 24),:,:,:],n_day_fit)

lat_alt_amp = np.mean(a,axis=0)

alt = 115000
alt_ind = np.argmin(np.abs(heights - alt))

a, p = m2_fft_lat_lon(wind_field[(day_s * 24):(day_e * 24),alt_ind,:,:],n_day_fit)

lat_lon_amp = np.mean(a,axis=0)

#%%

day_s = 0
day_e = 365

alt = 105000
alt_ind = np.argmin(np.abs(heights - alt))

a, p = migrating_m2_fft_lat_time(wind_field[(day_s * 24):(day_e * 24),alt_ind,:,:],n_day_fit)

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
nstep = 3
nstep_2 = 5

plot_max = 75

linefs = 1.0
l_fs = 17.5
contour_step = 1

lat_ticks = np.arange(-75,76,25)
lon_ticks = np.arange(30,351,60)

ts = 24
fs = 26
title_fs = 24

h_min = 70

# dy, dt = np.shape(lat_time_amp.T)
# amp_time = np.arange(dt) / 8

# dy, dt = np.shape(lat_time_pha.T)
# pha_time = np.arange(dt)

fig = plt.figure(figsize=(10.,10.))

# ax1 = fig.add_axes([0.0, 0.65, 0.75, 0.5]) # [xmin, ymin, dx, dy]
# ax3 = fig.add_axes([0.0, 0.00, 0.75, 0.5])

ax1 = fig.add_axes([0.0, 0.65, 0.75, 0.5]) # [xmin, ymin, dx, dy]
ax3 = fig.add_axes([0.95, 0.65, 0.75, 0.5])

ax1.contourf(lats,heights / 1000,lat_alt_amp[::1,:],np.arange(smin,smax,nstep),vmin=0,vmax=plot_max,cmap='Greys')
cs = ax1.contour(lats,heights / 1000,lat_alt_amp[::1,:], np.arange(smin,smax,nstep),linewidths=linefs,colors='black')
ax1.clabel(cs, np.arange(smin,smax,nstep)[::contour_step], inline=1, fontsize=l_fs,fmt='%1.0f')

ax1.set_xticks(lat_ticks)
ax1.set_xticklabels(lat_ticks,size=ts)

ax1.tick_params(axis='x', which='major', labelsize=ts)
ax1.tick_params(axis='y', which='major', labelsize=ts)
ax1.set_ylabel('Altitude (km)',fontsize=fs)
ax1.set_xlabel(r'Latitude ($^{\circ})$',fontsize=fs)
ax1.set_title('Jan migrating M2 Amp (m/s) in V',fontsize=title_fs)

ax3.contourf(lons,lats,lat_lon_amp,np.arange(smin,smax,nstep_2),vmin=0,vmax=plot_max,cmap='Greys')
cs = ax3.contour(lons,lats,lat_lon_amp, np.arange(smin,smax,nstep_2),linewidths=linefs,colors='black')
ax3.clabel(cs, np.arange(smin,smax,nstep_2)[::contour_step], inline=1, fontsize=l_fs,fmt='%1.0f')
ax3.set_ylim(lats[0],lats[-1])
ax3.set_xlim(lons[0],lons[-1])
ax3.set_xticks(lon_ticks)
ax3.set_xticklabels(lon_ticks,size=ts)
ax3.set_yticks(lat_ticks)
ax3.set_yticklabels(lat_ticks,size=ts)
ax3.tick_params(axis='x', which='major', labelsize=ts)
ax3.tick_params(axis='y', which='major', labelsize=ts)
ax3.set_ylabel(r'Latitude ($^{\circ})$',fontsize=fs)
ax3.set_xlabel(r'Longitude ($^{\circ})$',fontsize=fs)
ax3.set_title('Jan M2 Amp (m/s) at 115 km in V',fontsize=title_fs)

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










