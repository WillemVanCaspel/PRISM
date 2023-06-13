#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  6 15:02:04 2020

@author: wim
"""

import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import circmean

file_n_1 = 'test1.dat'

save_f = '/home/wim/Desktop/Projects/iSSW/data/sigprim_fits/'
model_dt = 1

altitudes = [97500]
lat_min = -90
lat_max = 90

ndays = 4 # sliding window length
tide = 2  # 1,2,3 for diurnal, semidiurnal, terdiurnal

full_v_1 = np.load(save_f + file_n_1 + '_u.npy')[:,:,:,:] 
height_min = 0
height_max = 100

dt, dz, dy, dx = np.shape(full_v_1)

height_arr = np.load(save_f + file_n_1 + '_height_array.npy')[height_min:height_max]
dope_lats = np.load(save_f + file_n_1 + '_lat.npy')[:]
dope_lons = np.load(save_f + file_n_1 + '_lon.npy')[:]

alt_inds = np.zeros(np.shape(altitudes)[0])
for index, item in enumerate(altitudes):
    alt_inds[index] = np.argmin(np.abs(height_arr - item))
    
lat_imin = np.argmin(np.abs(dope_lats - lat_min))
lat_imax = np.argmin(np.abs(dope_lats - lat_max))

full_v_1 = full_v_1[:,:,lat_imin:lat_imax,:]

lats = dope_lats[lat_imin:lat_imax]

ww = int(ndays * 24 / model_dt)
skim_step = int(24 / model_dt)

waves = np.arange(1,4) # zonal wavenumber 5 east to west
# waves = np.array([1,-4,-5])


amps = np.empty((dt - ww,np.shape(altitudes)[0],np.shape(lats)[0],np.shape(waves)[0]))
phas = np.empty((dt - ww,np.shape(altitudes)[0],np.shape(lats)[0],np.shape(waves)[0]))

amps2 = np.empty((dt - ww,np.shape(altitudes)[0],np.shape(lats)[0],np.shape(waves)[0]))
phas2 = np.empty((dt - ww,np.shape(altitudes)[0],np.shape(lats)[0],np.shape(waves)[0]))

amps3 = np.empty((dt - ww,np.shape(altitudes)[0],np.shape(lats)[0],np.shape(waves)[0]))
phas3 = np.empty((dt - ww,np.shape(altitudes)[0],np.shape(lats)[0],np.shape(waves)[0]))

for index, item in enumerate(alt_inds):
    
    for lat_i, lat_v in enumerate(lats):
        winds = full_v_1[:,int(item),lat_i,:]
            
        for t in range(dt - ww):
            wind_slice = winds[t:(t + ww),:]
            transform = np.fft.fft2(wind_slice)
                       
            for wave_ind, wave_number in enumerate(waves):
                wave = transform[ndays * tide,wave_number]
                amps[t,index,lat_i,wave_ind] = np.abs(wave) / (ww) / np.shape(winds)[1] * 2
                phas[t,index,lat_i,wave_ind] = circmean(np.angle(wave))
                
#%%

smin = 0
smax = 61
nstep = 5
plot_max = 57

smax2 = 28
nstep2 = 3
plot_max2 = 25.5

smax3 = 21
nstep3 = 2
plot_max3 = 19.

linefs = 0.3
l_fs = 14.
contour_step = 1

tmin = 20
tmax = 65

zmin = 20
zmax = 80

ts = 17.
fs = 20.
lw_ssw = 3.0

dx, dy = np.shape(amps[:,0,:,0].T)
dt = np.arange(dy) / 24 + int(ndays / 2)
dy = np.arange(dx)

lw_ssw = 3.0

wave1_ind = 1
wave2_ind = 0
wave3_ind = 2

cm = 'viridis'

fig = plt.figure(figsize=(6.5,6.5))

ax1 = fig.add_axes([0.0, 1.30, 1.0, 0.5]) # [xmin, ymin, dx, dy]
ax2 = fig.add_axes([0.0, 0.65, 1.0, 0.5]) 
ax3 = fig.add_axes([0.0, 0.00, 1.0, 0.5])

name = 'OnlyThermal'

im = ax1.contourf(dt,lats,amps[:,0,::1,wave1_ind].T,np.arange(smin,smax,nstep),vmin=0,vmax=plot_max,cmap=cm)
cs = ax1.contour(dt,lats,amps[:,0,::1,wave1_ind].T, np.arange(smin,smax,nstep),linewidths=linefs,colors='black')
ax1.clabel(cs, np.arange(smin,smax,nstep)[::contour_step], inline=1, fontsize=l_fs,fmt='%1.0f')
ax1.set_ylim(zmin,zmax)
ax1.set_xlim(tmin,tmax)
ax1.tick_params(axis='x', which='major', labelsize=ts)
ax1.tick_params(axis='y', which='major', labelsize=ts)
ax1.set_ylabel(r'Latitude ($^{\circ}$)',fontsize=fs)
ax1.set_title(name + ' SW2 at ' + "{:.1f}".format(altitudes[0] / 1000) + ' km in U',fontsize=fs)

im2 = ax2.contourf(dt,lats,amps[:,0,::1,wave2_ind].T,np.arange(smin,smax2,nstep2),vmin=0,vmax=plot_max2,cmap=cm)
cs = ax2.contour(dt,lats,amps[:,0,::1,wave2_ind].T, np.arange(smin,smax2,nstep2),linewidths=linefs,colors='black')
ax2.clabel(cs, np.arange(smin,smax2,nstep2)[::contour_step], inline=1, fontsize=l_fs,fmt='%1.0f')
ax2.set_ylim(zmin,zmax)
ax2.set_xlim(tmin,tmax)
ax2.tick_params(axis='x', which='major', labelsize=ts)
ax2.tick_params(axis='y', which='major', labelsize=ts)
ax2.set_ylabel(r'Latitude ($^{\circ}$)',fontsize=fs)
ax2.set_title(name + ' SW1 at ' + "{:.1f}".format(altitudes[0] / 1000) + ' km in U',fontsize=fs)

im3 = ax3.contourf(dt,lats,amps[:,0,::1,wave3_ind].T,np.arange(smin,smax3,nstep3),vmin=0,vmax=plot_max3,cmap=cm)
cs = ax3.contour(dt,lats,amps[:,0,::1,wave3_ind].T, np.arange(smin,smax3,nstep3),linewidths=linefs,colors='black')
ax3.clabel(cs, np.arange(smin,smax3,nstep3)[::contour_step], inline=1, fontsize=l_fs,fmt='%1.0f')
ax3.set_ylim(zmin,zmax)
ax3.set_xlim(tmin,tmax)
ax3.tick_params(axis='x', which='major', labelsize=ts)
ax3.tick_params(axis='y', which='major', labelsize=ts)
ax3.set_ylabel(r'Latitude ($^{\circ}$)',fontsize=fs)
ax3.set_title(name + ' SW3 at ' + "{:.1f}".format(altitudes[0] / 1000) + ' km in U',fontsize=fs)

cbar_ax = fig.add_axes([1.05, 1.3, 0.04, 0.5])
cb = fig.colorbar(im, cax=cbar_ax) 
cb.set_label(label='SW2 amplitude (ms$^{-1}$)', size=fs + 2)
cb.ax.tick_params(labelsize=ts + 2)

cbar_ax = fig.add_axes([1.05, 0.65, 0.04, 0.5])
cb = fig.colorbar(im2, cax=cbar_ax)
cb.set_label(label='SW1 amplitude (ms$^{-1}$)', size=fs + 2)
cb.ax.tick_params(labelsize=ts + 2)

cbar_ax = fig.add_axes([1.05, 0.0, 0.04, 0.5])
cb = fig.colorbar(im3, cax=cbar_ax)
cb.set_label(label=r'SW3 amplitude (ms$^{-1}$)', size=fs + 2)
cb.ax.tick_params(labelsize=ts + 2)

ssw_c = 'white'
ssw_l = 41
ax1.vlines(ssw_l,0,100,color=ssw_c,ls='--',lw=lw_ssw) # x, ymin, ymax
ax2.vlines(ssw_l,0,100,color=ssw_c,ls='--',lw=lw_ssw)
ax3.vlines(ssw_l,0,100,color=ssw_c,ls='--',lw=lw_ssw)

fw = 'bold'
ax1.text(17.5,85.,'(d)',fontsize=fs + 5,fontweight=fw)
ax2.text(17.5,85.,'(e)',fontsize=fs + 5,fontweight=fw)
ax3.text(17.5,85.,'(f)',fontsize=fs + 5,fontweight=fw)

ax3.set_xlabel('Days since Dec 1st 2012',fontsize=fs)

save_f = '/home/wim/Desktop/Projects/iSSW/paper/figures/'
# plt.savefig(save_f + 'PRISM_sw1_sw2_sw3_solar.png',dpi=300, bbox_inches="tight")

#%%

print(altitudes)