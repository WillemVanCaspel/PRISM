#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  6 15:02:04 2020

@author: wim
"""

import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import circmean

file_n_1 = 'test25.dat'
file_n_2 = 'test55.dat'
file_n_3 = 'test56.dat'

file_n_2 = 'test55.dat'
file_n_3 = 'test56.dat'

save_f = '/home/wim/Desktop/Projects/iSSW/data/sigprim_fits/'
model_dt = 1

altitudes = [97500]
lat_min = -90
lat_max = 90

ndays = 4 # sliding window length
tide = 2  # 1,2,3 for diurnal, semidiurnal, terdiurnal

full_v_1 = np.load(save_f + file_n_1 + '_u.npy')[:,:,:,:] 
full_v_2 = np.load(save_f + file_n_2 + '_u.npy')[:,:,:,:]
full_v_3 = np.load(save_f + file_n_3 + '_u.npy')[:,:,:,:]
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
full_v_2 = full_v_2[:,:,lat_imin:lat_imax,:]
full_v_3 = full_v_3[:,:,lat_imin:lat_imax,:]

lats = dope_lats[lat_imin:lat_imax]

ww = int(ndays * 24 / model_dt)
skim_step = int(24 / model_dt)

waves = np.arange(-5,6) # zonal wavenumber 5 east to west

amps = np.empty((dt - ww,np.shape(altitudes)[0],np.shape(lats)[0],np.shape(waves)[0]))
phas = np.empty((dt - ww,np.shape(altitudes)[0],np.shape(lats)[0],np.shape(waves)[0]))

amps2 = np.empty((dt - ww,np.shape(altitudes)[0],np.shape(lats)[0],np.shape(waves)[0]))
phas2 = np.empty((dt - ww,np.shape(altitudes)[0],np.shape(lats)[0],np.shape(waves)[0]))

amps3 = np.empty((dt - ww,np.shape(altitudes)[0],np.shape(lats)[0],np.shape(waves)[0]))
phas3 = np.empty((dt - ww,np.shape(altitudes)[0],np.shape(lats)[0],np.shape(waves)[0]))

for index, item in enumerate(alt_inds):
    
    for lat_i, lat_v in enumerate(lats):
        winds = full_v_1[:,int(item),lat_i,:]
        winds2 = full_v_2[:,int(item),lat_i,:]
        winds3 = full_v_3[:,int(item),lat_i,:]
            
        for t in range(dt - ww):
            wind_slice = winds[t:(t + ww),:]
            transform = np.fft.fft2(wind_slice)
            
            wind_slice2 = winds2[t:(t + ww),:]
            transform2 = np.fft.fft2(wind_slice2)
            
            wind_slice3 = winds3[t:(t + ww),:]
            transform3 = np.fft.fft2(wind_slice3)
           
            for wave_ind, wave_number in enumerate(waves):
                wave = transform[ndays * tide,wave_number]
                amps[t,index,lat_i,wave_ind] = np.abs(wave) / (ww) / np.shape(winds)[1] * 2
                phas[t,index,lat_i,wave_ind] = circmean(np.angle(wave))
                
                wave = transform2[ndays * tide,wave_number]
                amps2[t,index,lat_i,wave_ind] = np.abs(wave) / (ww) / np.shape(winds)[1] * 2
                phas2[t,index,lat_i,wave_ind] = circmean(np.angle(wave))
                
                wave = transform3[ndays * tide,wave_number]
                amps3[t,index,lat_i,wave_ind] = np.abs(wave) / (ww) / np.shape(winds)[1] * 2
                phas3[t,index,lat_i,wave_ind] = circmean(np.angle(wave))
                
# #%%

# lat_step = 1

# dx, dy = np.shape(amps[:,0,::1,wave_ind].T)

# dt = np.arange(dy)
# dy = np.arange(dx)

# plt.contourf(dt / 24,lats,amps[:,0,::1,wave_ind].T,vmin=0,levels=10)
# plt.colorbar()
# plt.title('wavenumber ' + str(waves[wave_ind]))
# # plt.title('SW1 amplitude')
# plt.ylabel('latitude')
# plt.xlabel('Days after 1st Dec 2012')
# # plt.yticks(np.arange(np.shape(lats)[0])[::lat_step],np.round(lats[::lat_step],1),size=7)
# plt.xlim(20,65)
# plt.ylim(20,80)
# # plt.ylim(30,80)
# print(file_n_1)

# save_f = '/home/wim/Desktop/'
# # plt.savefig(save_f + 'SW1_2013_ssw.png',dpi=300, bbox_inches="tight")

# #%%

#%%

smin = 0
smax = 56
nstep = 5

smax2 = 31
nstep2 = 3
plot_max2 = 26.5

# smax2 = 14
# nstep2 = 2
# plot_max2 = 14.5

plot_max = 51

linefs = 0.3
l_fs = 14.
contour_step = 1

tmin = 20
tmax = 65

zmin = 30
zmax = 80

ts = 17.
fs = 20.
lw_ssw = 3.0

dx, dy = np.shape(amps[:,0,:,0].T)
dt = np.arange(dy) / 24
dy = np.arange(dx)

wave1_ind = 7
wave2_ind = 6

cm = 'viridis'

fig = plt.figure(figsize=(6.5,6.5))

ax1 = fig.add_axes([0.0, 1.30, 1.0, 0.5]) # [xmin, ymin, dx, dy]
ax2 = fig.add_axes([0.0, 0.65, 1.0, 0.5]) 
ax3 = fig.add_axes([0.0, 0.00, 1.0, 0.5])

ax4 = fig.add_axes([1.4, 1.30, 1.0, 0.5])
ax5 = fig.add_axes([1.4, 0.65, 1.0, 0.5])
ax6 = fig.add_axes([1.4, 0.00, 1.0, 0.5])

ax1.contourf(dt,lats,amps[:,0,::1,wave1_ind].T,np.arange(smin,smax,nstep),vmin=0,vmax=plot_max,cmap=cm)
cs = ax1.contour(dt,lats,amps[:,0,::1,wave1_ind].T, np.arange(smin,smax,nstep),linewidths=linefs,colors='black')
ax1.clabel(cs, np.arange(smin,smax,nstep)[::contour_step], inline=1, fontsize=l_fs,fmt='%1.0f')
ax1.set_ylim(zmin,zmax)
ax1.set_xlim(tmin,tmax)
ax1.tick_params(axis='x', which='major', labelsize=ts)
ax1.tick_params(axis='y', which='major', labelsize=ts)
ax1.set_ylabel(r'Latitude ($^{\circ}$)',fontsize=fs)
ax1.set_title('OnlyThermal SW2 at ' + "{:.0f}".format(altitudes[0] / 1000) + ' km in U',fontsize=fs)

ax2.contourf(dt,lats,amps2[:,0,::1,wave1_ind].T,np.arange(smin,smax,nstep),vmin=0,vmax=plot_max,cmap=cm)
cs = ax2.contour(dt,lats,amps2[:,0,::1,wave1_ind].T, np.arange(smin,smax,nstep),linewidths=linefs,colors='black')
ax2.clabel(cs, np.arange(smin,smax,nstep)[::contour_step], inline=1, fontsize=l_fs,fmt='%1.0f')
ax2.set_ylim(zmin,zmax)
ax2.set_xlim(tmin,tmax)
ax2.tick_params(axis='x', which='major', labelsize=ts)
ax2.tick_params(axis='y', which='major', labelsize=ts)
ax2.set_ylabel(r'Latitude ($^{\circ}$)',fontsize=fs)
ax2.set_title('FixedForcing SW2 at ' + "{:.0f}".format(altitudes[0] / 1000) + ' km in U',fontsize=fs)

im2 = ax3.contourf(dt,lats,amps3[:,0,::1,wave1_ind].T,np.arange(smin,smax,nstep),vmin=0,vmax=plot_max,cmap=cm)
cs = ax3.contour(dt,lats,amps3[:,0,::1,wave1_ind].T, np.arange(smin,smax,nstep),linewidths=linefs,colors='black')
ax3.clabel(cs, np.arange(smin,smax,nstep)[::contour_step], inline=1, fontsize=l_fs,fmt='%1.0f')
ax3.set_ylim(zmin,zmax)
ax3.set_xlim(tmin,tmax)
ax3.tick_params(axis='x', which='major', labelsize=ts)
ax3.tick_params(axis='y', which='major', labelsize=ts)
ax3.set_ylabel(r'Latitude ($^{\circ}$)',fontsize=fs)
ax3.set_title('FixedAtmos SW2 at ' + "{:.0f}".format(altitudes[0] / 1000) + ' km in U',fontsize=fs)

ax4.contourf(dt,lats,amps[:,0,::1,wave2_ind].T,np.arange(smin,smax2,nstep2),vmin=0,vmax=plot_max2,cmap=cm)
cs = ax4.contour(dt,lats,amps[:,0,::1,wave2_ind].T, np.arange(smin,smax2,nstep2),linewidths=linefs,colors='black')
ax4.clabel(cs, np.arange(smin,smax2,nstep2)[::contour_step], inline=1, fontsize=l_fs,fmt='%1.0f')
ax4.set_ylim(zmin,zmax)
ax4.set_xlim(tmin,tmax)
ax4.tick_params(axis='x', which='major', labelsize=ts)
ax4.tick_params(axis='y', which='major', labelsize=ts)
ax4.set_ylabel(r'Latitude ($^{\circ}$)',fontsize=fs)
ax4.set_title('OnlyThermal SW1 at ' + "{:.0f}".format(altitudes[0] / 1000) + ' km in U',fontsize=fs)

ax5.contourf(dt,lats,amps2[:,0,::1,wave2_ind].T,np.arange(smin,smax2,nstep2),vmin=0,vmax=plot_max2,cmap=cm)
cs = ax5.contour(dt,lats,amps2[:,0,::1,wave2_ind].T, np.arange(smin,smax2,nstep2),linewidths=linefs,colors='black')
ax5.clabel(cs, np.arange(smin,smax2,nstep2)[::contour_step], inline=1, fontsize=l_fs,fmt='%1.0f')
ax5.set_ylim(zmin,zmax)
ax5.set_xlim(tmin,tmax)
ax5.tick_params(axis='x', which='major', labelsize=ts)
ax5.tick_params(axis='y', which='major', labelsize=ts)
ax5.set_ylabel(r'Latitude ($^{\circ}$)',fontsize=fs)
ax5.set_title('FixedForcing SW1 at ' + "{:.0f}".format(altitudes[0] / 1000) + ' km in U',fontsize=fs)

im = ax6.contourf(dt,lats,amps3[:,0,::1,wave2_ind].T,np.arange(smin,smax2,nstep2),vmin=0,vmax=plot_max2,cmap=cm)
cs = ax6.contour(dt,lats,amps3[:,0,::1,wave2_ind].T, np.arange(smin,smax2,nstep2),linewidths=linefs,colors='black')
ax6.clabel(cs, np.arange(smin,smax2,nstep2)[::contour_step], inline=1, fontsize=l_fs,fmt='%1.0f')
ax6.set_ylim(zmin,zmax)
ax6.set_xlim(tmin,tmax)
ax6.tick_params(axis='x', which='major', labelsize=ts)
ax6.tick_params(axis='y', which='major', labelsize=ts)
ax6.set_ylabel(r'Latitude ($^{\circ}$)',fontsize=fs)
ax6.set_title('FixedAtmos SW1 at ' + "{:.0f}".format(altitudes[0] / 1000) + ' km in U',fontsize=fs)

cbar_ax = fig.add_axes([2.45, 0.35, 0.04, 1.10])
cb = fig.colorbar(im, cax=cbar_ax) # ,pad=0.01,fraction=0.1
cb.set_label(label='SW1 amplitude (ms$^{-1}$)', size=fs + 2)
cb.ax.tick_params(labelsize=ts + 2)

cbar_ax = fig.add_axes([1.05, 0.35, 0.04, 1.10])
cb = fig.colorbar(im2, cax=cbar_ax) # ,pad=0.01,fraction=0.1
cb.set_label(label=r'SW2 amplitude (ms$^{-1}$)', size=fs + 2)
cb.ax.tick_params(labelsize=ts + 2)

ax3.set_xlabel('Days since Dec 1st 2012',fontsize=fs)
ax6.set_xlabel('Days since Dec 1st 2012',fontsize=fs)

save_f = '/home/wim/Desktop/Projects/iSSW/paper/figures/'
# plt.savefig(save_f + 'forcing_decomp.png',dpi=300, bbox_inches="tight")

#%%

print(altitudes)