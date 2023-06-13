#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  6 15:02:04 2020

@author: wim
"""

import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import circmean
from netCDF4 import Dataset
from wrf import interplevel
from datetime import datetime, timedelta

#%%

file_n_1 = 'test62.dat'

save_f = '/home/willemvc/Desktop/NTNU files/dope transfer MET/iSSW/data/sigprim_fits/'
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

amps = np.empty((dt - ww,np.shape(altitudes)[0],np.shape(lats)[0],np.shape(waves)[0]))
phas = np.empty((dt - ww,np.shape(altitudes)[0],np.shape(lats)[0],np.shape(waves)[0]))

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
                
#%

fil = 'long_run_11.25.2019.cam.h0.2013-01-05-00000.nc'
fil = 'long_run_11.25.2019.cam.h0.2013-01-12-00000.nc'
fol = '/media/willemvc/NieuwVolume/SD-WACCMX Ozone/'
fol = '/home/willemvc/Desktop/NTNU files/dope transfer MET/iSSW/data/ozone waccm/'
fil = 'waccm_ozone.nc'
wacc = Dataset(fol + fil)

t_ind = 0
alt = [40000]
time = datetime(2000,1,1) + timedelta(days = int(wacc['time'][t_ind]))

ozone = wacc['O3'][t_ind,:,:,:]
geoph = wacc['Z3'][t_ind,:,:,:] * 0.98 # to geometric height
lat = wacc['lat'][:]
lon = wacc['lon'][:]

ozone_field = interplevel(ozone,
                          geoph,
                          alt[0])[:,:]


#%%

# smin = 0
# smax = 5
# nstep = 0.5

# smino = 0
# smaxo = 8.5
# nstepo = 0.5

# smax2 = 19
# nstep2 = 1

# smin3 = 0
# smax3 = 4.5
# nstep3 = 0.5
# plot_max3 = 3.45

# plot_max = 4.1
# plot_max2 = 18
# plot_maxo = 8


# ozone plot
smino = 2
smaxo = 8.5
nstepo = 0.5

# SW2
smin2 = 0
smax2 = 19.6 + 5
nstep2 = 1
plot_max2 = 27

tstart = 0

# SW1 + SW3
smin = smin2
smax = smax2
nstep = 1
plot_max = plot_max2

smin3 = smin
smax3 = smax
nstep3 = nstep
plot_max3 = plot_max


plot_maxo = 8
plot_mino = 1.25

linefs = 0.3
l_fs = 14.
contour_step = 1
contour_step2 = 2

tmin = 20
tmax = 65

zmin = 20
zmax = 80

zmino = -20
zmaxo = 80

ts = 17.
fs = 20.
lw_ssw = 2.0

dx, dy = np.shape(amps[:,0,:,0].T)
dt = np.arange(dy) / 24 + int(ndays / 2)
dy = np.arange(dx)

wave1_ind = 0
wave2_ind = 1
wave3_ind = 2

# #%%

# plt.pcolormesh(amps[:,0,:,wave2_ind])
# plt.colorbar()

# #%%

cm = 'viridis'
cm2 = 'Greys'

fig = plt.figure(figsize=(6.5,6.5))

ax1 = fig.add_axes([0.0, 1.35, 1.0, 0.5]) # [xmin, ymin, dx, dy]
ax2 = fig.add_axes([0.0, 0.65, 1.0, 0.5]) 

ax4 = fig.add_axes([1.4, 1.35, 1.0, 0.5])
ax5 = fig.add_axes([1.4, 0.65, 1.0, 0.5])

im = ax1.contourf(lon,lat,ozone_field*1e6,np.arange(smino,smaxo,nstepo),vmin=plot_mino,vmax=plot_maxo,cmap=cm)
cs = ax1.contour(lon,lat,ozone_field*1e6,np.arange(smino,smaxo,nstepo),linewidths=linefs,colors='black')
ax1.clabel(cs, np.arange(smino,smaxo,nstepo)[::contour_step], inline=1, fontsize=l_fs,fmt='%1.1f')
ax1.set_ylim(zmino,zmaxo)
ax1.tick_params(axis='x', which='major', labelsize=ts)
ax1.tick_params(axis='y', which='major', labelsize=ts)
ax1.set_ylabel(r'Latitude ($^{\circ}$)',fontsize=fs)
ax1.set_xlabel(r'Longitude ($^{\circ}$)',fontsize=fs)
ax1.set_title(r'O$_3$ at 40 km altitude on Jan 10th 2013',fontsize=fs)

im2 = ax2.contourf(dt,lats,amps[:,0,::1,wave1_ind].T,np.arange(smin,smax,nstep),vmin=0,vmax=plot_max2,cmap=cm2)
cs = ax2.contour(dt,lats,amps[:,0,::1,wave1_ind].T, np.arange(smin,smax,nstep),linewidths=linefs,colors='black')
ax2.clabel(cs, np.arange(smin,smax,nstep)[tstart::contour_step], inline=1, fontsize=l_fs,fmt='%1.0f')
ax2.set_ylim(zmin,zmax)
ax2.set_xlim(tmin,tmax)
ax2.tick_params(axis='x', which='major', labelsize=ts)
ax2.tick_params(axis='y', which='major', labelsize=ts)
ax2.set_ylabel(r'Latitude ($^{\circ}$)',fontsize=fs)
ax2.set_title('WACStrat SW1 at 97 km in U',fontsize=fs)

im3 = ax4.contourf(dt,lats,amps[:,0,::1,wave2_ind].T,np.arange(smin2,smax2,nstep2),vmin=0,vmax=plot_max2,cmap=cm2)
cs = ax4.contour(dt,lats,amps[:,0,::1,wave2_ind].T, np.arange(smin2,smax2,nstep2),linewidths=linefs,colors='black')
ax4.clabel(cs, np.arange(smin2,smax2,nstep2)[tstart::contour_step], inline=1, fontsize=l_fs,fmt='%1.0f')
ax4.set_ylim(zmin,zmax)
ax4.set_xlim(tmin,tmax)
ax4.tick_params(axis='x', which='major', labelsize=ts)
ax4.tick_params(axis='y', which='major', labelsize=ts)
ax4.set_ylabel(r'Latitude ($^{\circ}$)',fontsize=fs)
ax4.set_title('WACStrat SW2 at 97 km in U',fontsize=fs)

#%

im4 = ax5.contourf(dt,lats,amps[:,0,::1,wave3_ind].T,np.arange(smin3,smax3,nstep3),vmin=0,vmax=plot_max2,cmap=cm2)
cs = ax5.contour(dt,lats,amps[:,0,::1,wave3_ind].T, np.arange(smin3,smax3,nstep3),linewidths=linefs,colors='black')
ax5.clabel(cs, np.arange(smin3,smax3,nstep3)[tstart::contour_step], inline=1, fontsize=l_fs,fmt='%1.0f')
ax5.set_ylim(zmin,zmax)
ax5.set_xlim(tmin,tmax)
ax5.tick_params(axis='x', which='major', labelsize=ts)
ax5.tick_params(axis='y', which='major', labelsize=ts)
ax5.set_ylabel(r'Latitude ($^{\circ}$)',fontsize=fs)
ax5.set_title('WACStrat SW3 at 97 km in U',fontsize=fs)

ssw_t = 41
ssw_t2 = 34.5
ssw_t3 = 53

lwp = 2.0

cmbar = 'black'
ax2.axvline(ssw_t,0,100,color=cmbar,ls='--',lw=lwp)
ax4.axvline(ssw_t,0,100,color=cmbar,ls='--',lw=lwp)
ax5.axvline(ssw_t,0,100,color=cmbar,ls='--',lw=lwp)

ax2.axvline(ssw_t2,0,100,color=cmbar,ls='--',lw=lwp)
ax4.axvline(ssw_t2,0,100,color=cmbar,ls='--',lw=lwp)
ax5.axvline(ssw_t2,0,100,color=cmbar,ls='--',lw=lwp)

ax2.axvline(ssw_t3,0,100,color=cmbar,ls='--',lw=lwp)
ax4.axvline(ssw_t3,0,100,color=cmbar,ls='--',lw=lwp)
ax5.axvline(ssw_t3,0,100,color=cmbar,ls='--',lw=lwp)

cbar_ax = fig.add_axes([1.05, 1.30, 0.05, 0.5])
cb = fig.colorbar(im, cax=cbar_ax) # ,pad=0.01,fraction=0.1
cb.set_label(label=r'O$_3$ ' + '[ppm]', size=fs + 2)
cb.ax.tick_params(labelsize=ts + 2)

cbar_ax = fig.add_axes([1.05, 0.65, 0.05, 0.5])
cb = fig.colorbar(im2, cax=cbar_ax) # ,pad=0.01,fraction=0.1
cb.set_label(label=r'Amplitude (ms$^{-1}$)', size=fs + 2)
cb.ax.tick_params(labelsize=ts + 2)

cbar_ax = fig.add_axes([2.45, 1.30, 0.05, 0.5])
cb = fig.colorbar(im3, cax=cbar_ax) # ,pad=0.01,fraction=0.1
cb.set_label(label=r'Amplitude (ms$^{-1}$)', size=fs + 2)
cb.ax.tick_params(labelsize=ts + 2)

cbar_ax = fig.add_axes([2.45, 0.65, 0.05, 0.5])
cb = fig.colorbar(im4, cax=cbar_ax) # ,pad=0.01,fraction=0.1
cb.set_label(label=r'Amplitude (ms$^{-1}$)', size=fs + 2)
cb.ax.tick_params(labelsize=ts + 2)

ax2.set_xlabel('Days since Dec 1st 2012',fontsize=fs)
ax5.set_xlabel('Days since Dec 1st 2012',fontsize=fs)

ssw_c = '0.675'

ax2.hlines(43.3,0,100,color=ssw_c,ls='--',lw=lw_ssw)
ax2.hlines(51.3,0,100,color=ssw_c,ls='--',lw=lw_ssw)
ax2.hlines(67.9,0,100,color=ssw_c,ls='--',lw=lw_ssw)

ax4.hlines(43.3,0,100,color=ssw_c,ls='--',lw=lw_ssw)
ax4.hlines(51.3,0,100,color=ssw_c,ls='--',lw=lw_ssw)
ax4.hlines(67.9,0,100,color=ssw_c,ls='--',lw=lw_ssw)

ax5.hlines(43.3,0,100,color=ssw_c,ls='--',lw=lw_ssw)
ax5.hlines(51.3,0,100,color=ssw_c,ls='--',lw=lw_ssw)
ax5.hlines(67.9,0,100,color=ssw_c,ls='--',lw=lw_ssw)


fw = 'bold'
ax1.text(-20,87.,'(a)',fontsize=fs + 5,fontweight=fw)
ax2.text(17.5,84.,'(b)',fontsize=fs + 5,fontweight=fw)
ax4.text(17.5,84.,'(c)',fontsize=fs + 5,fontweight=fw)
ax5.text(17.5,84.,'(d)',fontsize=fs + 5,fontweight=fw)

save_f = '/home/willemvc/Desktop/NTNU files/dope transfer MET/iSSW/paper/MET figures/'
# plt.savefig(save_f + 'StratForcing_decomp_HR.png',dpi=300, bbox_inches="tight")

#%%

print(altitudes)