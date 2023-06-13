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
file_n_2 = 'test58.dat'
file_n_3 = 'test25.dat'

file_n_2 = 'test60.dat'
file_n_3 = 'test59.dat'

file_n_2 = 'test65.dat'
file_n_3 = 'test66.dat'

# file_n_2 = 'test2.dat'
# file_n_3 = 'test66.dat'

save_f = '/home/willemvc/Desktop/NTNU files/dope transfer MET/iSSW/data/sigprim_fits/'
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

waves = np.arange(1,7) # zonal wavenumber 5 east to west

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
                
#%%

smin = 0
smax = 71
nstep = 5
plot_max = 67

# smax2 = 28
# nstep2 = 3
# plot_max2 = 25.5

# smax3 = 21
# nstep3 = 2
# plot_max3 = 19.

smax2 = smax
nstep2 = 5
plot_max2 = plot_max

smax3 = smax
nstep3 = 5
plot_max3 = plot_max

linefs = 0.3
l_fs = 14.
contour_step = 1

tmin = 20
tmax = 65



zmin = 20
zmax = 80

ts = 17.
fs = 20.

dx, dy = np.shape(amps[:,0,:,0].T)
dt = np.arange(dy) / 24 + int(ndays / 2)
dy = np.arange(dx)

wave1_ind = 1
wave2_ind = 0
wave3_ind = 2

# wave1_ind = 3
# wave2_ind = 4
# wave3_ind = 5

cm = 'Greys'
# cm = 'viridis'

fig = plt.figure(figsize=(6.75,5.75))

ax2 = fig.add_axes([0.0, 1.30, 0.95, 0.5]) # [xmin, ymin, dx, dy]
ax1 = fig.add_axes([0.0, 0.65, 0.95, 0.5]) 
ax3 = fig.add_axes([0.0, 0.00, 0.95, 0.5])

ax5 = fig.add_axes([1.35, 1.30, 0.95, 0.5])
ax4 = fig.add_axes([1.35, 0.65, 0.95, 0.5])
ax6 = fig.add_axes([1.35, 0.00, 0.95, 0.5])

name = 'FixedAtmos'

im = ax1.contourf(dt,lats,amps2[:,0,::1,wave1_ind].T,np.arange(smin,smax,nstep),vmin=0,vmax=plot_max,cmap=cm)
cs = ax1.contour(dt,lats,amps2[:,0,::1,wave1_ind].T, np.arange(smin,smax,nstep),linewidths=linefs,colors='black')
ax1.clabel(cs, np.arange(smin,smax,nstep)[::contour_step], inline=1, fontsize=l_fs,fmt='%1.0f')
ax1.set_ylim(zmin,zmax)
ax1.set_xlim(tmin,tmax)
ax1.tick_params(axis='x', which='major', labelsize=ts)
ax1.tick_params(axis='y', which='major', labelsize=ts)
ax1.set_ylabel(r'Latitude ($^{\circ}$)',fontsize=fs)
ax1.set_title(name + ' SW2 at ' + "{:.1f}".format(altitudes[0] / 1000) + ' km in U',fontsize=fs)
ax1.set_title(name + ' SW2 at 97 km in U',fontsize=fs)


im2 = ax2.contourf(dt,lats,amps2[:,0,::1,wave2_ind].T,np.arange(smin,smax2,nstep2),vmin=0,vmax=plot_max2,cmap=cm)
cs = ax2.contour(dt,lats,amps2[:,0,::1,wave2_ind].T, np.arange(smin,smax2,nstep2),linewidths=linefs,colors='black')
ax2.clabel(cs, np.arange(smin,smax2,nstep2)[::contour_step], inline=1, fontsize=l_fs,fmt='%1.0f')
ax2.set_ylim(zmin,zmax)
ax2.set_xlim(tmin,tmax)
ax2.tick_params(axis='x', which='major', labelsize=ts)
ax2.tick_params(axis='y', which='major', labelsize=ts)
ax2.set_ylabel(r'Latitude ($^{\circ}$)',fontsize=fs)
ax2.set_title(name + ' SW1 at ' + "{:.1f}".format(altitudes[0] / 1000) + ' km in U',fontsize=fs)
ax2.set_title(name + ' SW1 at 97 km in U',fontsize=fs)

im3 = ax3.contourf(dt,lats,amps2[:,0,::1,wave3_ind].T,np.arange(smin,smax3,nstep3),vmin=0,vmax=plot_max3,cmap=cm)
cs = ax3.contour(dt,lats,amps2[:,0,::1,wave3_ind].T, np.arange(smin,smax3,nstep3),linewidths=linefs,colors='black')
ax3.clabel(cs, np.arange(smin,smax3,nstep3)[::contour_step], inline=1, fontsize=l_fs,fmt='%1.0f')
ax3.set_ylim(zmin,zmax)
ax3.set_xlim(tmin,tmax)
ax3.tick_params(axis='x', which='major', labelsize=ts)
ax3.tick_params(axis='y', which='major', labelsize=ts)
ax3.set_ylabel(r'Latitude ($^{\circ}$)',fontsize=fs)
ax3.set_title(name + ' SW3 at ' + "{:.1f}".format(altitudes[0] / 1000) + ' km in U',fontsize=fs)
ax3.set_title(name + ' SW3 at 97 km in U',fontsize=fs)


name = 'FixedForcingZM'


im4 = ax4.contourf(dt,lats,amps3[:,0,::1,wave1_ind].T,np.arange(smin,smax,nstep),vmin=0,vmax=plot_max,cmap=cm)
cs = ax4.contour(dt,lats,amps3[:,0,::1,wave1_ind].T, np.arange(smin,smax,nstep),linewidths=linefs,colors='black')
ax4.clabel(cs, np.arange(smin,smax,nstep)[::contour_step], inline=1, fontsize=l_fs,fmt='%1.0f')
ax4.set_ylim(zmin,zmax)
ax4.set_xlim(tmin,tmax)
ax4.tick_params(axis='x', which='major', labelsize=ts)
ax4.tick_params(axis='y', which='major', labelsize=ts)
ax4.set_ylabel(r'Latitude ($^{\circ}$)',fontsize=fs)
ax4.set_title(name + ' SW2 at ' + "{:.1f}".format(altitudes[0] / 1000) + ' km in U',fontsize=fs)
ax4.set_title(name + ' SW2 at 97 km in U',fontsize=fs)


im5 = ax5.contourf(dt,lats,amps3[:,0,::1,wave2_ind].T,np.arange(smin,smax2,nstep2),vmin=0,vmax=plot_max2,cmap=cm)
cs = ax5.contour(dt,lats,amps3[:,0,::1,wave2_ind].T, np.arange(smin,smax2,nstep2),linewidths=linefs,colors='black')
ax5.clabel(cs, np.arange(smin,smax2,nstep2)[::contour_step], inline=1, fontsize=l_fs,fmt='%1.0f')
ax5.set_ylim(zmin,zmax)
ax5.set_xlim(tmin,tmax)
ax5.tick_params(axis='x', which='major', labelsize=ts)
ax5.tick_params(axis='y', which='major', labelsize=ts)
ax5.set_ylabel(r'Latitude ($^{\circ}$)',fontsize=fs)
ax5.set_title(name + ' SW1 at ' + "{:.1f}".format(altitudes[0] / 1000) + ' km in U',fontsize=fs)
ax5.set_title(name + ' SW1 at 97 km in U',fontsize=fs)

# im6 = ax6.contourf(dt,lats,amps3[:,0,::1,wave3_ind].T,np.arange(smin,smax3,nstep3),vmin=0,vmax=plot_max3,cmap=cm)
# cs = ax6.contour(dt,lats,amps3[:,0,::1,wave3_ind].T, np.arange(smin,smax3,nstep3),linewidths=linefs,colors='black')
# ax6.clabel(cs, np.arange(smin,smax3,nstep3)[::contour_step], inline=1, fontsize=l_fs,fmt='%1.0f')
# ax6.set_ylim(zmin,zmax)
# ax6.set_xlim(tmin,tmax)
# ax6.tick_params(axis='x', which='major', labelsize=ts)
# ax6.tick_params(axis='y', which='major', labelsize=ts)
# ax6.set_ylabel(r'Latitude ($^{\circ}$)',fontsize=fs)
# ax6.set_title(name + ' SW3 at ' + "{:.1f}".format(altitudes[0] / 1000) + ' km in U',fontsize=fs)
# ax6.set_title(name + ' SW3 at 97 km in U',fontsize=fs)


cbar_ax = fig.add_axes([1, 1.3, 0.04, 0.5])
cb = fig.colorbar(im, cax=cbar_ax) 
cb.set_label(label='SW1 amplitude (ms$^{-1}$)', size=fs + 0)
cb.ax.tick_params(labelsize=ts + 2)

cbar_ax = fig.add_axes([1, 0.65, 0.04, 0.5])
cb = fig.colorbar(im2, cax=cbar_ax)
cb.set_label(label='SW2 amplitude (ms$^{-1}$)', size=fs + 0)
cb.ax.tick_params(labelsize=ts + 2)

cbar_ax = fig.add_axes([1, 0.0, 0.04, 0.5])
cb = fig.colorbar(im3, cax=cbar_ax)
cb.set_label(label=r'SW3 amplitude (ms$^{-1}$)', size=fs + 0)
cb.ax.tick_params(labelsize=ts + 2)

cbar_ax = fig.add_axes([2.35, 1.3, 0.04, 0.5])
cb = fig.colorbar(im4, cax=cbar_ax) 
cb.set_label(label='SW1 amplitude (ms$^{-1}$)', size=fs + 0)
cb.ax.tick_params(labelsize=ts + 2)

cbar_ax = fig.add_axes([2.35, 0.65, 0.04, 0.5])
cb = fig.colorbar(im5, cax=cbar_ax)
cb.set_label(label='SW2 amplitude (ms$^{-1}$)', size=fs + 0)
cb.ax.tick_params(labelsize=ts + 2)

cbar_ax = fig.add_axes([2.35, 0.0, 0.04, 0.5])
cb = fig.colorbar(im6, cax=cbar_ax)
cb.set_label(label=r'SW3 amplitude (ms$^{-1}$)', size=fs + 0)
cb.ax.tick_params(labelsize=ts + 2)

# ssw_c = 'black'
# ssw_l = 41
# ax1.vlines(ssw_l,0,100,color=ssw_c,ls='--',lw=lw_ssw) # x, ymin, ymax
# ax2.vlines(ssw_l,0,100,color=ssw_c,ls='--',lw=lw_ssw)
# ax3.vlines(ssw_l,0,100,color=ssw_c,ls='--',lw=lw_ssw)

# ax4.vlines(ssw_l,0,100,color=ssw_c,ls='--',lw=lw_ssw) # x, ymin, ymax
# ax5.vlines(ssw_l,0,100,color=ssw_c,ls='--',lw=lw_ssw)
# ax6.vlines(ssw_l,0,100,color=ssw_c,ls='--',lw=lw_ssw)

lw_ssw = 2.0

ssw_c = 'black'
ssw_l = 41
ssw_l2 = 34.5
ssw_l3 = 53
ax1.vlines(ssw_l,0,100,color=ssw_c,ls='--',lw=lw_ssw) # x, ymin, ymax
ax2.vlines(ssw_l,0,100,color=ssw_c,ls='--',lw=lw_ssw)
ax3.vlines(ssw_l,0,100,color=ssw_c,ls='--',lw=lw_ssw)
ax4.vlines(ssw_l,0,100,color=ssw_c,ls='--',lw=lw_ssw)
ax5.vlines(ssw_l,0,100,color=ssw_c,ls='--',lw=lw_ssw)
ax6.vlines(ssw_l,0,100,color=ssw_c,ls='--',lw=lw_ssw)

ax1.vlines(ssw_l2,0,100,color=ssw_c,ls='--',lw=lw_ssw) # x, ymin, ymax
ax2.vlines(ssw_l2,0,100,color=ssw_c,ls='--',lw=lw_ssw)
ax3.vlines(ssw_l2,0,100,color=ssw_c,ls='--',lw=lw_ssw)
ax4.vlines(ssw_l2,0,100,color=ssw_c,ls='--',lw=lw_ssw)
ax5.vlines(ssw_l2,0,100,color=ssw_c,ls='--',lw=lw_ssw)
ax6.vlines(ssw_l2,0,100,color=ssw_c,ls='--',lw=lw_ssw)

ax1.vlines(ssw_l3,0,100,color=ssw_c,ls='--',lw=lw_ssw) # x, ymin, ymax
ax2.vlines(ssw_l3,0,100,color=ssw_c,ls='--',lw=lw_ssw)
ax3.vlines(ssw_l3,0,100,color=ssw_c,ls='--',lw=lw_ssw)
ax4.vlines(ssw_l3,0,100,color=ssw_c,ls='--',lw=lw_ssw)
ax5.vlines(ssw_l3,0,100,color=ssw_c,ls='--',lw=lw_ssw)
ax6.vlines(ssw_l3,0,100,color=ssw_c,ls='--',lw=lw_ssw)

ssw_c = '0.675'
ax1.hlines(43.3,0,100,color=ssw_c,ls='--',lw=lw_ssw)
ax1.hlines(51.3,0,100,color=ssw_c,ls='--',lw=lw_ssw)
ax1.hlines(67.9,0,100,color=ssw_c,ls='--',lw=lw_ssw)

ax2.hlines(43.3,0,100,color=ssw_c,ls='--',lw=lw_ssw)
ax2.hlines(51.3,0,100,color=ssw_c,ls='--',lw=lw_ssw)
ax2.hlines(67.9,0,100,color=ssw_c,ls='--',lw=lw_ssw)

ax3.hlines(43.3,0,100,color=ssw_c,ls='--',lw=lw_ssw)
ax3.hlines(51.3,0,100,color=ssw_c,ls='--',lw=lw_ssw)
ax3.hlines(67.9,0,100,color=ssw_c,ls='--',lw=lw_ssw)

ax4.hlines(43.3,0,100,color=ssw_c,ls='--',lw=lw_ssw)
ax4.hlines(51.3,0,100,color=ssw_c,ls='--',lw=lw_ssw)
ax4.hlines(67.9,0,100,color=ssw_c,ls='--',lw=lw_ssw)

ax5.hlines(43.3,0,100,color=ssw_c,ls='--',lw=lw_ssw)
ax5.hlines(51.3,0,100,color=ssw_c,ls='--',lw=lw_ssw)
ax5.hlines(67.9,0,100,color=ssw_c,ls='--',lw=lw_ssw)

ax6.hlines(43.3,0,100,color=ssw_c,ls='--',lw=lw_ssw)
ax6.hlines(51.3,0,100,color=ssw_c,ls='--',lw=lw_ssw)
ax6.hlines(67.9,0,100,color=ssw_c,ls='--',lw=lw_ssw)

fw = 'bold'
ax1.text(17.5,85.,'(b)',fontsize=fs + 5,fontweight=fw)
ax2.text(17.5,85.,'(a)',fontsize=fs + 5,fontweight=fw)
ax3.text(17.5,85.,'(c)',fontsize=fs + 5,fontweight=fw)

# ax4.text(17.5,85.,'(e)',fontsize=fs + 5,fontweight=fw)
ax5.text(17.5,85.,'(d)',fontsize=fs + 5,fontweight=fw)
ax6.text(17.5,85.,'(f)',fontsize=fs + 5,fontweight=fw)

ax3.set_xlabel('Days since Dec 1st 2012',fontsize=fs)
ax4.set_xlabel('Days since Dec 1st 2012',fontsize=fs)

save_f = '/home/willemvc/Desktop/NTNU files/dope transfer MET/iSSW/paper/MET figures/'
# plt.savefig(save_f + 'fixedZM_HR.png', dpi=300, bbox_inches="tight")

#%%

plt.contourf(dt,lats,amps3[:,0,::1,wave1_ind].T - amps2[:,0,::1,wave1_ind].T)
plt.xlim(20,65)
plt.ylim(20,80)
plt.colorbar()

#%%

print(altitudes)