#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 15:56:26 2020

@author: wim
"""

import numpy as np
import h5py  
from matplotlib import pyplot as plt
import tide_func
from scipy.interpolate import griddata
from scipy.stats import circmean

#%


n_days = 16 # sliding window size
n_days_plt = n_days
model_dt = 1

l_in = 0 # 0 for S2 1 for M2 tide in 16-day window fit

#% read in DOPE output and sample at AVAILABLE SD locations for set years
file_n_1 = 'test58.dat'
file_n_2 = 'test25.dat'
var = 'u'
var2 = 'u'

radars = ['CMO'] # Bea, CMO, Col, Kir, Sod, Trd
# radars = ['Kir','Sod','Trd']
trd_off = np.array([0,20,40]) * 0
cords = [r'CMOR (43.3$^{\circ}$N, 80.8$^{\circ}$W)',
         r'Collm (51.3$^{\circ}$N, 13.0$^{\circ}$E)',
         r'Kiruna (67.9$^{\circ}$N, 21.1$^{\circ}$E)']

months = ['201212','201301','201302','201303']

data_len = 81
alt_step = 1 # sub-sample altitude range radars
alt_step_model = 1

for itms, rads in enumerate(radars):
    file_path = '/home/willemvc/Desktop/NTNU files/dope transfer MET/meteor_radars/' + rads + '/' 
    file_name = 'Meteor_radar_' + rads + '_GW_' + months[0] + '.h5'
    f1 = h5py.File(file_path + file_name,'r+')
    
    altitudes = f1['info']['alt'][0,:]
    radar_min = 80
    radar_max = 100
    
    min_ind = np.argmin(np.abs(altitudes - radar_min))
    max_ind = np.argmin(np.abs(altitudes - radar_max)) + 1
    
    altitudes = altitudes[min_ind:max_ind][::alt_step]
    
    group = 'wind'
    var_obs = var
    var_err = var + '_err'
    
    model_offset = 0 # altitude offset model (km) if desired
        
    for index, item in enumerate(months):
        file_path = '/home/willemvc/Desktop/NTNU files/dope transfer MET/meteor_radars/' + rads + '/'
        file_name = 'Meteor_radar_' + rads + '_GW_' + item + '.h5'
        f1 = h5py.File(file_path + file_name,'r+')
        
        if index == 0:
            data = f1[group][var_obs][:,min_ind:max_ind]
            error = f1[group][var_err][:,min_ind:max_ind]
            data = data[:,::alt_step]
            error = error[:,::alt_step]
            
        else:
            temp_d = f1[group][var_obs][:,min_ind:max_ind]
            temp_e = f1[group][var_err][:,min_ind:max_ind]
            temp_d = temp_d[:,::alt_step]
            temp_e = temp_e[:,::alt_step]
            
            data = np.concatenate((data,temp_d))
            error = np.concatenate((error,temp_e))
            
    data = data[:(data_len * 24),:] # custom time range
    
    if rads == 'CMO':
        sd_lats = np.array('43.3'.split(),
                           dtype='float')
        sd_lons = np.array([(360 - 80.8)])
        title = 'CMOR (43.3N,80.8W)'
        
    elif rads == 'Trd':
        sd_lats = np.array('64.4'.split(),
                       dtype='float') 
        sd_lons = np.array([10.5]) +  trd_off[itms]
        title = 'Trondheim (64.4N,10.5E)'
        
    elif rads == 'Sva':
        sd_lats = np.array('79.59'.split(), 
                       dtype='float')
        sd_lons = np.array([15.51])    
        title = 'Svalbard (79.59N,15.51E)'
        
    elif rads == 'Sod':
        sd_lats = np.array('67.25'.split(),
                       dtype='float')
        sd_lons = np.array([26.35])
        title = 'Sodankyla (67.25N,26.35E)'
        
    elif rads== 'Kir':
        sd_lats = np.array('67.88'.split(),
                        dtype='float')
        sd_lons = np.array([21.07])
        title = 'Kiruna (67.9N,21.1E)'

        # sd_lats = -np.array('67.57'.split(),
        #         dtype='float')
        # print(sd_lats)
        # sd_lons = np.array([360 - 68.13])
        # title = 'Rothera (67.57S, 68.13W)'
        
    elif rads == 'Col':
        sd_lats = np.array('51.3'.split(),
                       dtype='float')
        sd_lons = np.array([13.00])
        title = 'Collm (51.3N,13.0E)'
    
    elif rads == 'Eureka':
        sd_lats = np.array('79.59'.split(),
                       dtype='float')
        sd_lons = np.array([360 - 56.27])
        title = 'Eureka (79.59N,56.27W)'
        
    elif rads == 'Bea':
        sd_lats = np.array('74.26'.split(),
                       dtype='float')
        sd_lons = np.array([19.02])
        title = 'Bear Island (74.26N, 19.02E)'
        
    f1.close()
    
    # from here on, 'data' holds meteor radar data in time-alt
    
    #%
     
    save_f = '/home/willemvc/Desktop/NTNU files/dope transfer MET/iSSW/data/sigprim_fits/'
    
    full_v_1_temp = np.load(save_f + file_n_1 + '_' + var + '.npy')[:,:,:,:]
    full_v_2_temp = np.load(save_f + file_n_2 + '_' + var + '.npy')[:,:,:,:]
    
    height_arr = np.load(save_f + file_n_1 + '_height_array.npy')[:]
    dope_lats = np.load(save_f + file_n_1 + '_lat.npy')[:]
    dope_lons = np.load(save_f + file_n_1 + '_lon.npy')[:]
    
    min_ind = np.argmin(np.abs(height_arr - radar_min * 1000))
    max_ind = np.argmin(np.abs(height_arr - radar_max * 1000)) + 1
    
    height_arr = height_arr[max_ind:min_ind]
    height_arr = height_arr[::alt_step_model]
    
    # read in only northern hemisphere above 30 degrees North to save space
    lat_min = 40
    lat_max = 80
    
    lats = dope_lats
    lons = dope_lons
    
    lat_ind_min = np.argmin(np.abs(lats - lat_min))
    lat_ind_max = np.argmin(np.abs(lats - lat_max))
    
    lats = lats[lat_ind_min:lat_ind_max] 
    
    full_v_1_temp = full_v_1_temp[:,max_ind:min_ind,lat_ind_min:lat_ind_max,:]
    full_v_1_temp = full_v_1_temp[:,::alt_step_model,:,:]
    full_v_1 = np.copy(full_v_1_temp) # dt, dz, dy, dx
    
    full_v_2_temp = full_v_2_temp[:,max_ind:min_ind,lat_ind_min:lat_ind_max,:]
    full_v_2_temp = full_v_2_temp[:,::alt_step_model,:,:]
    full_v_2 = np.copy(full_v_2_temp) # dt, dz, dy, dx
    
    del full_v_1_temp
    del full_v_2_temp
    
    #%
    
    dt, dz, dy, dx = np.shape(full_v_1)
    sd_interp = np.zeros((dt, dz))
    sd_interp2 = np.zeros((dt, dz))
    
    lon_x, lat_x = np.meshgrid(lons,lats)   
    
    # lat, lon indexing
    points = np.array((lat_x.flatten(),lon_x.flatten())).T
    sd_stations = list(zip(sd_lats,sd_lons)) # single radar locaction
    
    # interpolate to station site
    for z in range(dz):
        for t in range(dt):       
            temp = np.copy(full_v_1[t,z,:,:])
            temp2 = np.copy(full_v_2[t,z,:,:])
            values = temp.flatten()
            values2 = temp2.flatten()
            for index, item in enumerate(sd_stations):
                sd_interp[t,z] = griddata(points,
                                          values,
                                          item,
                                          method='cubic')
                sd_interp2[t,z] = griddata(points,
                                           values2,
                                           item,
                                           method='cubic')
                
        print('interpolating z-level ', z)
                
    #%
    
    dt, dz = np.shape(sd_interp)
    sd_interp_ext = np.zeros((dt * 3,dz))
    sd_interp_ext2 = np.zeros((dt * 3,dz))
    
    # interpolate from model time-grid to hourly
    if model_dt == 3:
        for z in range(dz):
            t_array = np.arange(np.shape(sd_interp[:,z])[0] * 3) / 3
            t_in = np.arange(np.shape(sd_interp[:,z])[0])
            sd_interp_ext[:,z] = np.interp(t_array,t_in,sd_interp[:,z])
            
    else: 
        sd_interp_ext = np.copy(sd_interp)
        sd_interp_ext2 = np.copy(sd_interp2)
        
    #%
    
    half_w = int(n_days / 2 * 24) 
    half_w_plt = int(n_days_plt / 2 * 24)
    dat_len = np.shape(data)[0]
    
    dt, dz = np.shape(data)
    vars_sd = np.zeros((dat_len - n_days_plt * 24,dz,9)) # DT, SDT, TDT, mean
    
    for z in range(dz): 
        vars_sd[:,z,:] = tide_func.trondheim_m2s2_fit(data[:,z],
                                                      n_days * 24,
                                                      2 * 24)[half_w_plt:(dat_len - half_w_plt),:]
        
        # trondheim_tide_fit(data[:,z],
          #                                            n_days * 24 ,
           #                                           2 * 24)[half_w:(dat_len - half_w),:]
        
        print('data tide z-level ', z)
    
    dt, dz = np.shape(sd_interp_ext)
    vars_na = np.zeros((dat_len - n_days_plt * 24,dz,9)) # DT, SDT, TDT, mean
    vars_na2 = np.zeros((dat_len - n_days_plt * 24,dz,9)) # DT, SDT, TDT, mean
    
    #%
    
    for z in range(dz):
        vars_na[:,z,:] = tide_func.trondheim_m2s2_fit(sd_interp_ext[:,z],
                                                      n_days * 24,
                                                      2 * 24)[half_w_plt:(dat_len - half_w_plt),:]
        vars_na2[:,z,:] = tide_func.trondheim_m2s2_fit(sd_interp_ext2[:,z],
                                                      n_days * 24,
                                                      2 * 24)[half_w_plt:(dat_len - half_w_plt),:]
        
        print('model tide z-level ', z)
    
    #%
    
    md = 1                      # number of days in running mean
    sample_r = int(12)          # phases sampled every sample_r hours
    
    center = np.pi
    
    # correct phase by 180 degrees when A < 0 and time offset to 00:00 local
    inds_sdt_sd = np.where(vars_sd[:,:,1] < 0)
    vars_sd_corr = np.copy(vars_sd) + (4 * np.pi / 360) * sd_lons[0]
    vars_sd_corr[:,:,4][inds_sdt_sd] += np.pi 
    
    inds_sdt_na = np.where(vars_na[:,:,1] < 0)
    vars_na_corr = np.copy(vars_na) + (4 * np.pi / 360) * sd_lons[0]
    vars_na_corr[:,:,4][inds_sdt_na] += np.pi
    
    inds_sdt_na2 = np.where(vars_na2[:,:,1] < 0)
    vars_na_corr2 = np.copy(vars_na2) + (4 * np.pi / 360) * sd_lons[0]
    vars_na_corr2[:,:,4][inds_sdt_na2] += np.pi
    
    # sdt_sd = tide_func.running_mean(np.abs(vars_sd[::sample_r,1]),md)
    # sdt_na = tide_func.running_mean(np.abs(vars_na[::sample_r,1]),md)
    
    if itms == 0:
        tot_sd_s2 = np.zeros((np.shape(vars_sd[:,:,1])[0],np.shape(vars_sd[:,:,1])[1],3))
        tot_sd_m2 = np.zeros((np.shape(vars_sd[:,:,1])[0],np.shape(vars_sd[:,:,1])[1],3))
        tot_na_s2 = np.zeros((np.shape(vars_na[:,:,1])[0],np.shape(vars_na[:,:,1])[1],3))
        tot_na_m2 = np.zeros((np.shape(vars_na[:,:,1])[0],np.shape(vars_na[:,:,1])[1],3))
        tot_na_s22 = np.zeros((np.shape(vars_na2[:,:,1])[0],np.shape(vars_na2[:,:,1])[1],3))
        tot_na_m22 = np.zeros((np.shape(vars_na2[:,:,1])[0],np.shape(vars_na2[:,:,1])[1],3))
        
        tot_sd_ph = np.zeros((np.shape(vars_sd_corr[:,:,1])[0],np.shape(vars_sd_corr[:,:,1])[1],3))
        tot_na_ph = np.zeros((np.shape(vars_na_corr[:,:,1])[0],np.shape(vars_na_corr[:,:,1])[1],3))
        tot_na_ph2 = np.zeros((np.shape(vars_na_corr2[:,:,1])[0],np.shape(vars_na_corr2[:,:,1])[1],3))
        
        tot_sd_s2[:,:,itms] = np.abs(vars_sd[:,:,0])
        tot_sd_m2[:,:,itms] = np.abs(vars_sd[:,:,1])
        tot_na_s2[:,:,itms] = np.abs(vars_na[:,:,0])
        tot_na_m2[:,:,itms] = np.abs(vars_na[:,:,1])
        tot_na_s22[:,:,itms] = np.abs(vars_na2[:,:,0])
        tot_na_m22[:,:,itms] = np.abs(vars_na2[:,:,1])
        
        tot_sd_ph[:,:,itms] = vars_sd_corr[:,:,4]
        tot_na_ph[:,:,itms] = vars_na_corr[:,:,4]
        tot_na_ph2[:,:,itms] = vars_na_corr2[:,:,4]
        
    else:
        tot_sd_s2[:,:,itms] = np.abs(vars_sd[:,:,0])
        tot_sd_m2[:,:,itms] = np.abs(vars_sd[:,:,1])
        tot_na_s2[:,:,itms] = np.abs(vars_na[:,:,0])
        tot_na_m2[:,:,itms] = np.abs(vars_na[:,:,1])
        tot_na_s22[:,:,itms] = np.abs(vars_na2[:,:,0])
        tot_na_m22[:,:,itms] = np.abs(vars_na2[:,:,1])
        
        tot_sd_ph[:,:,itms] = vars_sd_corr[:,:,4]
        tot_na_ph[:,:,itms] = vars_na_corr[:,:,4]
        tot_na_ph2[:,:,itms] = vars_na_corr2[:,:,4]
        
    print(itms, rads)
    
# #%%
# plt.pcolormesh(tot_sd_m2[:,:,0])
# #%%
#%

# UTC time of maximum
dx, dy, dz = np.shape(tot_sd_ph)
for x in range(dx):
    for y in range(dy):
        for z in range(dz):
            tot_sd_ph[x,y,z] = 12 / (2 * np.pi) * circmean(np.pi / 2 - tot_sd_ph[x,y,z])

dx, dy, dz = np.shape(tot_na_ph)
for x in range(dx):
    for y in range(dy):
        for z in range(dz):
            tot_na_ph[x,y,z] = 12 / (2 * np.pi) * circmean(np.pi / 2 - tot_na_ph[x,y,z])
            
#%%

smin = 0
smax = 51
nstep = 5

smax2 = 34
nstep2 = 4

linefs = 0.3
l_fs = 15.
contour_step = 1

tmin = 20
tmax = 65

zmin = 85
zmax = 97

ts = 17.3
fs = 21.
lw_ssw = 3.0
ssw_ls = '--'
ssw_c = 'white'

cm = 'viridis'
plot_max = 47.5
plot_max2 = 31.0

dy, dt = np.shape(tot_sd_s2[:,:,0].T)
sd_time = np.arange(dt) / 24 + int(n_days_plt / 2)

dy, dt = np.shape(tot_na_s2[:,:,0].T)
na_time = np.arange(dt) / 24 + int(n_days_plt / 2)

fig = plt.figure(figsize=(6.5,6.5))

ax1 = fig.add_axes([0.0, 1.30, 1.0, 0.5]) # [xmin, ymin, dx, dy]
ax2 = fig.add_axes([0.0, 0.65, 1.0, 0.5]) 
ax3 = fig.add_axes([0.0, 0.00, 1.0, 0.5])

ax4 = fig.add_axes([1.4, 1.30, 1.0, 0.5])
ax5 = fig.add_axes([1.4, 0.65, 1.0, 0.5])
ax6 = fig.add_axes([1.4, 0.00, 1.0, 0.5])

# ax4 = fig.add_axes([1.4, 1.30, 1.0, 0.5])
# ax5 = fig.add_axes([1.4, 0.65, 1.0, 0.5])

im = ax1.contourf(sd_time,altitudes,tot_sd_s2[:,:,0].T,np.arange(smin,smax,nstep),vmin=0,vmax=plot_max,cmap=cm)
cs = ax1.contour(sd_time,altitudes,tot_sd_s2[:,:,0].T, np.arange(smin,smax,nstep),linewidths=linefs,colors='black')
ax1.clabel(cs, np.arange(smin,smax,nstep)[::contour_step], inline=1, fontsize=l_fs,fmt='%1.0f')
ax1.set_ylim(zmin,zmax)
ax1.set_xlim(tmin,tmax)
ax1.tick_params(axis='x', which='major', labelsize=ts)
ax1.tick_params(axis='y', which='major', labelsize=ts)
ax1.set_ylabel('Altitude (km)',fontsize=fs)
ax1.set_title('Solar SDT ' + cords[0],fontsize=fs)

ax4.contourf(sd_time,altitudes,tot_sd_m2[:,:,0].T,np.arange(smin,smax2,nstep2),vmin=0,vmax=plot_max2,cmap=cm)
cs = ax4.contour(sd_time,altitudes,tot_sd_m2[:,:,0].T,np.arange(smin,smax2,nstep2),linewidths=linefs,colors='black')
ax4.clabel(cs, np.arange(smin,smax2,nstep2)[::contour_step], inline=1, fontsize=l_fs,fmt='%1.0f')
ax4.set_ylim(zmin,zmax)
ax4.set_xlim(tmin,tmax)
ax4.tick_params(axis='x', which='major', labelsize=ts)
ax4.tick_params(axis='y', which='major', labelsize=ts)
ax4.set_ylabel('Altitude (km)',fontsize=fs)
ax4.set_title('Lunar SDT ' + cords[0],fontsize=fs)

#%

# ax3.contourf(sd_time,altitudes,tot_sd[:,:,2].T,np.arange(smin,smax,nstep),vmin=0,vmax=plot_max,cmap=cm)
# cs = ax3.contour(sd_time,altitudes,tot_sd[:,:,2].T,np.arange(smin,smax,nstep),linewidths=linefs,colors='black')
# ax3.clabel(cs, np.arange(smin,smax,nstep)[::contour_step], inline=1, fontsize=l_fs,fmt='%1.0f')
# ax3.set_ylim(zmin,zmax)
# ax3.set_xlim(tmin,tmax)
# ax3.tick_params(axis='x', which='major', labelsize=ts)
# ax3.tick_params(axis='y', which='major', labelsize=ts)
# ax3.set_ylabel('Altitude (km)',fontsize=fs)
# ax3.set_title(cords[2],fontsize=fs)

ax2.contourf(na_time,height_arr / 1000,tot_na_s2[:,:,0].T,np.arange(smin,smax,nstep),vmin=0,vmax=plot_max,cmap=cm)
cs = ax2.contour(na_time,height_arr / 1000,tot_na_s2[:,:,0].T,np.arange(smin,smax,nstep),linewidths=linefs,colors='black')
ax2.clabel(cs, np.arange(smin,smax,nstep)[::contour_step], inline=1, fontsize=l_fs,fmt='%1.0f')
ax2.set_ylim(zmin,zmax)
ax2.set_xlim(tmin,tmax)
ax2.tick_params(axis='x', which='major', labelsize=ts)
ax2.tick_params(axis='y', which='major', labelsize=ts)
ax2.set_ylabel('Altitude (km)',fontsize=fs)
ax2.set_title('PRISM Solar SDT ' + cords[0],fontsize=fs)

im2 = ax5.contourf(na_time,height_arr / 1000,tot_na_m2[:,:,0].T,np.arange(smin,smax2,nstep2),vmin=0,vmax=plot_max2,cmap=cm)
cs = ax5.contour(na_time,height_arr / 1000,tot_na_m2[:,:,0].T,np.arange(smin,smax2,nstep2),linewidths=linefs,colors='black')
ax5.clabel(cs, np.arange(smin,smax2,nstep2)[::contour_step], inline=1, fontsize=l_fs,fmt='%1.0f')
ax5.set_ylim(zmin,zmax)
ax5.set_xlim(tmin,tmax)
ax5.tick_params(axis='x', which='major', labelsize=ts)
ax5.tick_params(axis='y', which='major', labelsize=ts)
ax5.set_ylabel('Altitude (km)',fontsize=fs)
ax5.set_title('PRISM Lunar SDT ' + cords[0],fontsize=fs)

ax3.contourf(na_time,height_arr / 1000,tot_na_s22[:,:,0].T,np.arange(smin,smax,nstep),vmin=0,vmax=plot_max,cmap=cm)
cs = ax3.contour(na_time,height_arr / 1000,tot_na_s22[:,:,0].T,np.arange(smin,smax,nstep),linewidths=linefs,colors='black')
ax3.clabel(cs, np.arange(smin,smax,nstep)[::contour_step], inline=1, fontsize=l_fs,fmt='%1.0f')
ax3.set_ylim(zmin,zmax)
ax3.set_xlim(tmin,tmax)
ax3.tick_params(axis='x', which='major', labelsize=ts)
ax3.tick_params(axis='y', which='major', labelsize=ts)
ax3.set_ylabel('Altitude (km)',fontsize=fs)
ax3.set_title('OnlySolar Solar SDT ' + cords[0],fontsize=fs)

ax6.contourf(na_time,height_arr / 1000,tot_na_m22[:,:,0].T,np.arange(smin,smax2,nstep2),vmin=0,vmax=plot_max2,cmap=cm)
cs = ax6.contour(na_time,height_arr / 1000,tot_na_m22[:,:,0].T,np.arange(smin,smax2,nstep2),linewidths=linefs,colors='black')
ax6.clabel(cs, np.arange(smin,smax2,nstep2)[::contour_step], inline=1, fontsize=l_fs,fmt='%1.0f')
ax6.set_ylim(zmin,zmax)
ax6.set_xlim(tmin,tmax)
ax6.tick_params(axis='x', which='major', labelsize=ts)
ax6.tick_params(axis='y', which='major', labelsize=ts)
ax6.set_ylabel('Altitude (km)',fontsize=fs)
ax6.set_title('OnlySolar Lunar SDT ' + cords[0],fontsize=fs)

#%
# ax6.contourf(na_time,height_arr / 1000,tot_na[:,:,2].T,np.arange(smin,smax,nstep),vmin=0,vmax=plot_max,cmap=cm)
# cs = ax6.contour(na_time,height_arr / 1000,tot_na[:,:,2].T,np.arange(smin,smax,nstep),linewidths=linefs,colors='black')
# ax6.clabel(cs, np.arange(smin,smax,nstep)[::contour_step], inline=1, fontsize=l_fs,fmt='%1.0f')
# ax6.set_ylim(zmin,zmax)
# ax6.set_xlim(tmin,tmax)
# ax6.tick_params(axis='x', which='major', labelsize=ts)
# ax6.tick_params(axis='y', which='major', labelsize=ts)
# ax6.set_ylabel('Altitude (km)',fontsize=fs)
# ax6.set_title('Model ' + cords[2],fontsize=fs)

fw = 'bold'
ax1.text(16.5,97.5,'(a)',fontsize=fs + 7,fontweight=fw)
ax2.text(16.5,97.5,'(b)',fontsize=fs + 7,fontweight=fw)
ax3.text(16.0,97.5,'(c)',fontsize=fs + 7,fontweight=fw)
ax4.text(16.5,97.5,'(d)',fontsize=fs + 7,fontweight=fw)
ax5.text(16.5,97.5,'(e)',fontsize=fs + 7,fontweight=fw)
ax6.text(16.0,97.5,'(f)',fontsize=fs + 7,fontweight=fw)

cbar_ax = fig.add_axes([2.45, 0.325, 0.04, 1.15])
cb = fig.colorbar(im2, cax=cbar_ax) 
cb.set_label(label=r'M2 amplitude (ms$^{-1}$) in U', size=fs + 2)
cb.ax.tick_params(labelsize=ts + 2)

cbar_ax = fig.add_axes([1.05, 0.325, 0.04, 1.15])
cb = fig.colorbar(im, cax=cbar_ax) 
cb.set_label(label=r'S2 amplitude (ms$^{-1}$) in U', size=fs + 2)
cb.ax.tick_params(labelsize=ts + 2)

ssw_t = 41
ssw_t2 = 34.5
ssw_t3 = 54
ax1.vlines(ssw_t,0,100,color=ssw_c,ls=ssw_ls,lw=lw_ssw) # x, ymin, ymax
ax2.vlines(ssw_t,0,100,color=ssw_c,ls=ssw_ls,lw=lw_ssw)
ax4.vlines(ssw_t,0,100,color=ssw_c,ls=ssw_ls,lw=lw_ssw)
ax5.vlines(ssw_t,0,100,color=ssw_c,ls=ssw_ls,lw=lw_ssw)
ax3.vlines(ssw_t,0,100,color=ssw_c,ls=ssw_ls,lw=lw_ssw)
ax6.vlines(ssw_t,0,100,color=ssw_c,ls=ssw_ls,lw=lw_ssw)

ax1.vlines(ssw_t2,0,100,color=ssw_c,ls=ssw_ls,lw=lw_ssw) # x, ymin, ymax
ax2.vlines(ssw_t2,0,100,color=ssw_c,ls=ssw_ls,lw=lw_ssw)
ax4.vlines(ssw_t2,0,100,color=ssw_c,ls=ssw_ls,lw=lw_ssw)
ax5.vlines(ssw_t2,0,100,color=ssw_c,ls=ssw_ls,lw=lw_ssw)
ax3.vlines(ssw_t2,0,100,color=ssw_c,ls=ssw_ls,lw=lw_ssw)
ax6.vlines(ssw_t2,0,100,color=ssw_c,ls=ssw_ls,lw=lw_ssw)

ax1.vlines(ssw_t3,0,100,color=ssw_c,ls=ssw_ls,lw=lw_ssw) # x, ymin, ymax
ax2.vlines(ssw_t3,0,100,color=ssw_c,ls=ssw_ls,lw=lw_ssw)
ax4.vlines(ssw_t3,0,100,color=ssw_c,ls=ssw_ls,lw=lw_ssw)
ax5.vlines(ssw_t3,0,100,color=ssw_c,ls=ssw_ls,lw=lw_ssw)
ax3.vlines(ssw_t3,0,100,color=ssw_c,ls=ssw_ls,lw=lw_ssw)
ax6.vlines(ssw_t3,0,100,color=ssw_c,ls=ssw_ls,lw=lw_ssw)

ax3.set_xlabel('Days since Dec 1st 2012',fontsize=fs)
ax6.set_xlabel('Days since Dec 1st 2012',fontsize=fs)

print(file_n_1,var)

save_f = '/home/wim/Desktop/Projects/iSSW/paper/figures/'
# plt.savefig(save_f + 's2m2_CMOR_solar.png',dpi=300, bbox_inches="tight")

#%%

smin = 0
smax = 13
nstep = 1

linefs = 0.3
l_fs = 14.
contour_step = 1

tmin = 20
tmax = 65

zmin = 85
zmax = 97

ts = 17.2
fs = 20.
lw_ssw = 3.0

cm = 'jet'
plot_max = 12

dy, dt = np.shape(tot_sd[:,:,0].T)
sd_time = np.arange(dt) / 24
sd_time = sd_time[::24]

dy, dt = np.shape(tot_na[:,:,0].T)
na_time = np.arange(dt) / 24
na_time = na_time[::24]

fig = plt.figure(figsize=(6.5,6.5))

ax1 = fig.add_axes([0.0, 1.30, 1.0, 0.5]) # [xmin, ymin, dx, dy]
ax2 = fig.add_axes([0.0, 0.65, 1.0, 0.5]) 
ax3 = fig.add_axes([0.0, 0.00, 1.0, 0.5])

ax4 = fig.add_axes([1.2, 1.30, 1.0, 0.5])
ax5 = fig.add_axes([1.2, 0.65, 1.0, 0.5])
ax6 = fig.add_axes([1.2, 0.00, 1.0, 0.5])

ax1.contourf(sd_time,altitudes,tot_sd_ph[::24,:,0].T,np.arange(smin,smax,nstep),vmin=0,vmax=plot_max,cmap=cm)
cs = ax1.contour(sd_time,altitudes,tot_sd_ph[::24,:,0].T, np.arange(smin,smax,nstep),linewidths=linefs,colors='black')
ax1.clabel(cs, np.arange(smin,smax,nstep)[::contour_step], inline=1, fontsize=l_fs,fmt='%1.0f')
ax1.set_ylim(zmin,zmax)
ax1.set_xlim(tmin,tmax)
ax1.tick_params(axis='x', which='major', labelsize=ts)
ax1.tick_params(axis='y', which='major', labelsize=ts)
ax1.set_ylabel('Altitude (km)',fontsize=fs)
ax1.set_title(cords[0],fontsize=fs)

ax2.contourf(sd_time,altitudes,tot_sd_ph[::24,:,1].T,np.arange(smin,smax,nstep),vmin=0,vmax=plot_max,cmap=cm)
cs = ax2.contour(sd_time,altitudes,tot_sd_ph[::24,:,1].T,np.arange(smin,smax,nstep),linewidths=linefs,colors='black')
ax2.clabel(cs, np.arange(smin,smax,nstep)[::contour_step], inline=1, fontsize=l_fs,fmt='%1.0f')
ax2.set_ylim(zmin,zmax)
ax2.set_xlim(tmin,tmax)
ax2.tick_params(axis='x', which='major', labelsize=ts)
ax2.tick_params(axis='y', which='major', labelsize=ts)
ax2.set_ylabel('Altitude (km)',fontsize=fs)
ax2.set_title(cords[1],fontsize=fs)

ax3.contourf(sd_time,altitudes,tot_sd_ph[::24,:,2].T,np.arange(smin,smax,nstep),vmin=0,vmax=plot_max,cmap=cm)
cs = ax3.contour(sd_time,altitudes,tot_sd_ph[::24,:,2].T,np.arange(smin,smax,nstep),linewidths=linefs,colors='black')
ax3.clabel(cs, np.arange(smin,smax,nstep)[::contour_step], inline=1, fontsize=l_fs,fmt='%1.0f')
ax3.set_ylim(zmin,zmax)
ax3.set_xlim(tmin,tmax)
ax3.tick_params(axis='x', which='major', labelsize=ts)
ax3.tick_params(axis='y', which='major', labelsize=ts)
ax3.set_ylabel('Altitude (km)',fontsize=fs)
ax3.set_title(cords[2],fontsize=fs)

ax4.contourf(na_time,height_arr / 1000,tot_na_ph[::24,:,0].T,np.arange(smin,smax,nstep),vmin=0,vmax=plot_max,cmap=cm)
cs = ax4.contour(na_time,height_arr / 1000,tot_na_ph[::24,:,0].T,np.arange(smin,smax,nstep),linewidths=linefs,colors='black')
ax4.clabel(cs, np.arange(smin,smax,nstep)[::contour_step], inline=1, fontsize=l_fs,fmt='%1.0f')
ax4.set_ylim(zmin,zmax)
ax4.set_xlim(tmin,tmax)
ax4.tick_params(axis='x', which='major', labelsize=ts)
ax4.tick_params(axis='y', which='major', labelsize=ts)
ax4.set_ylabel('Altitude (km)',fontsize=fs)
ax4.set_title('Model ' + cords[0],fontsize=fs)

im = ax5.contourf(na_time,height_arr / 1000,tot_na_ph[::24,:,1].T,np.arange(smin,smax,nstep),vmin=0,vmax=plot_max,cmap=cm)
cs = ax5.contour(na_time,height_arr / 1000,tot_na_ph[::24,:,1].T,np.arange(smin,smax,nstep),linewidths=linefs,colors='black')
ax5.clabel(cs, np.arange(smin,smax,nstep)[::contour_step], inline=1, fontsize=l_fs,fmt='%1.0f')
ax5.set_ylim(zmin,zmax)
ax5.set_xlim(tmin,tmax)
ax5.tick_params(axis='x', which='major', labelsize=ts)
ax5.tick_params(axis='y', which='major', labelsize=ts)
ax5.set_ylabel('Altitude (km)',fontsize=fs)
ax5.set_title('Model ' + cords[1],fontsize=fs)

ax6.contourf(na_time,height_arr / 1000,tot_na_ph[::24,:,2].T,np.arange(smin,smax,nstep),vmin=0,vmax=plot_max,cmap=cm)
cs = ax6.contour(na_time,height_arr / 1000,tot_na_ph[::24,:,2].T,np.arange(smin,smax,nstep),linewidths=linefs,colors='black')
ax6.clabel(cs, np.arange(smin,smax,nstep)[::contour_step], inline=1, fontsize=l_fs,fmt='%1.0f')
ax6.set_ylim(zmin,zmax)
ax6.set_xlim(tmin,tmax)
ax6.tick_params(axis='x', which='major', labelsize=ts)
ax6.tick_params(axis='y', which='major', labelsize=ts)
ax6.set_ylabel('Altitude (km)',fontsize=fs)
ax6.set_title('Model ' + cords[2],fontsize=fs)

fig.subplots_adjust(right=0.8)
# cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])

cbar_ax = fig.add_axes([2.25, 0.35, 0.04, 1.10])
cb = fig.colorbar(im, cax=cbar_ax) # ,pad=0.01,fraction=0.1
cb.set_label(label='SDT local time of maximum (hr) in U', size=fs + 2)
cb.ax.tick_params(labelsize=ts + 2)

# ax1.vlines(37,0,100,color='red',ls='--',lw=lw_ssw) # x, ymin, ymax
# ax2.vlines(37,0,100,color='red',ls='--',lw=lw_ssw)
# ax3.vlines(37,0,100,color='red',ls='--',lw=lw_ssw)
# ax4.vlines(37,0,100,color='red',ls='--',lw=lw_ssw)
# ax5.vlines(37,0,100,color='red',ls='--',lw=lw_ssw)
# ax6.vlines(37,0,100,color='red',ls='--',lw=lw_ssw)

ax3.set_xlabel('Days since Dec 1st 2012',fontsize=fs)
ax6.set_xlabel('Days since Dec 1st 2012',fontsize=fs)

print(file_n_1,var)

save_f = '/home/wim/Desktop/Projects/iSSW/paper/figures/'
# plt.savefig(save_f + 'U_pha_SDT.png',dpi=300, bbox_inches="tight")