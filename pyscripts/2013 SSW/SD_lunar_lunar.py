#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 16 11:24:54 2019

@author: wim
"""

import numpy as np
from matplotlib import pyplot as plt
import tide_func
from scipy.stats import circmean
import matplotlib.patches as mpatches
import seaborn as sns
from datetime import datetime, timedelta


"""
Input wave surfaces / climatology is generated at bottom of this file.

"""

#

year = '2013'
file_n_1 = 'test58.dat' # SIGPRIM-SD
file_n_2 = 'test25.dat' # onlythermal
file_n_3 = 'test26.dat' # onlylunar

file_n_1 = 'test58.dat' # SIGPRIM-SD
file_n_2 = 'test2.dat' # onlythermal
file_n_3 = 'test3.dat' # onlylunar

invar = 8 # 7, 8 = V,U

mean_height = 95000
FWHM_gauss = 25000

model_dt = 1

inc_sw1 = False
inc_sw1_sw3 = False
inc_sw1_se2 = False

inc_ind = 6 # 6 for first extra fit, 7 for 2nd

extended = False
ext = ''

n_days_fit = 81
offset = timedelta(days=-(121 - n_days_fit))
all_ind = int((24 / model_dt) * n_days_fit)

if inc_sw1:
    ext = 'sw1'
    extended = True
elif inc_sw1_sw3:
    ext = 'sw1_sw3'
    extended = True
elif inc_sw1_se2:
    ext = 'sw1_se2'
    extended = True

if invar == 7:
    var = 'v'
    var2 = 'v'
if invar == 8:
    var = 'u'
    var2 = 'u'

n_days = 4
window = int(n_days * 24) # superdarn
min_datapoints = 1 # 960 / 2
ww = 'ww_' + str(int(window / 24)) + 'd_'

if year == '2009':
    start_day = datetime(2008,12,1,0) 
    end_day = datetime(2009,4,1,0) 
    
if year == '2009_dart':
    start_day = datetime(2009,1,10) 
    end_day = datetime(2009,3,20) 
    
elif year == '2010':
    start_day = datetime(2009,12,1,0) 
    end_day = datetime(2010,4,1,0) 
    
elif year == '2013':
    start_day = datetime(2012,12,1,0) 
    end_day = datetime(2013,4,1,0) + offset
    # end_day = start_day + timedelta(days=45)

elif year == '2015':
    start_day = datetime(2014,12,1,0) 
    end_day = datetime(2015,4,1,0) 
    
elif year == '2016':
    start_day = datetime(2015,12,1,0) 
    end_day = datetime(2016,4,1,0) 
    
elif year == '2017':
    start_day = datetime(2016,12,1,0)
    end_day = datetime(2017,4,1,0) 
    
data, std, longitudes = tide_func.new_superdarn(start_day,
                                                end_day,
                                                invar,
                                                optimized_spread=True,
                                                n_skim='100')

if invar == 8:
    data = -data

# data, longitudes = tide_func.updated_superdarn_mc_skim(start_day,
#                                                         end_day,
#                                                         var=var,
#                                                         optimized_spread=True,
#                                                         skim_n='50')    

varibs2, dt, sdt, tdt, fit_er_v = tide_func.superdarn_dt_sdt_tdt_mean_fitter_std(data,
                                                                                 std,
                                                                                 longitudes,
                                                                                 32,
                                                                                 window,
                                                                                 24,
                                                                                 min_datapoints,
                                                                                 10)

soup_sample = int(24)  # match SIGPRIM output time resolution
smw = 4                # four day smoothing window

fit_v_amp = fit_er_v[::soup_sample,1]
fit_v_phas = fit_er_v[::soup_sample,4]    

if inc_sw1:
    varibs2, dt, sdt, tdt, fit_er_v = tide_func.superdarn_dt_sdt_tdt_mean_sw1_fitter(data,
                                                                                     longitudes,
                                                                                     32,
                                                                                     window,
                                                                                     soup_sample,
                                                                                     min_datapoints)
    
elif inc_sw1_sw3:
    varibs2, dt, sdt, tdt, fit_er_v = tide_func.superdarn_dt_sdt_tdt_mean_sw1_sw3_fitter(data,
                                                                                          longitudes,
                                                                                          32,
                                                                                          window,
                                                                                          soup_sample,
                                                                                          min_datapoints)
    
elif inc_sw1_se2:
    varibs2, dt, sdt, tdt, fit_er_v = tide_func.superdarn_dt_sdt_tdt_mean_sw1_se1_fitter(data,
                                                                                          longitudes,
                                                                                          32,
                                                                                          window,
                                                                                          soup_sample,
                                                                                          min_datapoints)
    print('hello')
    
else:
    varibs2, dt, sdt, tdt, fit_er_v = tide_func.superdarn_dt_sdt_tdt_mean_fitter(data,
                                                                                 longitudes,
                                                                                 32,
                                                                                 window,
                                                                                 soup_sample,
                                                                                 min_datapoints)
#%

# 5 additional days due to leap years before 2014
dt_14_15 = np.abs(varibs2[::soup_sample,0])
sdt_14_15 = np.abs(varibs2[::soup_sample,1])
tdt_14_15 = np.abs(varibs2[::soup_sample,2])

dt_14_15 = tide_func.running_mean(dt_14_15,smw)
sdt_14_15 = tide_func.running_mean(sdt_14_15,smw)
tdt_14_15 = tide_func.running_mean(tdt_14_15,smw)

# sub-sample starting from 00:00 hrs
phase_dt = np.copy(varibs2[::soup_sample,3])
phase_sdt = np.copy(varibs2[::soup_sample,4])
phase_tdt = np.copy(varibs2[::soup_sample,5])

# correct negative amplitude phase offset
inds_sdt = np.where(varibs2[::soup_sample,1] < 0)
phase_sdt[inds_sdt] += np.pi

center = np.pi # center phase-curve in middle of figure with [0, 2 pi] y-axis

mean_sdt_v = circmean(phase_sdt)
mean_lotm_sdt_v = circmean(np.pi / 2 - mean_sdt_v)

v_midphase = 6 # hours at center of phase-plot (LTOM)
sdt_offset_v = -(v_midphase - mean_lotm_sdt_v / (2 * np.pi) * 12) / 12 * (2 * np.pi)

for k in range(np.shape(phase_sdt)[0]):
        phase_sdt[k] = circmean(np.pi / 2 - phase_sdt[k])
        
#% dope output 1

output_f = '/home/wim/Desktop/Projects/iSSW/data/SD_sampled/'
n_days = 4
window = int(n_days * 24 / model_dt)

vars_sd = np.load(output_f + file_n_1 + '_' + str(mean_height) + '_' + str(FWHM_gauss) + '_' + str(window) + '_' + var2 + ext + '.npy')[:all_ind,:]

md = smw              # number of days in running mean
sample_r = int(24 / model_dt)          # phases sampled every sample_r hours
cut_sdt = 0

# correct phase by 180 degrees when fitted amplitude is negative
inds_sdt_sd = np.where(vars_sd[:,1] < 0)
vars_sd_corr = np.copy(vars_sd)
vars_sd_corr[:,4][inds_sdt_sd] += np.pi
sdt_sd = tide_func.running_mean(np.abs(vars_sd[::sample_r,1]),md)

abs_sdt_phase = circmean(vars_sd_corr[::sample_r,4])

# meridional phase
v_midphase = 6.0 # hours at center of phase-plot (LTOM)
dmean_lotm_sdt_v = circmean(np.pi / 2 - abs_sdt_phase)
dsdt_offset_v = -(v_midphase - dmean_lotm_sdt_v / (2 * np.pi) * 12) / 12 * (2 * np.pi) 
dphase_sdt_v = vars_sd_corr[::sample_r,4]

for k in range(np.shape(dphase_sdt_v)[0]):
        dphase_sdt_v[k] = circmean(np.pi / 2 - dphase_sdt_v[k])
        
#% dope output 2

output_f = '/home/wim/Desktop/Projects/iSSW/data/SD_sampled/'
n_days = 4
window = int(n_days * 24 / model_dt)

vars_sd = np.load(output_f + file_n_2 + '_' + str(mean_height) + '_' + str(FWHM_gauss) + '_' + str(window) + '_' + var2 + ext + '.npy')[:all_ind,:]

md = smw              # number of days in running mean
sample_r = int(24 / model_dt)          # phases sampled every sample_r hours
cut_sdt = 0

# correct phase by 180 degrees when fitted amplitude is negative
inds_sdt_sd = np.where(vars_sd[:,1] < 0)
vars_sd_corr = np.copy(vars_sd)
vars_sd_corr[:,4][inds_sdt_sd] += np.pi
sdt_sd2 = tide_func.running_mean(np.abs(vars_sd[::sample_r,1]),md)

abs_sdt_phase = circmean(vars_sd_corr[::sample_r,4])

# meridional phase
v_midphase = 6.0 # hours at center of phase-plot (LTOM)
dmean_lotm_sdt_v = circmean(np.pi / 2 - abs_sdt_phase)
dsdt_offset_v = -(v_midphase - dmean_lotm_sdt_v / (2 * np.pi) * 12) / 12 * (2 * np.pi) 
dphase_sdt_v2 = vars_sd_corr[::sample_r,4]

for k in range(np.shape(dphase_sdt_v)[0]):
        dphase_sdt_v2[k] = circmean(np.pi / 2 - dphase_sdt_v2[k])
        
#% dope output 3

output_f = '/home/wim/Desktop/Projects/iSSW/data/SD_sampled/'
n_days = 4
window = int(n_days * 24 / model_dt)

vars_sd = np.load(output_f + file_n_3 + '_' + str(mean_height) + '_' + str(FWHM_gauss) + '_' + str(window) + '_' + var2 + ext + '.npy')[:all_ind,:]

md = smw              # number of days in running mean
sample_r = int(24 / model_dt)          # phases sampled every sample_r hours
cut_sdt = 0

# correct phase by 180 degrees when fitted amplitude is negative
inds_sdt_sd = np.where(vars_sd[:,1] < 0)
vars_sd_corr = np.copy(vars_sd)
vars_sd_corr[:,4][inds_sdt_sd] += np.pi
sdt_sd3 = tide_func.running_mean(np.abs(vars_sd[::sample_r,1]),md)

abs_sdt_phase = circmean(vars_sd_corr[::sample_r,4])

# meridional phase
v_midphase = 6.0 # hours at center of phase-plot (LTOM)
dmean_lotm_sdt_v = circmean(np.pi / 2 - abs_sdt_phase)
dsdt_offset_v = -(v_midphase - dmean_lotm_sdt_v / (2 * np.pi) * 12) / 12 * (2 * np.pi)
dphase_sdt_v3 = vars_sd_corr[::sample_r,4]

for k in range(np.shape(dphase_sdt_v)[0]):
        dphase_sdt_v3[k] = circmean(np.pi / 2 - dphase_sdt_v3[k])
               
#%%

xmin = 20
xmax = 65

err_step_amp = 475 # errorbar frequency amplitude
err_s_phase = 20   # errorbar frequency (sampled) phase

lss = 11.
ts = 14.4   # tick size
fs = 15   # laxes abel size
lw = 2.0   # line width amplitude curves
ms = 6.0    # marker size phase dots
tfs = 10.5  # text font size
alph = 1.0
alph_l = 0.5

c1 = 'red'
c2 = 'blue'

c4 = sns.set_hls_values(color = 'black', h = None, l = 0.25, s = None)
c3 = sns.set_hls_values(color = 'green', h = None, l = 0.25, s = None)
c2 = sns.set_hls_values(color = 'red', h = None, l = 0.35, s = None)
c1 = sns.set_hls_values(color = 'blue', h = None, l = 0.35, s = None)

fig_sample = 1
fig = plt.figure(figsize=(8.0,4))

plot_len = len(sdt_14_15)
plot_len_2 = np.shape(sdt_14_15[::fig_sample])[0] 

ax1 = fig.add_axes([0.0, 0.9, 0.8, 0.5],
                   xticklabels=[], xlim=(0,plot_len), ylim=(-2.25, 26.25))
# ax2 = fig.add_axes([0.0, 0.3, 0.8, 0.5],
#                    xlim=(0,plot_len_2),ylim=(-0.5, 2 * np.pi + 0.5))

xskim = 0
      
plotstep = 5
xlocs = np.arange(0,plot_len + 5,plotstep) 
xlabs = xlocs

xlocs_phase = np.arange(0,plot_len_2 + 10,plotstep / fig_sample) 

offset_a = 0
alpha_l = 0.09

ax1.plot(sdt_14_15,c=c1,alpha=alph,lw=lw,label='SD')
ax1.fill_between(np.arange(np.shape(sdt_14_15)[0]), 
                    (sdt_14_15 - fit_v_amp * 2),(sdt_14_15 + fit_v_amp * 2),alpha=alpha_l,color=c1)

ax1.plot(sdt_sd - offset_a,c=c2,ls='-',alpha=alph,lw=lw,label='PRISM-SD')
# ax1.plot(sdt_sd2 - offset_a,c=c3,ls='--',alpha=alph,lw=lw,label='OnlyThermal-SD')
ax1.plot(sdt_sd3,c=c3,ls='--',alpha=alph,lw=lw,label='LunarShift-SD')

if extended:
    sw1_amp = np.abs(varibs2[::soup_sample,inc_ind])
    ax1.plot(sw1_amp,c=c1,ls='--')
    ax1.plot(np.abs(vars_sd[::sample_r,inc_ind]),c=c2,ls='--')

# ax2.errorbar(np.arange(len(phase_sdt[::fig_sample])),phase_sdt[::fig_sample],yerr=fit_v_phas[::fig_sample] * 2,
#           fmt='x',c=c1,ms=ms,label=r'SD',capsize=1,elinewidth=1.25,capthick=0,errorevery=1)
# ax2.plot(dphase_sdt_v[::fig_sample],
#           '.',c=c2,ms=ms,label=r'PRISM-SD')
# ax2.plot(dphase_sdt_v2[::fig_sample],
#           '2',c=c3,ms=ms,label=r'OnlyThermal-SD')
# ax2.plot(dphase_sdt_v3[::fig_sample],
#           '+',c=c4,ms=ms,label=r'OnlyLunar-SD')

ax1.set_ylabel(r'$\mathregular{A}_{\mathregular{SW2}}$ (ms$^{-1}$) in U',fontsize=fs,labelpad=2)
ax1.xaxis.set_ticks(xlocs)
ax1.set_ylim(0,25)
ax1.tick_params(axis='y',labelsize=ts)
ax1.grid(True)
ax1.set_xlim(xmin,xmax)

# ax2.yaxis.set_ticks([0,np.pi / 2, np.pi, 3 * np.pi / 2, 2 * np.pi])
# ax2.yaxis.set_ticklabels([r'00:00',r'03:00',r'06:00',r'09:00',r'12:00'],size=ts)
# ax2.set_ylabel(r'$\mathregular{LTOM}_{\mathregular{SW2}}$ (hr) in U',fontsize=fs,labelpad=4)
ax1.set_xlabel(r'Days since Dec 1st 2012',fontsize=fs)
# ax1.xaxis.set_ticks(xlocs)
ax1.xaxis.set_ticklabels(xlocs,size=ts)
# ax2.grid(True)
# ax2.set_xlim(xmin / fig_sample,xmax / fig_sample)
# ax2.set_ylim(0,2*np.pi)

ax1.legend(ncol=4,loc=1,prop={'size':lss},markerscale=2,bbox_to_anchor=(1.001, 1.0),framealpha=0.6)
# ax2.legend(ncol=4,loc=2,prop={'size':lss},markerscale=2,bbox_to_anchor=(-0.001, 0.2),framealpha=0.6)

fw = 'bold' # text font weight
# ax1.text(12.5,25.,'(a)',fontsize=fs + 4,fontweight=fw)
# ax2.text(12.5,2 * np.pi + .5,'(b)',fontsize=fs + 4,fontweight=fw)


ssw_t = 41
ssw_t2 = 34.5
ssw_t3 = 53
ssw_c = 'black'
ssw_ls = '--'
lw_ssw = 2
ax1.vlines(ssw_t,0,100,color=ssw_c,ls=ssw_ls,lw=lw_ssw) # x, ymin, ymax
ax1.vlines(ssw_t2,0,100,color=ssw_c,ls=ssw_ls,lw=lw_ssw) # x, ymin, ymax
ax1.vlines(ssw_t3,0,100,color=ssw_c,ls=ssw_ls,lw=lw_ssw) # x, ymin, ymax

save_f = '/home/wim/Desktop/Projects/iSSW/paper/figures/'
# plt.savefig(save_f + 'superdarn_lunarshift.png',dpi=300, bbox_inches="tight")

#%%

print(file_n_1)

print(0.025 * 60)

print(mean_height,FWHM_gauss)