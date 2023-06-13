import numpy as np
import model_funcs as mf
import tide_func as tf
from matplotlib import pyplot as plt
import seaborn as sns
from matplotlib import gridspec
from wrf import interplevel
import scipy.optimize
import os
from netCDF4 import Dataset
from datetime import datetime
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches

#%%

u_lat = 70
t_lat = 90
interp_levels = np.array([48000,40000,30000])

folder = '/home/wim/Desktop/Projects/iSSW/dope/output/'
file_n = 'test58.dat'

day_begin = 80
n_days = 81 

model_dt = 1 # output timestep (hr)

u_field, v_field, t_field, lats, lons = mf.uvt_read_altitude_grid(file_n,
                                                                  folder,
                                                                  day_begin,
                                                                  n_days,
                                                                  model_dt,
                                                                  interp_levels)

#%

dt_nav_1 = (datetime(2013,1,1,0) - datetime(2012,12,1,0)).days
dt_nav_2 = (datetime(2013,4,1,0) - datetime(2013,1,1,0)).days
dy_nav = 90

# first year
nav_u = np.zeros((dt_nav_1,
                  np.shape(interp_levels)[0],
                  dy_nav,
                  dy_nav))

nav_t = np.copy(nav_u)

year_1 = '2012'
year_2 = '2013'

f_base = '/media/wim/backupdisc/NAVGEM Netcdf/' + year_1
folder = f_base + '/' + year_1
months = ['12']

data = Dataset(f_base + '/' + year_1 + '12/ha_navgem_uvzt_' + year_1 + '1201.nc',mode='r')

lats_nav = np.array(data['lat'][::2],dtype=np.float64) # -89.3 to 89.3
lons_nav = np.array(data['lon'][::4],dtype=np.float64)

k = 0

for index, month in enumerate(months):
    directory = os.fsencode(folder + month + '/')
    for file in sorted(os.listdir(directory)):
        # import data
        filename = os.fsdecode(file)
        data = Dataset(folder + month + '/' + str(filename), mode='r')
        
        # daily means with 90 x 90 lon x lat points
        data_u = np.mean(data['u'][:,:,::2,::4],axis=0)
        data_t = np.mean(data['t'][:,:,::2,::4],axis=0)
        data_z = np.mean(data['z'][:,:,::2,::4],axis=0)
        
        # data_z = tf.geopotential_to_geometric(data_z)
        
        dz, dy, dx = np.shape(data_u)
        dq = np.shape(interp_levels)[0]
        
        alts = np.ones((dq,dy,dx))
        for q in range(dq):
            alts[q,:,:] = alts[q,:,:] * interp_levels[q]
         
        # interp to DOPE levels and latitudes        
        for y in range(dy):
            for x in range(dx):
                nav_u[k,:,y,x] = np.interp(alts[:,y,x],data_z[::-1,y,x],data_u[::-1,y,x])
                nav_t[k,:,y,x] = np.interp(alts[:,y,x],data_z[::-1,y,x],data_t[::-1,y,x])
                
        k += 1
        print('Constructing NAVGEM-HA mean day: ', k)
         
nav_u_1 = np.copy(nav_u)
nav_t_1 = np.copy(nav_t)
        
#% second year
nav_u = np.zeros((dt_nav_2,
                  np.shape(interp_levels)[0],
                  dy_nav,
                  dy_nav))

nav_t = np.copy(nav_u)

f_base = '/media/wim/backupdisc/NAVGEM Netcdf/' + year_2
folder = f_base + '/' + year_2
months = ['01','02','03']

k = 0

for index, month in enumerate(months):
    directory = os.fsencode(folder + month + '/')
    for file in sorted(os.listdir(directory)):
        # import data
        filename = os.fsdecode(file)
        data = Dataset(folder + month + '/' + str(filename), mode='r')
        
        # daily means with 90 x 90 lon x lat points
        data_u = np.mean(data['u'][:,:,::2,::4],axis=0)
        data_t = np.mean(data['t'][:,:,::2,::4],axis=0)
        data_z = np.mean(data['z'][:,:,::2,::4],axis=0)
        
        # data_z = tf.geopotential_to_geometric(data_z)
        
        dz, dy, dx = np.shape(data_u)
        dq = np.shape(interp_levels)[0]
        
        alts = np.ones((dq,dy,dx))
        for q in range(dq):
            alts[q,:,:] = alts[q,:,:] * interp_levels[q]
         
        # interp to DOPE levels and latitudes        
        for y in range(dy):
            for x in range(dx):
                nav_u[k,:,y,x] = np.interp(alts[:,y,x],data_z[::-1,y,x],data_u[::-1,y,x])
                nav_t[k,:,y,x] = np.interp(alts[:,y,x],data_z[::-1,y,x],data_t[::-1,y,x])
                
        k += 1
        print('Constructing NAVGEM-HA mean day: ', k)
        
nav_u_2 = np.copy(nav_u)
nav_t_2 = np.copy(nav_t)

nav_u = np.concatenate((nav_u_1,nav_u_2),axis=0)
nav_t = np.concatenate((nav_t_1,nav_t_2),axis=0)

#%%

def pw_fit(x,
           a_s1,
           a_s2,
           a_s3,
           p_s1,
           p_s2,
           p_s3):
    
    return a_s1 * np.sin(x + p_s1) + a_s2 * np.sin(x * 2 + p_s2) + a_s3 * np.sin(x * 3 + p_s3)
    
nav_min = 50
nav_max = 70

u_mini = np.argmin(np.abs(lats_nav - nav_min))
u_maxi = np.argmin(np.abs(lats_nav - nav_max))

s_mini = np.argmin(np.abs(lats - nav_min))
s_maxi = np.argmin(np.abs(lats - nav_max))

nav_av = np.mean(nav_u[:,0,u_mini:u_maxi,:],axis=1)
sig_av = np.mean(u_field[:,0,s_mini:s_maxi,:],axis=1)

#%

dt, dx = np.shape(nav_av)
sig_amps = np.zeros((n_days,3))
nav_amps = np.zeros((n_days,3))

d_mean = 4

for d in range(81):
    sig_mean = np.mean(sig_av[(d * 24):(d * 24 + d_mean * 24),:],axis=0)
    
    popt, pcov = scipy.optimize.curve_fit(pw_fit, (lons / 180 * np.pi), sig_mean)
    sig_amps[d,:] = np.abs(popt[:3])
    
    nav_avg = np.mean(nav_av[d:(d + d_mean),:],axis=0)
    popt, pcov = scipy.optimize.curve_fit(pw_fit, (lons_nav / 180 * np.pi), nav_avg)
    nav_amps[d,:] = np.abs(popt[:3])     
    
nav_amps = np.roll(nav_amps,int(d_mean / 2),axis=0)
sig_amps = np.roll(sig_amps,int(d_mean / 2),axis=0)

#%

# plt.plot(sig_amps[:,1])
# plt.plot(nav_amps[:,1])

# plt.plot(sig_amps[:,0])
# plt.plot(nav_amps[:,0])

#%

u_navind = np.argmin(np.abs(lats_nav - u_lat))
t_navind = np.argmin(np.abs(lats_nav - t_lat))

u_latind = np.argmin(np.abs(lats - u_lat))
t_latind = np.argmin(np.abs(lats - t_lat))

zm_u = np.mean(u_field,axis=3)
zm_t = np.mean(t_field,axis=3)

d_start = 20 
d_end = 65 

day_plt = 40
alt_ind = 0
linefs = 1
l_fs = 9

minn = -80
zlvls = np.arange(minn,-minn + 1,20)
vmax = 65
vmin = -vmax

#%

c1 = sns.set_hls_values(color = 'blue', h = None, l = 0.35, s = None)
c2 = sns.set_hls_values(color = 'red', h = None, l = 0.40, s = None)

lat_ticks = np.arange(-75,76,15)
lon_ticks = np.arange(30,351,60)

ts = 12.0
fs = 11
lw = 1.5

t_ax = np.arange((d_end - d_start) * 24) / 24 + d_start
t_nav = np.arange(121)
xticks = np.arange(d_start,d_end + 1,5)
xtick_lab = xticks
# xticks = [20,25,31,36,41,46,51,56,62,65]
# xtick_lab = ['Dec-20', '25','Jan-01','5','10','15','20','25','Feb-01','']

umin = -80
umax = 81
tmin = 195
tmax = 281
uticks = np.arange(-70,80,20)
tticks = np.arange(200,281,15)

uticks_pw = np.arange(0,71,10)

fig, (ax,ax2) = plt.subplots(2,1,figsize=(5,5.0),gridspec_kw={'height_ratios': [1,1.0]})
fig.tight_layout()
plt.subplots_adjust(hspace = 0.35)

ax.plot(t_ax,zm_u[(d_start * 24):(d_end * 24),0,u_latind],c=c1,lw=lw)
ax.plot(t_nav[d_start:d_end],np.mean(nav_u[d_start:d_end,0,u_navind,:],axis=1),c=c1,lw=lw,ls='--')

# ax.plot(t_ax,np.mean(nav_u[(d_start * 24):(d_end * 24),0,u_navind,:],axis=1),c=c1,lw=lw)
# ax.plot(t_nav[d_start:d_end],np.mean(nav_t[d_start:d_end,0,u_navind,:],axis=1),c=c1,lw=lw)

ax.set_ylabel(r'$\bar{U}$ (ms$^{-1}$) at 70$^{\circ}$N / 48 km',fontsize=fs)
ax.set_xlim(d_start,d_end)
ax.set_ylim(umin,umax - 1)
ax.set_xticks(xticks)
ax.set_xticklabels(xtick_lab,size=ts)
ax.set_yticks(uticks)
ax.set_yticklabels(uticks,color=c1,size=ts)
# ax.set_xlabel('Days since Dec 1st 2012',fontsize=fs)

ax3 = ax.twinx()
ax3.plot(t_ax,zm_t[(d_start * 24):(d_end * 24),1,t_latind],c=c2,lw=lw,label='PRISM')
ax3.plot(t_nav[d_start:d_end],np.mean(nav_t[d_start:d_end,1,t_navind,:],axis=1),c=c2,lw=lw,ls='--',label='NAVGEM-HA')
ax3.set_ylabel(r'$\bar{T}$ (K) at 90$^{\circ}$N / 40 km',fontsize=fs)
ax3.set_yticks(tticks)
ax3.set_yticklabels(tticks,color=c2,size=ts)
ax3.set_ylim(tmin,tmax - 1)

handles, labels = plt.gca().get_legend_handles_labels()

line1 = Line2D([0], [0], label='PRISM', color=c1)
line2 = Line2D([0], [0], label='NAVGEM-HA', color=c1,ls='--')

handles.extend([line1,line2])

ax3.legend(handles=handles,ncol=1,prop={'size': 7.8},loc=1,bbox_to_anchor=(0.8, 1.0))

pw1 = u'#1f77b4'
pw2 = u'#ff7f0e'
pw3 = 'black'

ax2.set_ylabel(r'Amp (ms$^{-1}$) in U at 50-70$^{\circ}$N / 48 km',fontsize=fs)
ax2.plot(sig_amps[:,0],color=pw1,label='PRISM PW1')
ax2.plot(nav_amps[:,0],color=pw1,ls='--',label='NAVGEM-HA PW1')
ax2.plot(sig_amps[:,1],color=pw2,label='PRISM PW2')
ax2.plot(nav_amps[:,1],color=pw2,ls='--',label='NAVGEM-HA PW2')
# ax2.plot(sig_amps[:,2],color=pw3,label='PRISM PW3')
# ax2.plot(nav_amps[:,2],color=pw3,ls='--',label='NAVGEM-HA PW3')
ax2.legend(prop={'size': 7.8})
ax2.set_xlim(d_start,d_end)
ax2.set_ylim(0,70 - 1)
ax2.set_xticks(xticks)
ax2.set_xticklabels(xtick_lab,size=ts)
ax2.set_yticks(uticks_pw)
ax2.set_yticklabels(uticks_pw,color='black',size=ts)
ax2.set_xlabel('Days since Dec 1st 2012',fontsize=fs)

# ct = ax2.contourf(lons,lats,u_field[(day_plt * 24),alt_ind,:,:],zlvls,vmax=vmax,vmin=vmin)
# cs = ax2.contour(lons,lats,u_field[(day_plt * 24),alt_ind,:,:], [0],linewidths=linefs,colors='white')
# ax2.clabel(cs, [0], inline=1, fontsize=l_fs,fmt='%1.0f')
# ax2.set_yticks(lat_ticks)
# ax2.set_yticklabels(lat_ticks,size=ts)
# ax2.set_ylim(-5,80)
# ax2.set_xticks(lon_ticks)
# ax2.set_xticklabels(lon_ticks,size=ts)
# ax2.set_xlabel(r'Longitude ($^{\circ}$)',size=fs + 0.5)
# ax2.set_ylabel(r'Latitude ($^{\circ}$)',size=fs + 0.5)

# cb = fig.colorbar(ct, ax=ax2,label='test',pad=0.035,fraction=0.035)
# cb.set_label(label='U (ms$^{-1}$) at 48 km / Jan 11th 2013', size=fs)
# cb.ax.tick_params(labelsize=ts)

fw = 'bold'
ax3.text(17.5,285,'(a)',fontsize=fs + 5,fontweight=fw)
ax2.text(17.5,75,'(b)',fontsize=fs + 5,fontweight=fw)
 
ssw_c = 'black'
ssw_ls = '--'
lw_ssw = 1.25
ax3.vlines(41,0,300,color=ssw_c,ls=ssw_ls,lw=lw_ssw) # x, ymin, ymax
ax3.vlines(34.5,0,300,color=ssw_c,ls=ssw_ls,lw=lw_ssw) # x, ymin, ymax
ax3.vlines(53,0,300,color=ssw_c,ls=ssw_ls,lw=lw_ssw) # x, ymin, ymax

ax2.vlines(41,0,300,color=ssw_c,ls=ssw_ls,lw=lw_ssw) # x, ymin, ymax
ax2.vlines(34.5,0,300,color=ssw_c,ls=ssw_ls,lw=lw_ssw) # x, ymin, ymax
ax2.vlines(53,0,300,color=ssw_c,ls=ssw_ls,lw=lw_ssw) # x, ymin, ymax

save_f = '/home/wim/Desktop/Projects/iSSW/paper/figures/'
# plt.savefig(save_f + 'ssw_model_validation_pw.png',dpi=300, bbox_inches="tight")

#%%














#%%

test = np.arange(18 * 31).reshape(18,31)
test = np.array(test,dtype='object')

x = np.arange(18,dtype='object')
y = np.arange(31,dtype='object')

#%%

plt.contourf(y,x,test)
plt.colorbar()

#%%

print(np.shape(test))

#%%


amplitudes = np.zeros((18,31))

for day in range(31):
        
    amplitudes[:,day] = your_amplitudes

#%%

np.save()

#%%

































