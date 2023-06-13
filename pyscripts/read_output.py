#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 19:43:09 2020

@author: wim
"""

#%% single timestep full plot

import numpy as np
import model_funcs as mf
from matplotlib import pyplot as plt

folder = '/home/wim/Desktop/Projects/iSSW/dope/output/'
file = 'test2.dat'

record = 1200

uf, vf, tf, lats, lons, spres, geob, tracers = mf.read_uv_grid(file,
                                                               folder,
                                                               record)

# zonal mean and flip altitude
field_1 = np.mean(uf[::-1,:,:],axis=2)

#%

plt.pcolormesh(field_1)
plt.colorbar()

#%% read and save model output

import numpy as np
import model_funcs as mf
from matplotlib import pyplot as plt

folder = '/home/wim/Desktop/Projects/iSSW/dope/output/'
file_n = 'test3.dat'

day_begin = 80
n_days = 81

model_dt = 1 # output timestep (hr)
interp_levels = np.arange(130000,69999,-2500)

u_field, v_field, t_field, lats, lons = mf.uvt_read_altitude_grid(file_n,
                                                                  folder,
                                                                  day_begin,
                                                                  n_days,
                                                                  model_dt,
                                                                  interp_levels)

#%

save_f = '/home/wim/Desktop/Projects/iSSW/data/sigprim_fits/'

np.save(save_f + file_n + '_u.npy', u_field)
# np.save(save_f + file_n + '_v.npy', v_field)
# np.save(save_f + file_n + '_t.npy', t_field)
np.save(save_f + file_n + '_lat.npy',lats)
np.save(save_f + file_n + '_lon.npy',lons)
np.save(save_f + file_n + '_height_array.npy', interp_levels)

#%%


print(360 / 24)

















#%%

folder = '/home/wim/Desktop/Projects/iSSW/dope/output/'
file = 'test1.dat'

f = FortranFile(folder + file,'r')

# read in grid dimensions / metadata gunk
md = f.read_record(dtype=np.int32)
(k1, k2, m, n, l) = (md[0], md[1], md[2], md[3], md[4])

slv = f.read_record(dtype=np.float32)
lat = f.read_record(dtype=np.float32)
lon = np.linspace(0,360,k1,endpoint=False)
zstdlvl = f.read_record(dtype=np.float32)

#%%

u = f.read_record(dtype=np.float32).reshape(k1,k2,l,order='F')
v = f.read_record(dtype=np.float32).reshape(k1,k2,l,order='F')
t = f.read_record(dtype=np.float32).reshape(k1,k2,l,order='F')
p = f.read_record(dtype=np.float32).reshape(k1,k2,order='F')
geob = f.read_record(dtype=np.float32).reshape(k1,k2,order='F')

#%%

u = np.swapaxes(u,0,2)
temp = np.mean(u,axis=2)
plt.pcolormesh(temp[:,:])
plt.colorbar()
