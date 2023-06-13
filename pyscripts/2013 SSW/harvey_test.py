#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 10 10:13:17 2021

@author: wim
"""

from netCDF4 import Dataset

folder = '/media/wim/NieuwVolume/WACCMX-DART/WACCMX DART 2009/'
file = 'WACCMX+DART_T_2009011000-2009011923.nc'

#%%

op = Dataset(folder + file)

print(op)
print(op['PRESSURE'][:])