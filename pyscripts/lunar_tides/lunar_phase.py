#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  5 12:34:22 2021

@author: wim
"""

import pylunar as pl

phase = pl.MoonInfo((0,0,0),(0,0,0))
phase.update((2014,1,1,0,0,0))

print(phase.age())

#%%