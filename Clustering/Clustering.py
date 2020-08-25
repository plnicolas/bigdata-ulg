#! /usr/bin/env python
# -*- coding: utf-8 -*-
# created by Evgeny Ivanov for Big Data course at ULg to prove that, in fact, the climate change is real, 17/11/2018
# execution time is ~ 1 min

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import os
from datetime import date, datetime, timedelta
from use_func import *
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable

lo2 = lambda x: datetime(1900,1,1,0,0,0) + timedelta(hours=x) # converts ESMWF time to our units

# 1) DEFINE THE PATH
path = '/CECI/home/ulg/mast/eivanov/Big_Data/Big_Data_Files_Air_Soil_Temp'
d = os.listdir(path)

# Read netcdf data files and build an array of air and soil temeprature
ta = np.zeros((39*12,61,120))
ts = np.zeros((39*12,61,120))

for i in d:
	nc = Dataset(path + '/' + i, 'r', format='NETCDF4')
	lt = nc.variables['latitude'][:]
	ln = np.roll(nc.variables['longitude'][:],60)
	ln_shifted = np.array([el-360 if  el>179.9 else el for el in ln])			# convert logitudes from 0-360 to -180-+180
	t = nc.variables['time'][:]
	for k in range(len(t)):
		ind = (lo2(int(t[k])).year-1979)*12+(lo2(int(t[k])).month-1)
		ta[ind] = nc.variables['t2m'][k]-273.15
		ts[ind] = nc.variables['stl4'][k]-273.15
	print(lo2(int(nc.variables['longitude'][0])))
	nc.close()

fil = open('Array_temp.txt', 'w')
for i in range(len(ta)):
	taa = np.round((ta[i,:,:]).flatten(),1)
	taa = [str(i) for i in taa]
	fil.write('  '.join(taa) + '\n')
fil.close()

lat_bel,lon_bel = 50+51/60, 4+21/60
fil = open('Array_temp_Belgium.txt', 'w')
for i in range(len(ta)):
	taa = np.round((ta[i,13,1:3]).flatten(),1)
	taa = [str(i) for i in taa]
	fil.write('  '.join(taa) + '\n')
fil.close()
