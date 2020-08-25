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

for i in d:
	nc = Dataset(path + '/' + i, 'r', format='NETCDF4')
	lt = nc.variables['latitude'][:]
	ln = np.roll(nc.variables['longitude'][:],60)
	ln_shifted = np.array([el-360 if  el>179.9 else el for el in ln])			# convert logitudes from 0-360 to -180-+180
	t = nc.variables['time'][:]
	for k in range(len(t)):
		ind = (lo2(int(t[k])).year-1979)*12+(lo2(int(t[k])).month-1)
		ta[ind] = nc.variables['t2m'][k]-273.15
	nc.close()

# interpolate them into grid cell centers
# longitude curviture is not taken into account, assumed that it is constant and equal to 111 km along its lenght
ta = np.roll(ta,60,axis=2)

cities = ['Brussels', 'New York', 'Los Angeles', 'Brasilia', 'Yellowknife', 'McMurdo', 'Johannesburg', 'Kinshasa', 'Valletta', 'Moscow', 'Tokio', 'Yakutsk', 'Singapore', 'Karachi', 'Melbourne']
city_lat = [50+51/60, 40+42/60, 34+3/60, -15-47/60, 62+26/60, -77-50/60, -26-12/60, -4-19/60, 35+53/60, 55+45/60, 35+41/60, 62+2/60, 1+17/60, 24+51/60, -37-48/60]
city_lon = [4+21/60, -74+00/60, -118-15/60, -47-52/60, -114-23/60, 166+40/60, 28+2/60, 15+19/20, 14+30/60, 37+37/60, 139+41/60, 129+44/60, 103+50/60, 67+0/60, 144+5]
year = np.arange(39)+1979

idx_lat = np.zeros((len(cities)))
idx_lon = np.zeros((len(cities)))
for i in range(len(cities)):
	idx_lat[i] = (np.abs(lt-city_lat[i])).argmin()
	idx_lon[i] = (np.abs(ln_shifted-city_lon[i])).argmin()

temp = np.zeros((len(cities),12))
for i in range(len(cities)):
	temp_add = [[] for x in range(12)]
	for j in range(len(ta)):
		temp_add[j%12].append(ta[j,int(idx_lat[i]),int(idx_lon[i])])
	for k in range(12):
		z = np.polyfit(year, temp_add[k], 1)
		pa = np.poly1d(z)
		temp[i,k] = pa(2017)-pa(1979)

months = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
colors = ['cyan','cyan','yellow','yellow','yellow','green','green','green','red','red','red','cyan']
for i in range(len(cities)):
	fig = plt.figure(figsize=(8,4))
	ax1 = fig.add_subplot(111)
	opacity = 0.4; error_config = {'ecolor': '0.3'}
	ax1.bar(range(len(months)), list(map(float,temp[i])),color=colors)
	plt.xticks(range(len(months)), months)
	plt.ylim(-2.5,2.5)
	ax1.set_ylabel('Trend [C]', fontweight='bold')
	ax1.set_title(cities[i])
	fig.savefig('/CECI/home/ulg/mast/eivanov/Big_Data/Monthly/%s.png' %(cities[i]), dpi=100, bbox_inches='tight')
