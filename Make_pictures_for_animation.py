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

ta_center = np.zeros((39*12,61-1,120))
area = np.zeros((61-1,120))

for i in range(len(lt)-1):
	for j in range(len(ln)-1):
		area[i,j] = get_area_from_geo_coord(lt[i],lt[i+1],ln_shifted[j],ln_shifted[j+1])

for i in range(len(lt)-1):
	area[i,-1] = get_area_from_geo_coord(lt[i],lt[i+1],177,180)

for i in range(len(lt)-1):
	for j in range(len(ln)-1):
		for k in range(len(ta)):
			ta_center[k,i,j] = (ta[k,i,j]+ta[k,i+1,j]+ta[k,i,j+1]+ta[k,i+1,j+1])/4

for i in range(len(lt)-1):
	for k in range(len(ta)):
		ta_center[k,i,-1] = (ta[k,i,-1]+ta[k,i+1,-1]+ta[k,i,0]+ta[k,i+1,0])/4

# calculate annual temperatures for each grid cell
ta_cell = np.zeros((39,61-1,120))
for i in range(len(lt)-1):
	for j in range(len(ln)):
		for k in range(len(ta_cell)):
			ta_cell[k,i,j] = np.sum(ta_center[k*12 : k*12 +12,i,j])/12

norm = mpl.colors.Normalize(vmin=-2,vmax=2)
year = np.arange(39)+1979

airT = np.zeros((len(ta)))
# calculate weighted global temperature
for i in range(len(ta)):
	airT[i] = np.sum(ta_center[i]*area)/np.sum(area)

air_year = []
# calculate annual temperature
# difference in month lenghts are not taken into account assumed that all months are equal
for i in range(int(len(ta)/12)):
	air_year.append(np.sum(airT[i*12 : i*12 +12])/12)

# calculate the global trend
zz = np.polyfit(year, air_year, 1)
zz_pa = np.poly1d(zz)


trend_air = np.zeros((39,61-1,120))
# now calculate a trend for each grid cell, and based on this trend calculate a difference between the year 2017 and the year 1979 
for i in range(len(lt)-1):
	for j in range(len(ln)):
		z_air = np.polyfit(year, ta_cell[:,i,j], 1)
		z_air_pa = np.poly1d(z_air)
		for k in range(len(trend_air)):
			trend_air[k,i,j] = z_air_pa(1979+k) - z_air_pa(1979) - (zz_pa(1979+k)) + (zz_pa(1979))
	#print(trend_air[i])



for k in range(len(trend_air)):
	if k%5==0:
		fig = plt.figure(figsize=(8, 6))
		ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
		m1 = Basemap(projection='cyl', ax=ax)
		m1.drawcoastlines()
		m1.drawmapboundary(fill_color='aqua')
		m1.drawcountries()
		#m1.fillcontinents(color='#ddaa66',lake_color='#9999FF')
		cax = make_axes_locatable(ax).append_axes("right", size=0.4, pad=0.15)
		cmap = mpl.cm.get_cmap('seismic')
		cb = mpl.colorbar.ColorbarBase(cax,cmap=cmap,norm = norm,orientation = 'vertical')
		x, y = m1(ln_shifted,lt)	
	for i in range(len(ln)):
		for j in range(len(lt)-1):
			patches = []
			if x[i] != 177:
				polygon = np.array([(x[i],y[j]),(x[i+1],y[j]),(x[i+1],y[j+1]),(x[i],y[j+1])])
			if x[i] == 177:
				polygon = np.array([(x[i],y[j]),(180,y[j]),(180,y[j+1]),(x[i],y[j+1])])
			patches.append(Polygon(polygon))
			ax.add_collection(PatchCollection(patches, facecolor=cmap(norm(trend_air[k,j,i])), edgecolor='none'))
	ann = ax.annotate(str(1979+k) + "  Temp:" + "{:+.2f}".format(zz_pa(year[k])-zz_pa(year[0])), xy=(0.7, 0.05), xycoords='axes fraction'); print(1979+k)
	fig.savefig("/CECI/home/ulg/mast/eivanov/Big_Data/Pics_for_animation/Spatial_trend_air_temperature_%s.png" %(1979+k), dpi=150, bbox_inches='tight')
	ann.remove()
