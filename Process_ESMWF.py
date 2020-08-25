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
	nc.close()


# interpolate them into grid cell centers
# longitude curviture is not taken into account, assumed that it is constant and equal to 111 km along its lenght
ta = np.roll(ta,60,axis=2)
ts = np.roll(ts,60,axis=2)

ta_center = np.zeros((39*12,61-1,120))
ts_center = np.zeros((39*12,61-1,120))
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
			ts_center[k,i,j] = (ts[k,i,j]+ts[k,i+1,j]+ts[k,i,j+1]+ts[k,i+1,j+1])/4

for i in range(len(lt)-1):
	for k in range(len(ta)):
		ta_center[k,i,-1] = (ta[k,i,-1]+ta[k,i+1,-1]+ta[k,i,0]+ta[k,i+1,0])/4
		ts_center[k,i,-1] = (ts[k,i,-1]+ts[k,i+1,-1]+ts[k,i,0]+ts[k,i+1,0])/4

airT = np.zeros((len(ta)))
soilT = np.zeros((len(ts)))

# calculate weighted global temperature
for i in range(len(ta)):
	airT[i] = np.sum(ta_center[i]*area)/np.sum(area)
	soilT[i] = np.sum(ts_center[i]*area)/np.sum(area)

air_year = []
soil_year = []
# calculate annual temperature
# difference in month lenghts are not taken into account assumed that all months are equal
for i in range(int(len(ta)/12)):
	air_year.append(np.sum(airT[i*12 : i*12 +12])/12)
	soil_year.append(np.sum(soilT[i*12 : i*12 +12])/12)

y = np.arange(39)+1979

# Plot the trends for the global temperature
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(y,soil_year,'black',label='soil')
ax.plot(y,air_year,'blue',label='air')
ax.set_xlim(1979,2017)
z_air = np.polyfit(y, air_year, 1)
pa = np.poly1d(z_air)
ax.plot(y,pa(y),"--",color='blue')
z_soil = np.polyfit(y, soil_year, 1)
ps = np.poly1d(z_soil)
ax.plot(y,ps(y),"--",color='black')
ax.text(2005,min(air_year),"y=%.6fx+%.6f" %(z_air[0],z_air[1]),color='blue', fontsize=7)
ax.text(2005,min(air_year)+0.1,"y=%.6fx+%.6f" %(z_soil[0],z_soil[1]),color='black', fontsize=7)
plt.legend(loc=2,prop={'size': 8})
plt.ylabel('Temperature [C]')
fig.savefig("/CECI/home/ulg/mast/eivanov/Big_Data/Results/Annual_air_and_soil_temperature_dynamics.png", dpi=200, bbox_inches='tight')


# calculate annual temperatures for each grid cell
ta_cell = np.zeros((39,61-1,120))
ts_cell = np.zeros((39,61-1,120))
for i in range(len(lt)-1):
	for j in range(len(ln)):
		for k in range(len(ta_cell)):
			ta_cell[k,i,j] = np.sum(ta_center[k*12 : k*12 +12,i,j])/12
			ts_cell[k,i,j] = np.sum(ts_center[k*12 : k*12 +12,i,j])/12



trend_air = np.zeros((61-1,120))
trend_soil = np.zeros((61-1,120))
# now calculate a trend for each grid cell, and based on this trend calculate a difference between the year 2017 and the year 1979 
for i in range(len(lt)-1):
	for j in range(len(ln)):
		z_air = np.polyfit(y, ta_cell[:,i,j], 1)
		z_air_pa = np.poly1d(z_air)
		trend_air[i,j] = z_air_pa(2017)-z_air_pa(1979)
		z_soil = np.polyfit(y, ts_cell[:,i,j], 1)
		z_soil_pa = np.poly1d(z_soil)
		trend_soil[i,j] = z_soil_pa(2017)-z_soil_pa(1979)

# plot a map for air temperature spatial trends

norm = mpl.colors.Normalize(vmin=-2,vmax=2)

fig = plt.figure(figsize=(10, 8))
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
m1 = Basemap(projection='cyl', ax=ax)
m1.drawcoastlines()
m1.drawmapboundary(fill_color='aqua')
m1.drawcountries()
#m1.fillcontinents(color='#ddaa66',lake_color='#9999FF')
cax = make_axes_locatable(ax).append_axes("right", size=0.4, pad=0.15)

cmap = mpl.cm.get_cmap('coolwarm')
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
		ax.add_collection(PatchCollection(patches, facecolor=cmap(norm(trend_air[j,i])), edgecolor='none'))
fig.savefig("/CECI/home/ulg/mast/eivanov/Big_Data/Results/Spatial_trend_air_temperature.png", dpi=200, bbox_inches='tight')

# plot a map for soil temperature spatial trends

cmap = mpl.cm.get_cmap('bwr')
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
		ax.add_collection(PatchCollection(patches, facecolor=cmap(norm(trend_soil[j,i])), edgecolor='none'))
fig.savefig("/CECI/home/ulg/mast/eivanov/Big_Data/Results/Spatial_trend_soil_temperature.png", dpi=200, bbox_inches='tight')
