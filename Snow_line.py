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
lo3 = lambda x: datetime(1979,1,1,0,0,0) + timedelta(hours=x)

# 1) DEFINE THE PATH
path = '/scratch/ulg/mast/eivanov/Files_for_Big_Data/'
d = os.listdir(path)

# Read netcdf data files and build an array of air and soil temeprature
ta = np.zeros((29*365*2+10*366*2,61,120))
k = 0

for i in sorted(d):
	nc = Dataset(path + '/' + i, 'r', format='NETCDF4')
	lt = nc.variables['latitude'][:]
	ln = np.roll(nc.variables['longitude'][:],60)
	ln_shifted = np.array([el-360 if  el>179.9 else el for el in ln])			# convert logitudes from 0-360 to -180-+180
	t = nc.variables['time'][:]
	for m in range(len(t)):
		ta[k] = nc.variables['sf'][m]
		k = k+1
	nc.close()

ta = ta*1000 # m to mm

"""
# visualize cliamtic line shift
ta_before = np.zeros((len(lt),len(ln)))
for i in range(0,(1997-1979)*2*365+5*2):
	ta_before += ta[i]
ta_before = ta_before / (i+1)

k = 0
ta_after = np.zeros((len(lt),len(ln)))
for i in range((1997-1979)*2*365+5*2,(2018-1997)*2*365+5*2):
	ta_after += ta[i]; k = k+1
ta_after = ta_after / k


diff = (ta_after - ta_before)*2 # m to mm; hakf-day to day
arr = np.array([diff.min(),diff.max()])

norm = mpl.colors.Normalize(vmin=-1*abs(arr).max(),vmax=abs(arr).max())


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
for i in range(len(ln_shifted)):
	for j in range(len(lt)-1):
		patches = []
		if x[i] != 177:
			polygon = np.array([(x[i],y[j]),(x[i+1],y[j]),(x[i+1],y[j+1]),(x[i],y[j+1])])
		if x[i] == 177:
			polygon = np.array([(x[i],y[j]),(180,y[j]),(180,y[j+1]),(x[i],y[j+1])])
		patches.append(Polygon(polygon))
		try:
			ax.add_collection(PatchCollection(patches, facecolor=cmap(norm(diff[j,i+60])), edgecolor='none'))
		except:
			ax.add_collection(PatchCollection(patches, facecolor=cmap(norm(diff[j,-60+i])), edgecolor='none'))
fig.savefig("/CECI/home/ulg/mast/eivanov/Big_Data/Results/Snowfall_difference.png", dpi=200, bbox_inches='tight')
"""

index = (1997-1979)*2*365+5*2
x_before = [[] for i in range(12)]
for i in range(0,index):
	x_before[int(lo3(int(i*12)).month)-1].append(ta[i])

x_after = [[] for i in range(12)]
for i in range(index,len(ta)):
	x_after[int(lo3(int(i*12)).month)-1].append(ta[i])

x_total = [[] for i in range(12)]
for i in range(len(ta)):
	x_total[int(lo3(int(i*12)).month)-1].append(ta[i])


def count_event(pp,rr):
	count_before,count_after = 0,0
	for i in range(len(ta[0:index])):
		if ta[i,pp,rr] > var99_before[mm[i],pp,rr]:
		#if ta[i,pp,rr] > var99_total[mm[i],pp,rr]:
			count_before = count_before + 1
	for i in range(len(ta[index:])):
		if ta[index + i,pp,rr] > var99_after[mm[i+index],pp,rr]:
		#if ta[index + i,pp,rr] > var99_total[mm[i+index],pp,rr]:
			count_after = count_after + 1
	return count_before,count_after

mm = []
for i in range(len(ta)):
	mm.append(int(lo3(int(i*12)).month)-1)

var99_before,var99_after,var99_total = np.zeros((12,len(lt),len(ln))),np.zeros((12,len(lt),len(ln))),np.zeros((12,len(lt),len(ln)))
for i in range(12):
	mean_before = np.mean(np.array(x_before[i]),axis=0)
	mean_after = np.mean(np.array(x_after[i]),axis=0)
	mean = np.mean(np.array(x_total[i]),axis=0)
	stdev_before = np.std(np.array(x_before[i]),axis=0)
	stdev_after = np.std(np.array(x_after[i]),axis=0)
	stdev = np.std(np.array(x_total[i]),axis=0)
	var99_before[i] = mean_before + ( stdev_before * 3 )
	var99_after[i] = mean_after + ( stdev_after * 3 )
	var99_total[i] = mean_after + ( stdev * 3 )

ta_before, ta_after = np.zeros((len(lt),len(ln))),np.zeros((len(lt),len(ln)))
for i in range(len(lt)):
	for j in range(len(ln)):
		ta_before[i,j],ta_after[i,j] = count_event(i,j)
		print(i,j)
ta_after = ta_after / ((2018-1997)/(1997-1979)) # eliminate difference between lengths of time-series

ta_before = ta_before / 12
ta_after = ta_after / 12

extreme_diff = ta_after - ta_before

norm = mpl.colors.Normalize(vmin=-1*abs(extreme_diff).max(),vmax=abs(extreme_diff).max())


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
for i in range(len(ln_shifted)):
	for j in range(len(lt)-1):
		patches = []
		if x[i] != 177:
			polygon = np.array([(x[i],y[j]),(x[i+1],y[j]),(x[i+1],y[j+1]),(x[i],y[j+1])])
		if x[i] == 177:
			polygon = np.array([(x[i],y[j]),(180,y[j]),(180,y[j+1]),(x[i],y[j+1])])
		patches.append(Polygon(polygon))
		try:
			ax.add_collection(PatchCollection(patches, facecolor=cmap(norm(extreme_diff[j,i+60])), edgecolor='none'))
		except:
			ax.add_collection(PatchCollection(patches, facecolor=cmap(norm(extreme_diff[j,-60+i])), edgecolor='none'))
fig.savefig("/CECI/home/ulg/mast/eivanov/Big_Data/Results/Snowfall_extreme.png", dpi=200, bbox_inches='tight')

#ta_before = np.zeros((len(lt),len(ln)))
#for i in range(0,(1997-1979)*2*365+5*2):
#	ta[i][ta[i]>0] = 1
#	ta[i][ta[i]<0] = 0
#	ta_before += ta[i]
#ta_before[ta_before>30]=30

#ta_after = np.zeros((len(lt),len(ln)))
#for i in range((1997-1979)*2*365+5*2,(2018-1997)*2*365+5*2):
#	ta[i][ta[i]>0] = 1
#	ta[i][ta[i]<0] = 0
#	ta_after += ta[i]
#ta_after = ta_after / ((2018-1997)-(1997-1979))
#ta_after[ta_after>30]=30

