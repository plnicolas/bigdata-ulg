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
	nc.close()

classify,k = [],1
for line in open('Nodes_temp.txt','r').readlines()[:]:
	classify = [float(i) for i in line.split( )]
	temp = np.array(classify).reshape(len(lt),len(ln_shifted))
	# plot a map for air temperature spatial trends
	temp2 = np.zeros((len(temp),len(temp.T)))
	zero = np.where(ln_shifted==0)[0][0]
	for i in range(len(ln_shifted)):
		temp2[:,i] = temp[:,int((i - zero ) % len(ln_shifted))] # phase shift
	norm = mpl.colors.Normalize(vmin=-30,vmax=30)

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
			ax.add_collection(PatchCollection(patches, facecolor=cmap(norm(temp2[j,i])), edgecolor='none'))
	ann1 = ax.annotate("Mean Temp, N Hemisphere:" + "{:+.2f}".format(np.mean(temp2[0:30])), xy=(0.6, 0.07), xycoords='axes fraction')
	ann2 = ax.annotate("Mean Temp, S Hemisphere:" + "{:+.2f}".format(np.mean(temp2[31:])), xy=(0.6, 0.02), xycoords='axes fraction')
	fig.savefig("/CECI/home/ulg/mast/eivanov/Big_Data/Results/Case_%s.png" %(k), dpi=200, bbox_inches='tight'); k=k+1; ann1.remove(); ann2.remove()

