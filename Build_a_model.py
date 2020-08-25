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

z = np.polyfit(y, air_year, 2)
f = np.poly1d(z)
x_new = np.arange(2017,2031)
y_old = f(y)
y_new = f(x_new)

# simple prediction (based on polynomial of 2 order)
fig = plt.figure()
ax = fig.add_subplot(111)
ax.scatter(y,air_year,c='black',label='Observations')
ax.plot(y,list(y_old),c='blue',linestyle='-',label='2nd order polinomial approximation')
ax.plot(x_new,list(y_new),c='blue',linestyle='--',label='Prediction based on 2nd order polinomial')
ax.set_xlim(1979,2030)
ax.set_ylim(13.5,15)
plt.ylabel('Temperature [C]')



# The zero-dimensional model of the Earth
def co2_over_time(year):
	# emissivity of the carbon dioxide for its current concentration in the atmosphere
	# http://www.biocab.org/ECO2.pdf
	total_emiss=0.6175 	# (default)		 effective emissivity
	current_emissivity = 0.0017
	current_co2 = 400
	total_ppa = 100000
	co2 = [300,600]
	aa = [1.69963e-03, 1.6988e-03]
	z_air = np.polyfit(co2, aa, 1)
	pa = np.poly1d(z_air)
	co2_c = 1.4967 * year - 2623.7025 # calculate CO2 concentration based on the trend calculated by Lionel  # 0.18
	co2_c_low = (1.4967) * 2030 - 2623.7025 - 0.18 # calculate CO2 concentration based on the trend calculated by Lionel  # 0.18
	co2_c_high = (1.4967) * 2030 - 2623.7025 + 0.18 # calculate CO2 concentration based on the trend calculated by Lionel  # 0.18
	co2_emiss_predict = pa(co2_c)
	other_emissivity = (total_ppa*total_emiss - current_co2*current_emissivity) / (total_ppa - current_co2) # calculate combined emissivity of the other gases
	print(other_emissivity,co2_emiss_predict,co2_c)
	emiss_predict = (co2_emiss_predict*co2_c + other_emissivity*(total_ppa-co2_c)) / total_ppa
	emiss_predict_high = (co2_emiss_predict*co2_c_high + other_emissivity*(total_ppa-co2_c_high)) / total_ppa
	emiss_predict_low = (co2_emiss_predict*co2_c_low + other_emissivity*(total_ppa-co2_c_low)) / total_ppa
	return emiss_predict, emiss_predict_high, emiss_predict_low

def ice_over_time(year):
	# simple trend
	#1979 - 28.3 mln km2
	#2017 - 26.8 mln km2
	ice = -0.03120038 * year + 90.03236626 # ice extent prediction calculated by Pierre-Loup  # +-1.86
	ice_high = (-0.03120038) * 2030 + 90.03236626 + 1.86 # ice extent prediction calculated by Pierre-Loup  # +-1.86
	ice_low = (-0.03120038) * 2030 + 90.03236626 - 1.86 # ice extent prediction calculated by Pierre-Loup  # +-1.86
	ice_2017 = -0.03120038 * 2017 + 90.03236626
	# calculate change of emissivity based on change of ice extent
	es = 510.1 # earth surface km2
	es_al = 0.3 #  earth albedo on 2017
	ice_al = 0.6 # ice albedo
	polar_water_al = 0.35 # water albedo in polar region (corresponding to 80 degreees parallel)
	x = (es*es_al - ice_2017*ice_al) / (es - ice_2017) # calculate albedo of a part of Earth which is not covered by ice
	print(x, ice)
	new_albedo = (ice * ice_al + (es - ice_2017)*x + (ice_2017-ice)*polar_water_al) / es
	new_albedo_high = (ice_high * ice_al + (es - ice_2017)*x + (ice_2017-ice_high)*polar_water_al) / es
	new_albedo_low = (ice_low * ice_al + (es - ice_2017)*x + (ice_2017-ice_low)*polar_water_al) / es
	return new_albedo, new_albedo_high, new_albedo_low

def climate_model(year):
	if year<1979:
		alpha=0.3 		# (default) earth surface albedo
		emiss=0.6175 	# (default)		 effective emissivity
	elif year > 1978 and year < 2100 :
		emiss, emiss_high, emiss_low = co2_over_time(year)
		alpha, alpha_high, alpha_low = ice_over_time(year)
		print(emiss, alpha)
	else:
		print("the future is ambigious, no predictions can be made for this year")		
	S = 1367 		# solar constant
	sigma = 5.67 * 1 / (10**8)			# stephan-boltsman constant
	T = ( ((1-alpha)*S) / (4*emiss*sigma) )**0.25
	T_low = ( ((1-alpha_high)*S) / (4*emiss_high*sigma) )**0.25
	T_high = ( ((1-alpha_low)*S) / (4*emiss_low*sigma) )**0.25
	return T-273.15, T_low - 273.15, T_high - 273.15

model_pred = np.zeros((14))
for i in range(2017,2031):
	if i != 2030:
		model_pred[i-2017] = climate_model(i)[0]
	else:
		model_pred[i-2017], low, high = climate_model(i)

ax.plot(np.array([2017,2031]),np.array([y_old[-1],high]),c='red',linestyle='-.',label='high case scenario (+' + '\u03C3'+')')
ax.plot(np.arange(2017,2031),model_pred,c='green',linestyle='-.',label='0D model prediction')
ax.plot(np.array([2017,2031]),np.array([y_old[-1],low]),c='m',linestyle='-.',label='low case scenario  (-' + '\u03C3'+')')  # plt.xlabel(u'\u03bc = 50')

ax.text(2005,min(air_year),"y=%.2fx2+%.2fx+%.2f" %(f.coef[0],f.coef[1],f.coef[2]),color='blue', fontsize=7)
plt.legend(loc=2,prop={'size': 8})
fig.savefig("/CECI/home/ulg/mast/eivanov/Big_Data/Results/Statistical_approximation_and_prediction.png", dpi=200, bbox_inches='tight')
