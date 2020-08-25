from netCDF4 import Dataset
import numpy as np
from decimal import Decimal
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.basemap import Basemap


data = Dataset('jb_iby_sry_gtn_giy.nc', 'r', format='NETCDF4')

seaLevelVariable = data.variables['height']

# Extract longitude and latitude
lons = data.variables['lon'][:]
lats = data.variables['lat'][:]
# Extract sea level in millimeters
seaLevel = seaLevelVariable[:,:,:].data
# Extract time periods
months = data.variables['month'][:]
years = data.variables['year'][:]


minValue = 99999
maxValue = -99999
seaLevelDFList = []
for i in range(311):
	level = seaLevel[i]
	levelDF = pd.DataFrame(level)

	# Purge values == 32700 (land)
	levelDF = levelDF[levelDF != 32700]

	# Compute min and max sea level values (for the colorbar)
	if min(levelDF.min()) < minValue:
		minValue = min(levelDF.min())

	if max(levelDF.max()) > maxValue:
		maxValue = max(levelDF.max())

	seaLevelDFList.append(levelDF)


sea_units = seaLevelVariable.units

# Get some parameters for the Stereographic Projection
lon_0 = lons.mean()
lat_0 = lats.mean()

m = Basemap(resolution='l',projection='mill',\
	            lat_0=lat_0,lon_0=lon_0)

# Because our lon and lat variables are 1D,
# use meshgrid to create 2D arrays
lon, lat = np.meshgrid(lons, lats)
xi, yi = m(lon, lat)

"""
for i in range(311):

	# Plot Data
	cs = m.pcolor(xi,yi,np.squeeze(seaLevelDFList[i].values), vmin=minValue , vmax=maxValue , cmap="jet")


	# Add Coastlines, States, and Country Boundaries
	m.drawcoastlines()
	m.drawstates()
	m.drawcountries()

	# Add Colorbar
	cbar = m.colorbar(cs, location='bottom', pad="10%")
	cbar.set_label(sea_units)

	plt.text(0.4, 1.2, str(years[i]) + "-" + str(months[i]))
	# Add Title
	plt.title('Sea level w.r.t reference point')
	#plt.savefig("test{}.jpg".format(i))
	plt.close()

"""

meanSeaLevels = np.zeros((311, 1))
for i in range(311):
	meanSeaLevels[i] = np.nanmean(seaLevelDFList[i].values)

print(meanSeaLevels)

x = np.linspace(1, 311, 311)
plt.plot(x, meanSeaLevels)
plt.xlabel("Time period")
plt.ylabel("Mean sea level (w.r.t. reference level)")
plt.title("Evolution of the mean sea level")
plt.savefig("evol_meansealevel.png")
plt.show()

# Mean variation between two months
meanSeaLevelDF = pd.DataFrame(meanSeaLevels)
print("Mean change in mean sea level (monthly):", np.nanmean(meanSeaLevelDF.diff()))