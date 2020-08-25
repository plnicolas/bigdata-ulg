import numpy as np

def get_area_from_geo_coord(lat_up,lat_down,lon_left,lon_right):

	lat_center_deg = lat_down + (lat_up - lat_down)/2
	lon_center_deg = lon_left + (lon_right - lon_left)

	lat_center_km = (lon_right-lon_left)/360 * 40075 * np.cos(np.radians(lat_center_deg))
	lon_center_km = (lat_up - lat_down) * 111

	return lat_center_km*lon_center_km
