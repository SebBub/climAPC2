#%%
import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import cmocean as cmo
import os
from pathlib import Path
import cartopy.crs as ccrs

parameters = {'axes.titlesize': 11, 'font.family':'Calibri', 'font.size':11, #'axes.titlesize': 'large'
              'axes.titleweight': 'regular', 'legend.fontsize': 11, 'axes.labelsize':11} # 'mathtext.default':''}
plt.rcParams.update(parameters)

#%% read data
data_folder = Path(os.getcwd()) / 'Data'
file_name = 'air.mon.mean.nc'
airtemp_da_rm = xr.open_dataset(data_folder / file_name).air.rolling(time=3).mean()
# append '.air' to directly read the data array with the variable air, not the dataset
# directly apply rolling filter to get monthly timeseries of three-months rolling average

#%% Part 1
# Which location was warmest/coldest/hottest/driest within a three-month period? What
# country or ocean is at that location? At what time did the maximum (minimum) occur?
t_max = airtemp_da_rm.max(dim='time')
t_min = airtemp_da_rm.min(dim='time')

#TODO build this code into a reusable function that can be applied to any (time, lat, lon) dataarray, code in Helper_functions.py
t_max_lat = t_max.max(dim='lon').idxmax('lat').values.item()
t_max_lon = t_max.max(dim='lat').idxmax('lon').values.item()
t_max_time = airtemp_da_rm.sel(lat=t_max_lat, lon=t_max_lon).idxmax('time').values #http://xarray.pydata.org/en/stable/generated/xarray.DataArray.sel.html#xarray.DataArray.sel
print('The maximum occured at lat ' + str(t_max_lat) + ' and lon ' + str(t_max_lon) + ' at on ' + str(t_max_time)[0:10])

t_min_lat = t_min.min(dim='lon').idxmin('lat').values.item()
t_min_lon = t_min.min(dim='lat').idxmin('lon').values.item()
t_min_time = airtemp_da_rm.sel(lat=t_min_lat, lon=t_min_lon).idxmin('time').values #http://xarray.pydata.org/en/stable/generated/xarray.DataArray.sel.html#xarray.DataArray.sel
print('The minimum occured at lat ' + str(t_min_lat) + ' and lon ' + str(t_min_lon) + ' at on ' + str(t_min_time)[0:10])


# 1. Create a map plot of the minimum and maximum temperature (precipitation):
# we also mark the location of the max and min with a cross so that we can locate the country

fig1, ax1 = plt.subplots(nrows=2, subplot_kw = {'projection':ccrs.PlateCarree()})
fig1.subplots_adjust(hspace=.4)

t_max.plot.contourf(ax=ax1[0], vmin=-75, vmax=50, cmap='coolwarm', cbar_kwargs = {'label':'Air temperature'})
ax1[0].coastlines(resolution='110m')
ax1[0].set(title = 'Maximum air temperature')
gl1 = ax1[0].gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='gray', alpha=0.7, linestyle='--')
gl1.top_labels = False
gl1.right_labels = False
ax1[0].scatter(x=t_max_lon, y=t_max_lat, s=70, c='magenta', marker='x')

t_min.plot.contourf(ax=ax1[1], vmin=-75, vmax=50, cmap='coolwarm', cbar_kwargs = {'label':'Air temperature'})
ax1[1].coastlines(resolution='110m')
ax1[1].set(title = 'Minimum air temperature')
gl2 = ax1[1].gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='gray', alpha=0.7, linestyle='--')
gl2.top_labels = False
gl2.right_labels = False
ax1[1].scatter(x=t_min_lon, y=t_min_lat, s=70, c='magenta', marker='x')
fig1.tight_layout()
fig1.show()
fig1.savefig('min_max_locations.png', dpi=400)

#%% Part 2 #TODO
# Compute the monthly climatologies and anomalies with 'climanom.m' and the range of the
# mean seasonal cycle.
# Plot the mean temperature (precipitation), averaged over all longitudes, as a function of
# latitude. Do the same for the range of the monthly climatologies.
# Now, subtract the latitudinal means and create a map plot of the result. Can you explain the
# results?