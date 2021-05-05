#TODO filter for North America extent

#%%
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import cmocean as cmo
import scipy
import os
from pathlib import Path
import cartopy.crs as ccrs
from statsmodels.tsa.seasonal import seasonal_decompose

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
t_max = airtemp_da_rm.max(dim='time') #reduce time dimension and select maximum value over all timsteps at each gridpoint
t_min = airtemp_da_rm.min(dim='time')

t_max_lat = t_max.max(dim='lon').idxmax('lat').values.item() #reduce the longitude dimension and select the maximum for
# each gridpoint, then determine latitude value of the maximum along latitude with idxmax() function, extract corres-
# ponding Numpy array by appending .values and convert 1D array to float by appending .item()
t_max_lon = t_max.max(dim='lat').idxmax('lon').values.item() #same but for longitude
t_max_time = airtemp_da_rm.sel(lat=t_max_lat, lon=t_max_lon).idxmax('time').values # see http://xarray.pydata.org/en/stable/generated/xarray.DataArray.sel.html#xarray.DataArray.sel
print('The maximum occured at lat ' + str(t_max_lat) + ' and lon ' + str(t_max_lon) + ' at on ' + str(t_max_time)[0:10])

#same for minimum
t_min_lat = t_min.min(dim='lon').idxmin('lat').values.item()
t_min_lon = t_min.min(dim='lat').idxmin('lon').values.item()
t_min_time = airtemp_da_rm.sel(lat=t_min_lat, lon=t_min_lon).idxmin('time').values #http://xarray.pydata.org/en/stable/generated/xarray.DataArray.sel.html#xarray.DataArray.sel
print('The minimum occured at lat ' + str(t_min_lat) + ' and lon ' + str(t_min_lon) + ' at on ' + str(t_min_time)[0:10])


# 1. Create a map plot of the minimum and maximum temperature (precipitation):
# we also mark the location of the max and min with a cross so that we can locate the country

fig1, ax1 = plt.subplots(nrows=2, subplot_kw = {'projection':ccrs.PlateCarree()})
fig1.subplots_adjust(hspace=.4)

t_max.plot.contourf(ax=ax1[0], vmin=-75, vmax=50, cmap='coolwarm', cbar_kwargs = {'label':'Air temperature in °C'})
ax1[0].coastlines(resolution='110m')
ax1[0].set(title = 'Maximum air temperature')
gl1 = ax1[0].gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='gray', alpha=0.7, linestyle='--')
gl1.top_labels = False
gl1.right_labels = False
ax1[0].scatter(x=t_max_lon, y=t_max_lat, s=70, c='magenta', marker='x')

t_min.plot.contourf(ax=ax1[1], vmin=-75, vmax=50, cmap='coolwarm', cbar_kwargs = {'label':'Air temperature in °C'})
ax1[1].coastlines(resolution='110m')
ax1[1].set(title = 'Minimum air temperature')
gl2 = ax1[1].gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='gray', alpha=0.7, linestyle='--')
gl2.top_labels = False
gl2.right_labels = False
ax1[1].scatter(x=t_min_lon, y=t_min_lat, s=70, c='magenta', marker='x')
fig1.tight_layout()
fig1.show()
fig1.savefig('min_max_locations.png', dpi=400)

#%% Part 2
# Compute the monthly climatologies and anomalies with 'climanom.m' and the range of the mean seasonal cycle.

# see here: http://xarray.pydata.org/en/stable/time-series.html - "Resampling and grouped operations"
# monthly climatology:
airtemp_da_rm_climM = airtemp_da_rm.groupby('time.month').mean()
# monthly anomaly:
airtemp_da_rm_anomM = airtemp_da_rm.groupby('time.month') - airtemp_da_rm_climM #TODO not entirely sure this is the correct way.. needs some plausibility check of the results

# Plot the mean temperature (precipitation), averaged over all longitudes, as a function of latitude.
airtemp_da_rm_lat = airtemp_da_rm.mean(dim='lon').mean(dim='time')

fig2, ax2 = plt.subplots()
airtemp_da_rm_lat.plot(ax=ax2, label='Time and zonal mean temperature')
airtemp_da_rm.max(dim='time').mean(dim='lon').plot(ax=ax2, label='Maximum values')
airtemp_da_rm.min(dim='time').mean(dim='lon').plot(ax=ax2, label='Minimum values')

ax2.set(xlabel='Latitude in degrees North', ylabel='Air temperature in °C', title='Zonal and time mean temperature as a function of latitude')
ax2.grid()
ax2.legend(loc='lower right')
fig2.show()
fig2.savefig('mean_t_lat.png', dpi=500)

# Do the same for the range of the monthly climatologies.
airtemp_da_rm_ClimRange = (airtemp_da_rm_climM.max('month') - airtemp_da_rm_climM.min('month')).mean(dim='lon')

fig3, ax3 = plt.subplots()
airtemp_da_rm_ClimRange.plot(ax=ax3, label='Zonal mean of range ')
# airtemp_da_rm_lat.plot(ax=ax3, label='Time and zonal mean temperature')
ax3.set(xlabel='Latitude in degrees North', ylabel='Air temperature in °C', title='Absolute maximum range of zonal averaged temperature ranges')
ax3.grid()
ax3.legend(loc='lower right')
fig3.show()
fig3.savefig('mean_range_lat.png', dpi=500)

# Now, subtract the latitudinal means and create a map plot of the result. Can you explain the results?
fig4, ax4 = plt.subplots(subplot_kw = {'projection':ccrs.PlateCarree()})

#TODO wie kann man nun bei dem Term der abgezogen wird das wieder um eine longitude erweitern?! Damit wir es
# von unserer 2D Karte abziehen können...? Wir erwarten z.B. dass die Kontinente wärmer sind als die Ozeane...
(airtemp_da_rm.mean(dim='time') - airtemp_da_rm.mean(dim='time').mean(dim='lon')).plot(ax=ax4)
ax4.coastlines()
gl = ax4.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='gray', alpha=0.7, linestyle='--')
gl.top_labels = False
gl.right_labels = False
fig4.show()


#%% Part 3
# Pick a location near Kiel. Alternatively, you can use a 2x2 box around Kiel and compute the average or use the MATLAB
# function 'interp2.m' to interpolate to Kiel. For the purpose of this exercise, all methods are fine but be aware that
# generally, interpolation in lat-lon coordinates should be weighted to take their true distances into account.

# we start by visualising the grid of the .nc file to better understand what interpolation/regridding is necessary, if any
fig5, ax5 = plt.subplots(subplot_kw={'projection' : ccrs.PlateCarree()})
gl = ax5.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='k', alpha=1, linestyle='--')
gl.xlocator = mticker.FixedLocator(np.round(airtemp_da_rm.lon.values, 2))
gl.ylocator = mticker.FixedLocator(np.round(airtemp_da_rm.lat.values, 2))
ax5.coastlines(resolution='10m')
ax5.scatter(10.1228, 54.3233, s=50, marker='x', color='r')
ax5.set(ylabel='Latitude in deg North', xlabel='Longitude in deg East', title='Grid of source file')
ax5.set_extent([5,15,52,57])
ax5.annotate('Kiel', (.53,.44), xycoords='axes fraction', weight='bold', c='r')
fig5.show()

# for simplification we chose the gridpoint at 55°N 10°E, no regridding or interpolation, can be added later
# Plot the temporal evolution of the temperature (precipitation) for this location or box, with the seasonal cycle subtracted
# (i.e. use the result obtained from 'climanom.m'). Was there a significant trend over the time period?

airtemp_Kiel = airtemp_da_rm.sel(lat=55, lon=10)
################################################################################################
# this is just for a brief check... statsmodels has a nice decomposition into trends, etc...
decompose_result = seasonal_decompose(airtemp_Kiel.to_dataframe().drop(labels=['lat','lon'], axis='columns').dropna(),
                                      model="additive")
decompose_result.plot()
plt.show()
################################################################################################

# good tutorial on deseasonalization: https://machinelearningmastery.com/time-series-seasonality-with-python/
# see here: https://xarray.pydata.org/en/v0.5/examples/weather-data.html#calculate-monthly-anomalies !
clim = airtemp_Kiel.groupby('time.month').mean()
anom = airtemp_Kiel.groupby('time.month') - clim

# we now look at a 1D timeseries at point lat,lon, hence we extract the data from the multidimensional data array
# into a pandas dataframe (with a datetime index for easy datetime functionality)
Kiel_anom = anom.to_dataframe().drop(labels=['lat','lon'], axis='columns').dropna() #also drop NaNs

#now we continue with this nice and clean 1D Pandas timeseries and plot for Kiel:
fig5, ax5 = plt.subplots()
Kiel_anom['air'].plot(ax=ax5)
ax5.set(xlabel='Time', ylabel='Temperature anomaly', title='Kiel monthly temperature anomaly')
fig5.show()

#%%
# Create a map of the trends using the function 'sigtrendmap.m'. Was there a significant trend anywhere else?
# If you are in the precipitation group, please use precipitable water. Feel free to follow 'example_3.m' to complete this task.

# http://atedstone.github.io/rate-of-change-maps/

vals = airtemp_da_rm.dropna(dim='time').values #here we also delete any timestep with NaNs included in the array because
# polyfit cannot handle NaNs
months = range(0, len(airtemp_da_rm.dropna(dim='time').time.values)) #because the months are stored as numpy datetime we need some ascending list of integers
# as our "x-coordinate" for performing a linear regression over time
vals2 = vals.reshape(len(months), -1)
regressions = np.polyfit(months, vals2, 1)
trends_month = regressions[0,:].reshape(vals.shape[1], vals.shape[2])
trends_year = trends_month * 12

# now rebuild an Xarray dataarray with the numpy array that contains the trends
trends_array = xr.DataArray(data=trends_year,
                            coords=[('lat',airtemp_da_rm.lat.values), ('lon',airtemp_da_rm.lon.values)],
                            name='Trends in air temperature',
                            attrs=dict(
                                description="Linear air temperature trends",
                                units="degC / year"
                            ))

fig6, ax6 = plt.subplots(subplot_kw={'projection' : ccrs.PlateCarree()})
ax6.coastlines()
trends_array.plot(ax=ax6)
gl = ax6.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=1, color='k', alpha=1, linestyle='--')
gl.right_labels  = gl.top_labels = False
ax6.set(title=trends_array.name)
fig6.show()


# Voluntary: Create maps of the trends for only winter or only summer or investigate trends
# for a specific region as a function of month.