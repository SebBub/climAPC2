# Example script to find locations and times of minimum and maximum values

#Import
from pathlib import Path
import os
import xarray as xr
import matplotlib.pyplot as plt

#%% Import the data
data_folder = Path(os.getcwd()) / 'Data'
file_name='sst.mnmean.nc'
sst = xr.open_dataset(data_folder / file_name).sst
mask = xr.open_dataarray(data_folder / 'lsmask_sst.nc')
sst_masked = sst.where(mask.values == 1)

#now apply a running mean filter to smooth the dataset
sst_rm = sst_masked.rolling(time=3).mean()
# Identify the maximum at each location over time
max_sst = sst_rm.max(dim = 'time')
# Find the maximum over all longitudes as a function of latitude:
max_sst_lat = max_sst.max(dim = 'lon')
# Get the latitude of the maximum temperature:
max_lat = max_sst_lat.idxmax('lat').values.item() # see https://stackoverflow.com/questions/40179593/how-to-get-the-coordinates-of-the-maximum-in-xarray
# note that .item converts 1D ndarray into scalar, see here: https://numpy.org/doc/stable/reference/generated/numpy.ndarray.item.html

# now find the maximum longitude
max_lon = max_sst.max(dim = 'lat').idxmax('lon').values.item()

#now find the maximum over time at the identified location
#TODO
max_sst_time = sst_rm[:,max_lat,max_lon].max(dim='time').values.item()

# See what country or ocean is at that location:

#%% Some random plots
# Maximum SST as a function of latitude
fig1, ax1 = plt.subplots(nrows=2)
fig1.subplots_adjust(hspace=.3)
max_sst_lat.plot(ax=ax1[0])
max_sst.max(dim = 'lat').plot(ax=ax1[1])
fig1.show()

