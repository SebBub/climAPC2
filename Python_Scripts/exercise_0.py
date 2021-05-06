#%% import of modules needed
import xarray as xr #needed for .nc file data exploration and plotting
from pathlib import Path #really nice library for working with file paths in Python
import matplotlib.gridspec as gridspec #needed to combine cartopy and lineplot, see here https://stackoverflow.com/questions/61433850/combination-of-normal-and-cartopy-subplots-within-the-same-figure
import matplotlib.pyplot as plt # really important module for plotting, there are two main
# ways of how to plot with matplotlib, one is very similar to Matlab, you basically do everything with
# 'plt.*'. The other is object-oriented and allows more flexibility... see here for a nice overview:
# https://matplotlib.org/matplotblog/posts/pyplot-vs-object-oriented-interface/
import cartopy.crs as ccrs #this is for making nice map plots in python, it works really well with xarray dataarrays
import os #this is a Python module for system file operations

# this is a concise way of updating the matplotlib configuration file, which sets fonts, etc.
# see here: https://matplotlib.org/stable/tutorials/introductory/customizing.html
parameters = {'axes.titlesize': 11, 'font.family':'Calibri', 'font.size':11, #'axes.titlesize': 'large'
              'axes.titleweight': 'regular', 'legend.fontsize': 11, 'axes.labelsize':11} # 'mathtext.default':''}
plt.rcParams.update(parameters)

#%% Import the data
data_folder = Path(os.getcwd()) / 'Data' #this sets the Pathlib path for data_folder to the folder "Data"
# where the file is currently stored, this is useful for sharing the file compared to specifying an absolute
# path because the location on the computer, disk name, etc. may vary between systems
file_name='sst.mnmean.nc'
sst_ds = xr.open_dataset(data_folder / file_name) #open the dataset and load it into a variable, you can combine different


# pathlib Path objects using the '/'
print('here we actually pass:')
print(data_folder / file_name)
print('which is the file path to the .nc file we want to open... the argument is a Pathlib filepath object, which almost all '
      'modules like Xarray, Pandas, etc. understand. \n'
      'I guess you already see why Pathlib is so nice. :-) \n')

# because the file sst.mnmean.nc has multipe data variables, we have opened it as an xarray dataset,
# not a dataarray ... we can look at all the variables in the dataset with:
print(sst_ds.data_vars)

# now we read the specific data variable that we are interested in (in this case SST) into the memory
# as a dataarray, built from the dataset, it's as simple as...
sst = sst_ds.sst
# note that writing the variables lat, lon, time
print(sst.coords)
# into separate variables as was done in the Matlab script is not necessary... also coversions of time
# are not necessary because Xarray automatically reads in the time dimension in Python datetime format:
# check the previous print command and look for '* time', there you can see the datatype

# now we proceed to masking the SST data over land, as there SST does not make sense
mask = xr.open_dataarray(data_folder / 'lsmask_sst.nc')
print(type(mask.values)) #here you can see that data inside Xarray dataarrays is actually stored in
# nested Numpy ndarrays... let's have a look at the data, we can see that the mask contains a grid
# of 0 and 1, relating to gridpoints over land (0) and gridpoints over sea (1)
print(mask.values)

#let us check if our mask and the SST dataset have the same dimensions and gridsize... only if this
# is the case, we can use the mask for the dataarray... otherwise we need to regrid the mask to the
# same dimensions and gridspacing as the dataarray

print(mask.lat.shape)
print(sst.lat.shape)

# OK, both datasets have 180 gridpoints for latitude... to be sure, let us check:
print(mask.lat == sst.lat) #checks element-wise, if the data is identical
print(mask.lon == sst.lon) # also looks good...

# Now we can easily use the .where - operation, see http://xarray.pydata.org/en/stable/generated/xarray.DataArray.where.html#xarray.DataArray.where
# this also returns a dataarray, but onyl the values / the data at the gridpoints
# where the conditions is True

sst_masked = sst.where(mask.values == 1) #write masked dataarray to new variable

# I think we can also reorder the dimensions of the main file as in Matlab
# but I do not think this is necessary in Python ...

#%% Plotting
# plot the first timestep of the dataset to check everything looks good
fig1, ax1 = plt.subplots() # create a figure and an axes instance with the function subplots() and
# assign them to two variables for the way forward
sst_masked[0,:,:].plot.contourf(ax=ax1) # see http://xarray.pydata.org/en/stable/indexing.html for indexing and selecting data
fig1.show() #you can see that Xarray with its builtin plotting functions (called in the line above) does a really nice job
# in labelling based on the data contained in the dataarray ... everything can be customized of course

# now plot a timeseries at a selected point to also check the dataset is OK w.r.t. time:
fig2, ax2 = plt.subplots()
sst_masked[:,80,200].plot(ax=ax2)
fig2.show() #looks good!

#%%
# now let us make the plots pretty... consider that the global changes for fonts, font sizes etc. were all made already
# at the beginning of the script by customizing matplotlib rcparams
fig3 = plt.figure()
fig3.subplots_adjust(hspace=.8)
gs = gridspec.GridSpec(3,3,fig3) #see here https://matplotlib.org/stable/api/_as_gen/matplotlib.gridspec.GridSpec.html#matplotlib.gridspec.GridSpec
ax3 = fig3.add_subplot(gs[0:2,:], projection = ccrs.PlateCarree()) # the last is needed to
# specify which map projection should be used by cartopy for the plot, here we use plate carree, see
# https://scitools.org.uk/cartopy/docs/latest/crs/projections.html also we can now call all cartopy functions etc.

ax3.coastlines(resolution='auto', color='k')
ax3.gridlines(color='lightgrey', linestyle='-', draw_labels=True)


ax4 = fig3.add_subplot(gs[2, :])

# 1st row plotting
sst_masked[0,:,:].plot.contourf(ax=ax3)
ax3.scatter(x=200,y=80, s=70, c='black', marker='x') #it is important to draw the point after the contourf plot because
# otherwise the contourf would hide the point

# 2nd row plotting
sst_masked[:,80,200].plot(ax=ax4)

#Labels, etc.
ax3.set(xlabel='Longitude', ylabel='Latitude', title=('SST during ' + sst_masked.time.values[0].astype('str')[0:7]))
ax4.set(xlabel='Time', ylabel='SST in degC', title='SST timeseries')

fig3.show()
# fig3.savefig('Nice_figure.png', dpi=300) #save to .png file

###############################################
# now follow the examples from exercise_0.m

#%%
fig4, ax5 = plt.subplots()

# profile of SST along near the equator at the first time step
sst_masked[0,:,:].where(sst_masked.lat == .5, drop=True).plot(ax=ax5, c='r', label=sst_masked.time.values[0].astype('str')[0:7]) # see here about the selecting http://xarray.pydata.org/en/stable/indexing.html
ax5.grid()
ax5.set(title='Equatorial SST during ' + sst_masked.time.values[0].astype('str')[0:7])

# Let's overlay the mean profile averaged throughout the entire data set
sst_masked.where(sst_masked.lat == .5, drop=True).mean(dim='time').plot(ax=ax5, c='k', label='Time mean')
ax5.legend(loc='lower left')
fig4.show()

#%% Hovmoller diagram
fig5, ax6 = plt.subplots()
sst_masked[1:120,:,:].mean(dim='lon').plot.contourf(ax=ax6)
fig5.show()