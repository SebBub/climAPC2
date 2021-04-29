% Using MATLAB to work with gridded climate data in NetCDF format

% NetCDF stands for Network Common Data Form. NetCDF is a ubiquitous file
% format for storing gridded data (observational data sets, model output,
% etc.) in the Earth sciences, especially in the realm of oceanography,
% atmospheric science, and climate science. Gridded means the data are
% arranged spatially on a 2D (e.g., lat x lon) or 3D (e.g., lat x lon x
% depth) grid, and can have a time dimension as well. NetCDF files commonly
% use the file extension .nc or .cdf (e.g., hadslp2.nc).

% There are at least three different ways to open NetCDF data in MATLAB.
%
% (1) Install and use a third-party toolbox such as mexnc/mexcdf.
% (2) Use MATLAB's built-in low-level NetCDF functions, which have the
%     convention netcdf.xxxx (e.g., netcdf.open, netcdf.getAtt, etc.)
% (3) Use MATLAB's built-in high-level NetCDF functions, which have the
%     convention ncxxxx (e.g., ncread, ncreadatt, etc.)
%
% Method (1) can be very useful but requires a fair amount of expertise in
% a UNIX-like environment to install and relies on other technical software
% libraries being installed first. Method (2) allows much control over the
% importing process and writing your own NetCDF files. Method (3) is
% probably the simplest and requires nothing more than MATLAB as installed
% "out of the box."
%
% In case you are interested:
%
% Variables in NetCDF files are often stored as nonsense numbers that need
% to be adjusted in order to take physically meaningful values. For
% example, an actual value of SST at some place and time might be 26.42
% (deg C). In the NetCDF file, this might have been stored as -861858. To
% convert all of the data in the NetCDF file to physically meaningful SST
% values, one must apply the scale_factor (an "attribute" of the variable,
% 0.01 in this case) and the add_offset (another "attribute" of the
% variable, 8645 in this case). In addition, bad or missing data might be
% flagged with a constant like -9999. Those values must be converted to NaN
% ("not a number") so MATLAB knows not to plot it or let it affect any
% calculations. Method (3) takes care of all of this automatically.

% A note about the order of dimensions in a variable:
%
% Take for example a variable that has just been imported into MATLAB, and
% the dimensions have been pre-arranged like this: SST(lon,lat,time). Many
% operations within MATLAB (e.g., plotting data) and most of the scripts I
% have written work with the least amount of fuss if the dimensions are
% ordered like this: SST(time,lat,lon). This may seem strange at first, but
% life will be simplified if one always makes a point to have all variables
% stored in the MATLAB workspace with the dimensional arrangement f(t,y,x)
% or, if there is a vertical dimension, f(t,z,y,x).
%
% Fortunately, rearranging the order of the dimensions is very simple and
% achieved with a single call of the permute function. The syntax is:
% permute(variable,[order]);
% where order would be 3 2 1 to rearrange SST(lon,lat,time) into
% SST(time,lat,lon).

% A note about dealing with time:
%
% Obviously, time is important in studying climate variability. That a data
% set contains monthly mean values of some variable spanning. say, Jan.
% 1854 through Dec. 2011 is not arbitrary. Dealing with the concept of time
% in MATLAB can be frustrating. Time in MATLAB is represented as fractional
% days since Jan. 1, 0000. For example, at the time of writing this, MATLAB
% says the time this very instant is 734895.4970484802, because that is how
% many days have passed since the instant the the Earth began its first
% rotation on Jan. 1, 0000. If you want to know what the wind field looked
% like during Dec. 1982, how would you know that Dec. 1982 is the 1548th
% time step in the data set? Or that the "MATLAB time" for Dec. 1, 1982 is
% 724246? To make matters worse, most of the NetCDF climate data sets I
% have come across express time in units of days since Jan. 1, 1800 (not
% 0000). This probably has to do with the instrumental record generally
% beginning in the mid-1800s. So, you have to manually confirm that, add on
% the number of days between Jan. 1, 0000 and Jan. 1, 1800, and THEN use
% one of MATLAB's functions to convert that fractional year into something
% readable (to us) like "Dec 1982." Fortunately, this can all be achieved
% with a few simple lines of code which are included in the example below.


% First, let's clean up the workspace. Be careful, as this will wipe out
% all variables currently in the workspace and close all figure windows.

close all, clear all

% Now, let's check out the contents of the NetCDF file that contains some
% data we are interested in working with. To do that, we will use the the
% built-in high-level NetCDF functions ncdisp.

% Enter the name of the file you would like to read and make sure it is in
% the Matlab seach path. To do so, you can go to the corresponding folder
% in the Matlab editor, righ-click on it and select "Add to path ->
% Selected Folder (and Sub-Folders)"

file_name='sst.mnmean.nc';

ncdisp(file_name,'/','min')

% So, the file contains lat, lon, time, time_bnds (ignore this), and sst.
% sst is a function of lon,lat, and time. The number of lats, lons, and
% time steps are also indicated.

% We know this is a monthly data set, but let's just confirm that the time
% variable in the file is expressed as days since Jan. 1, 1800. To do that,
% we'll use another high-level NetCDF function ncreadatt to query the units
% attribute of the time variable. Display it to the screen by wrapping the
% MATLAB function disp around it.

disp( ncreadatt(file_name,'time','units') )

% Now that we know what's what, let's read in the main variable and the
% dimensions that it is a function of. This is easy with ncread.

sst=ncread(file_name,'sst');
lon_sst=ncread(file_name,'lon');
lat_sst=ncread(file_name,'lat');
time_sst=ncread(file_name,'time');

% Note: We could have just as well called our dimension variables lon, lat,
% time, or anything whatsoever. However, it is good practice to give a
% fairly specific name to the dimensions. This is in case we eventually
% want to open another data set that is on a different spatial grid or
% spans a different interval of time. When we do that, we'll have multiple
% lat, lon, and time dimensions - one approprite to the dimensions of each
% variable.

% Let's convert the time we just imported to the time that MATLAB
% understands, and make an extra version of time that WE can understand.

time_sst=time_sst+datenum(1800,1,1);    % Add all of the days between 0000 and 1800
time2_sst=datestr(time_sst,'mmm yyyy'); % Make a version we can understand

% To be really sure, let's just check what is the first and last time.

disp( time2_sst(1,:) )
disp( time2_sst(end,:) )


% Since SST Data over land does not make sense, let's remove it. To do so,
% we first make sure that the mask the same dimension as the SST. There are
% many more efficient ways to remove the data over land, but this one is
% quick for now...

mask=ncread('lsmask_sst.nc','mask');
mask3=repmat(mask,[1 1 size(sst,3)]);
sst(mask3==0)=NaN; clear mask3



% When we imported the sst variable, it was originally a function of lon,
% lat and time, in that order. Let's reorder the dimensions of the main
% variable to be time x lat x lon. We want what is currently the third
% dimension (time) to come first, then what is currently the second
% dimension (lat), then what is currently the first dimension (lon)... or
% 3, 2, then 2.

sst=permute(sst,[3 2 1]);

% All set! Now we have everything we need in our workspace, and nothing
% more. We are ready to begin working with the data. Let's make a quick
% lat-lon plot of our main variable at the first time step, just to be sure
% it looks right.

figure(1)
contourf(lon_sst,lat_sst,squeeze(sst(10,:,:)),30,'linestyle','none');
colorbar

% Yep, that looks like an SST field, the continents appear to be in the
% right place (and do not have SST values, of course!), and the numbers
% along the color scale seem reasonable for SST in degrees C.

% Now let's make sure everything also looks OK w.r.t. time. We can do this
% by plotting a time series of SST at some random location in the central
% Pacific.


figure(2)
plot(sst(:,80,200))



% Yep, looks like a time series of SST somewhere in the tropics. Zoom in on
% the time series and see that the values oscillate strongly with a period
% of 12 time steps. This makes sense since we are working with a monthly
% data set and SST has an annual cycle in most regions of the world.
% Zooming back out of course we see there is more going on than just the
% annual cycle, such as interannual, decadal and multidecadal variability,
% as well as perhaps a trend.

% Our plots aren't very pretty. Let's make a publication-quality figure
% containing both plots, including a marker on the map indicating where our
% time series is from, and a readable time axis on the time series.



close 1 2

figure(1)
subplot(1,2,1)
    [~,h]=contourf(lon_sst,lat_sst,squeeze(sst(1,:,:)),50);
    set(h,'EdgeColor','none')
    colorbar
    caxis([0 28])
    set(gca,'FontSize',14)
    hold on
    contourf(lon_sst,lat_sst,-squeeze(mask)',[0 0],'facecolor',[0.7 0.7 0.7],'color','k','linewidth',2)
    %contour(lon_sst,lat_sst,-squeeze(mask)',[0 0],'color','k','linewidth',2)
    plot(lon_sst(200),lat_sst(80),'ro','MarkerFacecolor','r')
    hold off
    xlabel('Longitude (\circE)')
    ylabel('Latitude (\circN)')
    title(['SST during ',time2_sst(1,:)])
subplot(1,2,2)
    plot(time_sst,sst(:,80,200))
    set(gca,'FontSize',14)
    datetick('x','yyyy')
    box on
    grid on
    xlabel('Year')
    ylabel('SST (\circC)')
    title(['SST at ',num2str(lat_sst(80)),'\circN, ',num2str(lon_sst(200)),'\circE'])


% In order to make nice map plots, I recommend downloading the Mapping
% Toolobox M_map (https://www.eoas.ubc.ca/~rich/mapug.html). It is well
% documented and lets you choose from many different map projections,
% applies proper gridding to the latitude, longititude coordinates and
% draws pretty coastlines. Even though installing the toolbox is not
% required for this class, it would save you much trouble later ;-)
    
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% Examples


% In MATLAB, you have the flexibility to analyze anything w.r.t. any
% dimension or combination of dimensions you can think of. Let's make a
% profile of SST along near the equator at the first time step.

figure(2)
plot(lon_sst,squeeze(sst(1,lat_sst==0.5,:)),'color','k','linewidth',2)
box on
grid on
set(gca,'FontSize',14)
xlabel('Longitude (\circE)')
ylabel('SST (\circC)')
title(['Equatorial SST during ',time2_sst(1,:)])

% Let's overlay the mean profile averaged throughout the entire data set

figure(2)
hold on
plot(lon_sst,squeeze(mean(sst(:,lat_sst==0.5,:),1)),'color','r','linewidth',2)
hold off
legend(['SST during ',time2_sst(1,:)],'Mean over the full period')
title('Equatorial SST profiles')

% Rather than looking at equatorial SST profiles at individual times or
% static profiles averaged over an interval of time, let's make another
% figure where we bring in the time dimension, i.e. longitude still along
% the x-axis, but now time along the y-axis and the values will be
% represented by colors. This is called a Hovmoller diagram. Let's restrict
% the time dimension to the first 10 years (120 time steps).

figure(3)
[~,h]=contourf(lon_sst,time_sst(1:120),squeeze(sst(1:120,lat_sst==0.5,:)),50);
set(h,'EdgeColor','none')
colorbar
caxis([20 30])
set(gca,'FontSize',14)
datetick('y','yyyy')
xlabel('Longitude (\circE)')
ylabel('Year')
title(['Equatorial SST Hovmoller'])
ylim([time_sst(1) time_sst(120)])



% If you would like to go through another example, you can run
% 'example_quiver', which also shows how to include velocity vectors for
% wind data using the MATLAB function 'quiver'.
