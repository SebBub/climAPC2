% Example script to plot the mean sea level pressure and surface winds

clear all, close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load the data


file_name_1='slp.mon.mean.nc';
file_name_2='uwnd.mon.mean.nc';
file_name_3='vwnd.mon.mean.nc';
file_name_4='land.nc';

slp=ncread(file_name_1,'slp');
u=ncread(file_name_2,'uwnd');
v=ncread(file_name_3,'vwnd');
land=ncread(file_name_4,'land');

% in this example, we all variables have the same coordinates, so we do not
% need to distinguish them...

lon=ncread(file_name_1,'lon');
lat=ncread(file_name_1,'lat');
time=ncread(file_name_1,'time');

% Add the days between 0000 and 1800 and convert hours to days:

time=time/24+datenum(1800,1,1);
time2=datestr(time,'mmm yyyy');


slp=permute(slp,[3 2 1]);
u=permute(u,[3 2 1]);
v=permute(v,[3 2 1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute the means

slp_mean=squeeze(mean(slp,1));
u_mean=squeeze(mean(u,1));
v_mean=squeeze(mean(v,1));

% compute the speed of the mean winds. Note that this is different from the
% mean wind speed!!!

spd_mean=sqrt(u_mean.^2 + v_mean.^2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We can only show a few arrows, so we smooth the velocity field...

win_lon=7; win_lat=3;
u_mean_rm=smooth_map(u_mean,win_lon,win_lat);
v_mean_rm=smooth_map(v_mean,win_lon,win_lat);


% ... and only keep the central vectors:

N_lat=length(lat); N_lon=length(lon);

u_mean_rm=u_mean_rm(1+(win_lat-1)/2:win_lat:N_lat-(win_lat-1)/2,...
    1+(win_lon-1)/2:win_lon:N_lon-(win_lon-1)/2);
v_mean_rm=v_mean_rm(1+(win_lat-1)/2:win_lat:N_lat-(win_lat-1)/2,...
    1+(win_lon-1)/2:win_lon:N_lon-(win_lon-1)/2);

lat_rm=lat(1+(win_lat-1)/2:win_lat:N_lat-(win_lat-1)/2);
lon_rm=lon(1+(win_lon-1)/2:win_lon:N_lon-(win_lon-1)/2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make a plot using quiver

% In map projections, we need to scale the eastward velocity with cosd(lat)
% because because Matlab does not know that the distance between longitudes
% decreases towards the poles. Here we do not use any projection...

[lonn_rm,latt_rm]=meshgrid(lon_rm,lat_rm);


figure(1)
subplot(1,2,1)
    [~,h]=contourf(lon,lat,slp_mean,50);
    set(h,'EdgeColor','none')
    colorbar
    caxis([990 1020])
    set(gca,'FontSize',14)
    hold on
    contour(lon,lat,squeeze(land)',[1 1],'color','k','linewidth',2)
    hold off
    xlabel('Longitude (\circE)')
    ylabel('Latitude (\circN)')
    title('Mean SLP (hPa)')
subplot(1,2,2)
    [~,h]=contourf(lon,lat,spd_mean,50);
    set(h,'EdgeColor','none')
    colorbar
    caxis([0 8])
    set(gca,'FontSize',14)
    hold on
    quiver(lonn_rm(4:1:end-3,1:1:end),latt_rm(4:1:end-3,1:1:end),...
        u_mean_rm(4:1:end-3,1:1:end),v_mean_rm(4:1:end-3,1:1:end),...
        'color','k','linewidth',2)
    contour(lon,lat,-squeeze(land)',[0 0],'color','k','linewidth',2)
    hold off
    xlabel('Longitude (\circE)')
    ylabel('Latitude (\circN)')
    title('Very coarse wind field (m s^{-1})')






