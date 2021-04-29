% Example script to find locations and times of minimum and maximum values

clear all, close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load the data


file_name='sst.mnmean.nc';

sst=ncread(file_name,'sst');
lon_sst=ncread(file_name,'lon');
lat_sst=ncread(file_name,'lat');
time_sst=ncread(file_name,'time');

time_sst=time_sst+datenum(1800,1,1);
time2_sst=datestr(time_sst,'mmm yyyy');

mask=ncread('lsmask_sst.nc','mask');
mask3=repmat(mask,[1 1 size(sst,3)]);
sst(mask3==0)=NaN; clear mask3

sst=permute(sst,[3 2 1]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find the maximum. You can easily modify the code below to find the
% location and time of the minimum...


% Average over three months:
sst_rm=runmean_maps(sst,3);

% Identify the maximum at each location
[max_sst] = squeeze(max(sst_rm,[],1));

% Find the maximum over all longitudes as a function of latitude:
[max_sst_lat] = max(max_sst,[],2);

% Get the latitude of the maximum temperature:
[~,ilat]=max(max_sst_lat); max_lat=lat_sst(ilat);

% Then you find the maximum longitude with:
max_sst_lon=squeeze(max_sst(ilat,:));
[~,ilon]=max(max_sst_lon); max_lon=lon_sst(ilon);

% ... and the time:
max_sst_time=squeeze(sst_rm(:,ilat,ilon));
[~,itime]=max(max_sst_time); max_time=time_sst(itime);

% check if you have really found the right location and time:
check=isequal(squeeze(max(max(max(sst_rm)))),...
    squeeze(sst_rm(itime,ilat,ilon)));

if check==1
    disp(['The maximum occurred at Latitude = ' num2str(max_lat) ' ' ...
        char(176) 'N and Longitude = ' ...
        num2str(max_lon) ' ' char(176) 'E around ' time2_sst(itime,:) '.'])
elseif check==0
    disp('Something went wrong.')
end


% See what country or ocean is at that location:


figure(1)
[~,h]=contourf(lon_sst,lat_sst,max_sst,50);
set(h,'EdgeColor','none')
colorbar
caxis([0 30])
set(gca,'FontSize',14)
hold on
contourf(lon_sst,lat_sst,-squeeze(mask)',[0 0],'facecolor',[0.7 0.7 0.7],'color','k','linewidth',2)
%contour(lon_sst,lat_sst,-squeeze(mask)',[0 0],'color','k','linewidth',2)
plot(max_lon,max_lat,'rx','linewidth',4,'markersize',12)
hold off
xlabel('Longitude (\circE)')
ylabel('Latitude (\circN)')
title('Maximum SST (\circC)')


% If you would like to add a coastline to the air temperature and
% precipitation plots, you can use the 'land.nc' data. If you need help,
% just ask me.





