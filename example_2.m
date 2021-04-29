% Example script to subtract the mean over a latitude circule from the time
% series

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
% Compute the range as a function of latitude.

% Compute the climatological means and anomalies
[sst_seas,sst_ano]=climanom(sst,time_sst);

% Compute the range or the mean:
sst_mean=squeeze(mean(sst,1));
sst_range=squeeze(range(sst_seas,1));

% Average over all longitudes
sst_mean_lat = nanmean(sst_mean,2);
sst_range_lat = nanmean(sst_range,2);

% Subtract the mean as a function of latitude from  SST
sst_lat_ano=sst_mean-repmat(sst_mean_lat,[1 length(lon_sst)]);
sst_lat_ano_range=sst_range-repmat(sst_range_lat,[1 length(lon_sst)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figures

figure(1)
[~,h]=contourf(lon_sst,lat_sst,sst_range,50);
set(h,'EdgeColor','none')
colorbar
caxis([0 12])
set(gca,'FontSize',14)
hold on
contourf(lon_sst,lat_sst,-squeeze(mask)',[0 0],'facecolor',[0.7 0.7 0.7],'color','k','linewidth',2)
%contour(lon_sst,lat_sst,-squeeze(mask)',[0 0],'color','k','linewidth',2)
hold off
xlabel('Longitude (\circE)')
ylabel('Latitude (\circN)')
title('Seasonal SST range (\circC)')



figure(2)
subplot(1,2,1)
    plot(sst_mean_lat,lat_sst,'linewidth',2)
    %plot(sst_range_lat,lat_sst,'linewidth',2)
    set(gca,'FontSize',14)
    box on
    grid on
    xlabel('Mean SST (\circC)')
    ylabel('Latitude (\circN)')
    %title('Mean latitudinal SST range')
    title('Mean latitudinal SST')
subplot(1,2,2)
    [~,h]=contourf(lon_sst,lat_sst,sst_lat_ano,50);
    %[~,h]=contourf(lon_sst,lat_sst,sst_lat_ano_range,50);
    set(h,'EdgeColor','none')
    colorbar
    caxis([-6 6])
    set(gca,'FontSize',14)
    hold on
    contourf(lon_sst,lat_sst,-squeeze(mask)',[0 0],'facecolor',[0.7 0.7 0.7],'color','k','linewidth',2)
    %contour(lon_sst,lat_sst,-squeeze(mask)',[0 0],'color','k','linewidth',2)
    hold off
    xlabel('Longitude (\circE)')
    ylabel('Latitude (\circN)')
    title('SST with the mean latitudinal SST subtracted')


% If you are already very familiar with Matlab you can use cool diverging
% colorbars for this. You can find many nice colormaps online...


