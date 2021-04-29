% Example script to compute trends

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

% Pick a location in the Baltic (ignoring the coarse resolution...)
lat_baltic=54.31; lon_baltic=360-10.13;
[~,ilat]=min(abs(lat_sst-lat_baltic));
[~,ilon]=min(abs(lon_sst-lon_baltic));

sst_baltic=squeeze(sst_ano(:,ilat,ilon));

% While there are many ways to compute trends (e.g. by using a linear fit),
% one convenient way in Matlab is to just regress the sst on time. (We will
% talk more about regressions in the next class...)

X=ones(length(time_sst),2); X(:,2)=time_sst;
[sst_trend,trend_int]=regress(sst_baltic,X);

%If there is a significant trend at the 95% confidence level, both values
%in trend_int(2,:) should have the same sign.

pm=(trend_int(2,2)-trend_int(2,1))/2;

if ( abs(sst_trend(2)) <= abs(pm) )
    disp('There was no significant SST trend at the 95% confidence level.')
else
    sst_trend_full=sst_trend(2)*(time_sst(end)-time_sst(1));
    switch sign(sst_trend(2))
        case 1
            disp(['The SST increased by ' sprintf('%1.2e',sst_trend(2)) ...
                ' ' char(176) 'C per day or by ' sprintf('%1.2f',...
                sst_trend_full) ' ' char(176) 'C between ' ...
                time2_sst(1,:) ' and ' time2_sst(end,:) '.'])
        case -1
            disp(['The SST decreased by ' sprintf('%1.2e',sst_trend(2)) ...
                ' ' char(176) 'C per day or by ' sprintf('%1.2f',...
                sst_trend_full) ' ' char(176) 'C between ' ...
                time2_sst(1,:) ' and ' time2_sst(end,:) '.'])
    end
end


% Show a map of the trends
sst_trend_map=sigtrendmap(sst_ano,time_sst);

figure(1)
[~,h]=contourf(lon_sst,lat_sst,sst_trend_map,50);
set(h,'EdgeColor','none')
colorbar
caxis([-1 1]*1e-4)
set(gca,'FontSize',14)
hold on
contourf(lon_sst,lat_sst,-squeeze(mask)',[0 0],'facecolor',[0.7 0.7 0.7],'color','k','linewidth',2)
%contour(lon_sst,lat_sst,-squeeze(mask)',[0 0],'color','k','linewidth',2)
hold off
xlabel('Longitude (\circE)')
ylabel('Latitude (\circN)')
title('SST trends (\circC day^{-1})')

