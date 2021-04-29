% Correlation, linear regression and composite examples with ENSO and SST

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
% compute the enso index

[sstc,ssta]=climanom(sst,time_sst);

enso=mean(mean(ssta(:,lat_sst>=-6&lat_sst<=6,lon_sst>=180&lon_sst<=270),2),3);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute correlation, regression and composite

[cormap,cormap_sig]=cormap_ts(ssta,enso,0.05);
[regmap,regmap_sig]=regmap_ts(ssta,enso,0.05);

thresh=2*std(enso);
[pos_comp,neg_comp,comp_diff,comp_diff_sig]=composite(sst,enso,thresh);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figures


figure(1)
subplot(1,2,1)
    [~,h]=contourf(lon_sst,lat_sst,cormap,50);
    set(h,'EdgeColor','none')
    colorbar
    caxis([-1 1])
    set(gca,'FontSize',14)
    hold on
    contourf(lon_sst,lat_sst,-squeeze(mask)',[0 0],'facecolor',...
        [0.7 0.7 0.7],'color','k','linewidth',2)
    hold off
    xlabel('Longitude (\circE)')
    ylabel('Latitude (\circN)')
    title('Correlation')
subplot(1,2,2)
    [~,h]=contourf(lon_sst,lat_sst,regmap,50);
    set(h,'EdgeColor','none')
    colorbar
    caxis([-1 1])
    set(gca,'FontSize',14)
    hold on
    contourf(lon_sst,lat_sst,-squeeze(mask)',[0 0],'facecolor',...
        [0.7 0.7 0.7],'color','k','linewidth',2)
    hold off
    xlabel('Longitude (\circE)')
    ylabel('Latitude (\circN)')
    title('Regression')

    
    
    
figure(2)
subplot(1,2,1)
    [~,h]=contourf(lon_sst,lat_sst,comp_diff,50);
    set(h,'EdgeColor','none')
    colorbar
    caxis([-1 3])
    set(gca,'FontSize',14)
    hold on
    contourf(lon_sst,lat_sst,-squeeze(mask)',[0 0],'facecolor',...
        [0.7 0.7 0.7],'color','k','linewidth',2)
    hold off
    xlabel('Longitude (\circE)')
    ylabel('Latitude (\circN)')
    title('SST composite (\circC)')
subplot(1,2,2)
    [~,h]=contourf(lon_sst,lat_sst,comp_diff_sig,50);
    set(h,'EdgeColor','none')
    colorbar
    caxis([-1 3])
    set(gca,'FontSize',14)
    hold on
    contourf(lon_sst,lat_sst,-squeeze(mask)',[0 0],'facecolor',...
        [0.7 0.7 0.7],'color','k','linewidth',2)
    hold off
    xlabel('Longitude (\circE)')
    ylabel('Latitude (\circN)')
    title('2-sigma significant SST composite (\circC)')

