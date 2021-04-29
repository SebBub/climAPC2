function M_rm=smooth_map(M,win_lon,win_lat)

% SMOOTHED_DATA=SMOOTH_MAP(DATA,LONGITUDE_WINDOW,LATITUDE_WINDOW)
%
% Smoothes data over latitudes and longitudes
%
% M:    Input data field <lat x lon>
% win_lon: Longitude window (must not be even)
% win_lat: Latitude window (must not be even)
%
% M_rm: Smoothed data <lat x lon>


N_lon=size(M,2); N_lat=size(M,1);


M_rm=NaN(size(M));
for ilon=1+(win_lon-1)/2:N_lon-(win_lon-1)/2
    for ilat=1+(win_lat-1)/2:N_lat-(win_lat-1)/2
        M_rm(ilat,ilon)=squeeze(nanmean(nanmean(...
            M(ilat-(win_lat-1)/2:ilat+(win_lat-1)/2,...
            ilon-(win_lon-1)/2:ilon+(win_lon-1)/2),2),1));
    end
end