function [Mc,Ma]=climanom(M,time)

% [CLIM,ANOM]=CLIMANOM(MONTHLY_DATA);
%
% Computes mean climatology and anomalies
%
% M:    Input monthly data matrix <time x lat x lon>
% time: Time vector of the input data
%
% Mc: Mean climatology <12 x lat x lon>
% Ma: Anomalies <time x lat x lon>


% compute monthly climatology
dt=datevec(time); months=dt(:,2);
Mc=zeros(12,size(M,2),size(M,3));
Ma=zeros(size(M));
for icalmo=1:12
    ft=find(months==icalmo);
    Mc(icalmo,:,:)=nanmean(M(ft,:,:),1);
    for iy=1:length(ft)
        Ma(ft(iy),:,:)=squeeze(M(ft(iy),:,:))-squeeze(Mc(icalmo,:,:));
    end
end