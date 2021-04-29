function [rmap,pmap,smap]=cormap_ts(field,ts,siglev)

% This function computes the correlation between a time-varying field
% <time x lat x lon> and a single time series.
% The time (1st) dimension of the field must be of the same length
% as the length of the time series.
% Use siglev=0.05 for 95% confidence.

rmap=zeros(size(field,2),size(field,3));
pmap=zeros(size(field,2),size(field,3));
smap=zeros(size(field,2),size(field,3));

for y=1:size(field,2)
    for x=1:size(field,3)
        [cor,pval]=corrcoef(field(:,y,x),ts);
        rmap(y,x)=cor(1,2);
        pmap(y,x)=pval(1,2);
        if pmap(y,x)>=siglev
            smap(y,x)=NaN;
        else
            smap(y,x)=rmap(y,x);
        end
    end
end
