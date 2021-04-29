function [regmap,regmap_sig]=regmap_ts(field,ts,alpha)

% This function regresses a time-varying field <time x lat x lon> onto a
% time series. The time dimension of the field must be the same length as
% the time series. Use alpha=0.05 for 95% confidence.

ts=normalize(ts);
ts=[ones(size(ts)) ts];

field(isnan(field))=0;

regmap=zeros(size(field,2),size(field,3));
regmap_sig=zeros(size(regmap));

for y=1:size(field,2)
    for x=1:size(field,3)
        [b,bint]=regress(field(:,y,x),ts,alpha);
        regmap(y,x)=b(2);
        if sign(bint(2,1))==sign(bint(2,2))
            regmap_sig(y,x)=b(2);
        else
            regmap_sig(y,x)=NaN;
        end
    end
end
