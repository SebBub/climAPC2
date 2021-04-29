function M_rm=runmean_maps(M,win)

% M is input data matrix <time x lat x lon>
% win is window length (MUST be odd number...
% e.g., for 7-yr running mean in monthly data, use
% 7*12-1=83)
%
% M_rm is running mean of window length win applied to M
% <time x lat x lon>

nt=size(M,1);
M_rm=zeros(size(M));
for t=1+(win-1)/2:nt-(win-1)/2
    M_rm(t,:,:)=nanmean(M(t-(win-1)/2:t+(win-1)/2,:,:),1);
end
M_rm(1:(win-1)/2,:,:)=NaN;
M_rm(nt-(win-1)/2+1:nt,:,:)=NaN;

