function M_rm=runmean(M,win)

% M is input data data vector
% win is window length (MUST be odd number...
% e.g., for 7-yr running mean in monthly data, use
% 7*12-1=83)
%
% M_rm is running mean of window length win applied to M

nt=length(M);
M_rm=zeros(nt,1);
for t=1+(win-1)/2:nt-(win-1)/2
    M_rm(t)=squeeze(nanmean(M(t-(win-1)/2:t+(win-1)/2)));
end
M_rm(1:(win-1)/2)=NaN;
M_rm(nt-(win-1)/2+1:nt)=NaN;

