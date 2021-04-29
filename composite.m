function [pos_comp,neg_comp,comp_diff,comp_diff_sig]=composite(M,index,thresh)

% [POS_COMP,NEG_COMP,COMP_DIFF,COMP_DIFF_SIG] = COMPOSITE(M,INDEX,THRESH)
%
% M is < time x lat x lon > input data matrix
% index is < time > 1D "time series" array (not timeseries class)
% thres is the +/- threshold; just give the positive (or larger) value

pos_comp=squeeze(mean(M(index>thresh,:,:),1));
neg_comp=squeeze(mean(M(index<thresh,:,:),1));
comp_diff=pos_comp-neg_comp;

% two-sigma significance

pos_comp_se=squeeze(std(M(index>thresh,:,:),1))/sqrt(length(M(index>thresh,:,:)));
neg_comp_se=squeeze(std(M(index<thresh,:,:),1))/sqrt(length(M(index<thresh,:,:)));

comp_diff_sig=comp_diff;
comp_diff_sig( pos_comp-2*pos_comp_se <= neg_comp+2*neg_comp_se ) = NaN;

