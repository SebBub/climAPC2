function [ts_normed] = normalize(ts)

ts_normed = (ts-mean(ts))/std(ts);