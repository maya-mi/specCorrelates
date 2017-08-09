function [ns, scal, incpt] = normSeries(vec)
incpt = nanmean(vec);
vec = vec - incpt;
scal = max(abs(vec));
ns = vec /scal;  
end