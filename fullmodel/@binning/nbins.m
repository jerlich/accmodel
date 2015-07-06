% [n] = nbins(bo)
%
%  Given a binning object, returns a scalar that's the total number of bins
%

function [n] = nbins(bo)

n = numel(bo.x);

