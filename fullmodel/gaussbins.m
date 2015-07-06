% function [deltas ps] = gaussbins(dx, sigma2)
%
% returns slices of a gaussian with sigma2 variance, sampled in a way
% that's appropriate given bin size dx
%
% if sigma2 is small enough, then 50 subdivisions is enough
% if sigma2 is too big, then the subdivisions defining the gaussian need
% to be closer together than dx

function [deltas ps] = gaussbins(dx, sigma2)

if sigma2 < eps,
    deltas = 0;
    ps = 1;
else
    n = ceil(2*5*sqrt(sigma2)/dx);
    n = max([20 n]);

    deltas = (-n:n)*(5*sqrt(sigma2))/n;
    ps = exp(-deltas.^2/(2*sigma2));
    ps = ps/sum(ps);
end;