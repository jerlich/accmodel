% [Pd] = one_sided_posterior(bo, Pf, cutpoint, [sign=+1])
%
%  Given a binning object and a probability distribution, zeros out the
%  probability distribution that lies to one side of a scalar cutpoint.
%  DOES NOT NORMALIZE AFTER THAT, NORMALIZATION IS LEFT FOR THE USER TO DO
%  IF DESIRED. Deals gracefully with bins and continuous values of
%  cutpoint. For example, if cutpoint is exactly at a bin center, then that
%  bin is not zeroed out, but halved (i.e.
%
%  PARAMETERS:
%  -----------
%
%  bo         a binning object
%
%  Pf         a vector representing a probability distribution. Must be the
%             same length as bins(bo)
%
%  cutpoint   a scalar. Pf to one side of this will be zeroed out.
%
%  sign       Dy default, +1. This means zero out bins that are more
%             negative than cutpoint. If sign=-1, then bins that are more
%             positive than cutpoint are zeroed out.
%


function [Pd] = one_sided_posterior(bo, Pf, cutpoint, sign)

x = bo.x;
if nargin < 4,
   sign = 1;
end;

if (cutpoint > x(end) && sign==1) || (cutpoint < x(1) && sign==0),
   % We're either asking for only bins more positive than the most positive
   % bin; or we're asking for only bins more negative than the most
   % negative bin. That means there would be nothing left!
   Pd = 0*Pf;
 %  return;
end

Pd = zeros(size(Pf));

[hp, lp] = ceil_and_floor(bo, cutpoint);


if sign==1,
    Pd(hp+1:end) = Pf(hp+1:end);
else
    Pd(1:lp-1)   = Pf(1:lp-1);
end;

if lp==hp,
    Pd(lp) = Pf(lp)/2;
else
    dh = x(hp) - cutpoint;
    dl = cutpoint - x(lp);
    dd = dh + dl;
    if sign==1,
        Pd(hp) = Pf(hp) * (1/2 + dh/dd/2);
        Pd(lp) = Pf(lp) * (dh/dd/2);
    else
        Pd(hp) = Pf(hp) * (dl/dd/2);
        Pd(lp) = Pf(lp) * (1/2 + dl/dd/2);
    end;
end;


% if sign==1,
%    Pd(1:lp-1) = 0;
% else
%    Pd(hp+1:end) = 0;
% end;
% 
% if lp==hp,
%    Pd(lp) = Pd(lp)/2;
% else
%    dd = x(hp) - x(lp);
%    Pd(lp) = Pd(lp) * 0.5*(x(hp)-cutpoint)/dd;
%    Pd(hp) = Pd(hp) * (1 - 0.5*(cutpoint - x(lp))/dd);
% end;




