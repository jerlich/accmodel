% [hp, lp] = ceil_and_floor(bo, s)
%
% Given a binning object bo and a vector of scalar values s, returns the
% bin numbers for the bin just above and just below each element of s.  If
% s is exactly at a bin center, then the two bin ids returned are equal
%
% PARAMETERS:
% -----------
% 
%  bo          a binning object
%
%  s           a vector of scalars
%
%
% RETURNS:
% --------
%
% hp           same size as s
%              the bin number that has the smallest bin value that is
%              larger or equal to s
%
% lp           same size as s
%              the bin number that has the largest bin value that is
%              smaller or equal to s
%


function [hp, lp] = ceil_and_floor(bo, s)

x  = bo.x;
dx = bo.dx;
n  = numel(x);

hp = ceil( (s-x(2))/dx)+2; 
lp = floor((s-x(2))/dx)+2; 

hp(hp<1) = 1; hp(hp>n) = n;
lp(lp<1) = 1; lp(lp>n) = n;

hp(x(1)<s     & s<x(2))   = 2;
lp(x(end-1)<s & s<x(end)) = n-1;

hp(x(end)<s) = n; lp(x(end)<s) = n;
hp(s<x(1)) = 1;   lp(s<x(1)) = 1;


% 
% if s <= x(1), 
%    hp = 1;
%    lp = 1;
% elseif x(end) <= s,
%    hp = numel(x);
%    lp = numel(x);
% elseif x(1) < s && s < x(2),
%    lp = 1;          hp = 2;
% elseif x(end-1) < s && s < x(end),
%    lp = numel(x)-1; hp = numel(x);
% else
%    hp = ceil( (s-x(2))/dx) + 2;
%    lp = floor((s-x(2))/dx) + 2;
% end;
% 
