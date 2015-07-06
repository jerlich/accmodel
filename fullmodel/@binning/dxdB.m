% [dx] = dxdb(bo)
%
% Returns the derivative of each bin center position w.r.t. the Bound value
% B
%
% PARAMETERS:
% -----------
%
%   bo     A binning object
%
%
% RETURNS:
% --------
%
%   dx      a vector equal in length to bins(bo), representing the
%           derivative of each bin center position w.r.t. B
%
%
% EXAMPLE CALL:
% -------------
%
%    >> bo = binning(1, 2.3); dxdB(bo)
%
%         ans = [-2 0 0 0 0 0 2]
%



function [dx] = dxdB(bo)

x  = bo.x;
dx = zeros(size(x));

dx(end) = 2;
dx(1)   = -2;

% for the case where x(end-1)==B, and x(end)-x(end-1) == dx, there's no way
% to return the positive-side derivative: a new bin will be created, and
% that's discontinuous no matter how we cut it, so there it is. We return
% only the negative-sided derivative.



