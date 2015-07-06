% [bo] = make_point_bins(dx, B)
%
% Returns a binning object.
%
% Makes a series of points that will indicate bin centers. The first and
% last points will indicate sticky bins. No "bin edges" are made-- the edge
% between two bins is always implicity at the halfway point between their
% corresponding centers. The center bin is always at x=0; bin spacing
% (except for last and first bins) is always dx; and the position
% of the first and last bins is chosen so that |B| lies exactly at the
% midpoint between 1st (sticky) and 2nd (first real) bins, as well as
% exactly at the midpoint between last but one (last real) and last
% (sticky) bins.
%
% PARAMETERS:
% -----------
%
%   dx     The standard desired spacing between bin centers.
%
%   B      The desired boundary
%
%
% RETURNS:
% --------
%
%   bo      A binning object
%
%
% EXAMPLE CALL:
% -------------
%
%    >> struct(binning(1, 2.3))
%
%         x: [-2.6000 -2 -1 0 1 2 2.6000]
%        dx: 1
%         B: 2.3000
%



function [bo] = binning(dx, B)

bo = class(struct('x', [], 'dx', [], 'B', []), mfilename);

x = 0:dx:B;
if x(end)==B,
   x(end) = x(end)+dx;
else
   x = [x 2*B-x(end)];
end;
x = [-x(end:-1:2) x];

bo.x  = x;
bo.dx = dx;
bo.B  = B;

