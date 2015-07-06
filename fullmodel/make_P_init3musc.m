function [P_init W_i] = make_P_init3musc(bo, sigma2_i, inatt,biased_inatt)

% Makes the t=0 probability distribution.

x = bins(bo);
P_0 = zeros(size(x))';

P_0(ceil_and_floor(bo, 0)) = 1-inatt/2-biased_inatt/2;

% P_0(1) correspond to ipsi, or negative input
% P_0(end) correspond to contra, or posistive input

% Our bias leads to more ipsi choices so we need it to be an increase in P_0(1)

P_0(1) = biased_inatt/2;
P_0(end)   = inatt/2;

W_i = make_F_struct(x, dxdB(bo), sigma2_i, 0, 0, 0, binsize(bo)); 

P_init = W_i.F * P_0;
W_i.P0 = P_0;