function [P_init W_i] = make_P_init3(bo, sigma2_i, inatt)

x = bins(bo);
P_0 = zeros(size(x))';

P_0(ceil_and_floor(bo, 0)) = 1-inatt;
P_0(1)   = inatt/2;
P_0(end) = inatt/2;

W_i = make_F_struct(x, dxdB(bo), sigma2_i, 0, 0, 0, binsize(bo)); 

P_init = W_i.F * P_0;
W_i.P0 = P_0;