function [P_init W_i] = make_P_init(bo, sigma2_i, inatt)

x = bins(bo);
P_0 = zeros(size(x))';

P_0(ceil_and_floor(bo, 0)) = 1-inatt;
P_0(1)   = inatt/2;
P_0(end) = inatt/2;

[F_init, ~, dFdsig2, dFdB]= make_F(bo, sigma2_i, 1, 0);

P_init = F_init * P_0;

W_i.P0      = P_0;
W_i.F       = F_init;
W_i.dFdsig2 = dFdsig2;
W_i.dFdB    = dFdB;