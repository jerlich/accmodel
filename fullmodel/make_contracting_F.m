function [F dFds dFdsig2 dFdB dFdl dFdh] = make_contracting_F(bo, bo2)

x_exp = bins(bo2);
n_real_bins = nbins(bo);
n_extra_bins = (nbins(bo2) - n_real_bins)/2;

F  = eye(numel(x_exp), numel(x_exp));
F(:  , 1:n_extra_bins+1) = 0;
F(1  , 1:n_extra_bins+1) = 1;
F(:  , n_extra_bins+n_real_bins:end) = 0;
F(end, n_extra_bins+n_real_bins:end) = 1;

dFds    = 0;
dFdsig2 = 0;
dFdB    = 0;
dFdl    = 0;
dFdh    = 0;