function P = collapse_P(P_exp, bo, bo_exp)

n_extra_bins = (nbins(bo_exp)-nbins(bo))/2;
P = zeros(nbins(bo), 1);
P(1)       = sum(P_exp(1:n_extra_bins+1));
P(2:end-1) = P_exp(n_extra_bins+2:n_extra_bins+nbins(bo)-1);
P(end)     = sum(P_exp(n_extra_bins+nbins(bo):end));