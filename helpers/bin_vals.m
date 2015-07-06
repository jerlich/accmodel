function bin = bin_vals(vals, n_bins, varargin)
% function [bin, bin_list] = bin_vals(vals, n_bins, varargin)
%
% Bins input in quantiles
% 
% INPUT:
%
%   vals:       Values to bin
%   n_bins:     Number of bins
%
% OUTPUT:
%
%   bin:        vector of bin numbers, one for each element vals input

%% Default Parameters
split_zero = false;         % separately bin vals on either side of zero (values == 0 at randomly split between the two)

% override based on varargin
opts = overridedefaults(who, varargin);

if split_zero && mod(n_bins,2)==1
    error('Must be even number of bins to split at zero')
end

if split_zero
    n_bins = n_bins/2;
    % indices for positive and negative vals
    pos_ind = vals>0;
    zero_find = find(vals==0);
    
    % randomly assign half zero vals to positive and half to negative
    % first, reset random number generator to same start for consistency
    % across runs
%     RandStream.setDefaultStream(RandStream('mrg32k3a','Seed',10));
    RandStream.setGlobalStream(RandStream('mrg32k3a','Seed',10));
    
    % shuffle zero positions
    zero_shuffle = zero_find(randperm(length(zero_find)));
    % set first half after shuffling to positive
    zero_to_pos = zero_shuffle(1:round(end/2));
    pos_ind(zero_to_pos) = true;
    
    % make negative the rest of the trials
    neg_ind = ~pos_ind;
    
    % separately bin positive and negative
    signed_vals = vals(neg_ind);
    neg_bin = ceil(n_bins * untiedrank(signed_vals) / length(signed_vals));
    signed_vals = vals(pos_ind);
    pos_bin = ceil(n_bins * untiedrank(signed_vals) / length(signed_vals)) + n_bins;
    
    % put them together
    bin = nan(size(vals));
    bin(neg_ind) = neg_bin;
    bin(pos_ind) = pos_bin;
    
else
    % puts same number of trials into each bin
    bin = ceil(n_bins * untiedrank(vals) / length(vals));
end    
        
        