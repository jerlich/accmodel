% function [W] = make_each_F2(bo, sigma2_a, sigma2_s, lambda, dt, varargin)

function [W] = make_each_F2(bo, sigma2_a, sigma2_s, lambda, dt, varargin)

pairs = { ...
    'net_input'     0  ; ...
    'nclicks'       0  ; ...
    'W0'            [] ; ...
    };
parseargs(varargin, pairs);

W = cell(numel(net_input), 1);

for k = 1:numel(net_input),
    if nclicks(k) == 0 && ~isempty(W0),
        W{k} = W0;
    else
        sigma2_tot = tot_var(sigma2_a, sigma2_s, nclicks(k), dt);

        [F dFds dFdsig2 dFdB dFdl dFdh] = make_F2_mex( ...
                bins(bo), dxdB(bo), sigma2_tot, lambda, ...
                net_input(k)/dt, dt, binsize(bo)); %#ok<ASGLU>

        thisW.F = F;
        thisW.dFdsig2 = dFdsig2;
        thisW.dFdB    = dFdB;
        thisW.dFdl    = dFdl;
        thisW.dFdh    = dFdh;

        W{k} = thisW;
    end;
end;