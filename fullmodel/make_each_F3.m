% function [W] = make_each_F3(bo, sigma2_a, sigma2_s, lambda, dt, varargin)

function [W] = make_each_F3(bo, sigma2_a, sigma2_s, lambda, dt, varargin)

pairs = { ...
    'net_input'     0  ; ...
    'nclicks'       0  ; ...
    'W0'            [] ; ...
    };
parseargs(varargin, pairs);

W = struct('F', [], ...
           'dFdsig2',  [], 'dFdB', [], ...
           'dFdl', [], 'dFdh', []);
x = bins(bo);
dB = dxdB(bo);
dx = binsize(bo);

for k = 1:numel(net_input),
    if nclicks(k) == 0 && ~isempty(W0),
        W(k) = W0;
    else
        sigma2_tot = tot_var(sigma2_a, sigma2_s, nclicks(k), dt);

        W(k) = make_F_struct(x, dB, sigma2_tot, lambda, net_input(k)/dt, dt, dx); 
    end;
end;