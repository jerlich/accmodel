function [W] = make_each_F35(bo, sigma2_a, sigma2_s, lambda, dt, varargin)
% June 2012: var_s scaled by magnitude of each click

pairs = { ...
    'net_input'     0  ; ...
    'tot_input'     0  ; ...
    'nclicks'       0  ; ...
    'total_rate'   40  ; ...
    'var_mode'      1  ; ... % if 0, sensory variance proportional to nclicks; if 1, sensory variance proportional to total input
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
    if tot_input(k) == 0 && ~isempty(W0),
        W(k) = W0;
    else
        var_a = sigma2_a*dt;
        
        if var_mode == 1,
            var_s = tot_input(k)*sigma2_s/total_rate;
        else
            var_s = nclicks(k)*sigma2_s/total_rate;
        end;

        sigma2_tot = var_a + var_s;

        W(k) = make_F_struct(x, dB, sigma2_tot, lambda, net_input(k)/dt, dt, dx); 
    end;
end;