function [W] = make_each_F35musc(bo, sigma2_a, sigma2_s, biased_sigma2_s, biased_input, lambda, dt, varargin)
% June 2012: var_s scaled by magnitude of each click
% July 2013, include biased_sigma2_s

% The current code will have less noise for more depressed clicks or more "silenced" clicks.

pairs = { ...
    'net_input'     0  ; ...
    'tot_input'     0  ; ...
    'nclicks'       0  ; ...  % ipsi clicks
     'total_rate'   40  ; ...
    'var_mode'      1  ; ... % if 0, sensory variance proportional to nclicks; if 1, sensory variance proportional to total input
    'W0'            [] ; ...
    'gain_first'     0; ...
    };
parseargs(varargin, pairs);

W = struct('F', [], ...
           'dFdsig2',  [], 'dFdB', [], ...
           'dFdl', [], 'dFdh', []);
x = bins(bo);
dB = dxdB(bo);  
dx = binsize(bo);


% net_input   = net_c_input - net_i_input
% tot_input   = net_c_input + net_i_input

% net_i_input = net_c_input - net_input
% tot_input   = net_c_input + net_c_input - net_input
% net_c_input = (tot_input + net_input)/2


net_c_input = (tot_input+net_input)/2;
net_i_input = (tot_input-net_c_input);

net_input = net_c_input*biased_input-net_i_input;  % This is used in the make_F_struct but not to calculate variance

for k = 1:numel(net_input),

    if tot_input(k) == 0 && ~isempty(W0),
        W(k) = W0;
    else
        var_a = sigma2_a*dt;
        
        if var_mode == 1,
            if gain_first==1,
                % variance is based on the click input after scaling
                var_s = net_i_input(k)*sigma2_s/total_rate + net_c_input(k)*biased_input*biased_sigma2_s/total_rate;
             else
                % The variance is based on the click input without scaling.
                var_s = net_i_input(k)*sigma2_s/total_rate + net_c_input(k)*biased_sigma2_s/total_rate;
            end
        else
            error('not implemented');
            var_s = nclicks(k)*sigma2_s/total_rate;
        end;

        sigma2_tot = var_a + var_s;

        W(k) = make_F_struct(x, dB, sigma2_tot, lambda, net_input(k)/dt, dt, dx); 
    end;
end;


