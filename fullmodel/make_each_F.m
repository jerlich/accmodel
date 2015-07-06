% function [W] = make_each_F(bo, sigma2_a, sigma2_s, lambda, dt, varargin)

function [W ticking] = make_each_F(bo, sigma2_a, sigma2_s, lambda, dt, varargin)

pairs = { ...
    'net_input'     0  ; ...
    'nclicks'       0  ; ...
    'W0'            [] ; ...
    };
parseargs(varargin, pairs);

W = cell(numel(net_input), 1);

% if isempty(W0), %#ok<NODEF>
%     W0 = make_each_F(bo, sigma2_a, sigma2_s, lambda, dt, 'net_input', 0, 'nclicks', 0);
%     W0 = W0{1};
% end;

% ticking = zeros(size(net_input));
for k = 1:numel(net_input),
%     tic;
    if nclicks(k) == 0 && ~isempty(W0),
        W{k} = W0;
    else
        sigma2_tot = tot_var(sigma2_a, sigma2_s, nclicks(k), dt);

        % gamma = exp(lambda*dt)
        % phi = c/lambda = net_input/dt/lambda
        [F dFds dFdsig2 dFdB dFdg dFdphi] = make_F(bo, sigma2_tot, exp(lambda*dt), net_input(k)/dt/lambda);

        thisW.F = F;
        thisW.dFds    = dFds;
        thisW.dFdsig2 = dFdsig2;
        thisW.dFdB    = dFdB;
        thisW.dFdg    = dFdg;
        thisW.dFdphi  = dFdphi;

        W{k} = thisW;
    end;
%     ticking(k) = toc;
end;