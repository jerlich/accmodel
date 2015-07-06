function [LL dLL likey output] = single_trial35(param, mydata, varargin)
% this version scales the sensory noise with the magnitude of the clicks,
% June 2012, BWB


pairs = { ...
    'dx'          0.25  ; ...
    'dt'          0.02  ; ...
    'bo'            []  ; ...
    'P_init'        []  ; ...
    'W_i'           []  ; ...
    'W0'            []  ; ...
    'ndelta'      1e-4  ; ...
    'total_rate'    40  ; ... % total Poisson event rate per second
    'var_mode'       1  ; ... % if 0, sensory variance proportional to nclicks; if 1, sensory variance proportional to total input
    'adapt_mode'     0  ; ... % if 0, same-side adaptation; if 1, across-side adaptation
    'kappa'      [1 1]  ; ... % the biased weight of [left right] clicks
    'kappa_timerange'   [0 1000] ; ... % the time range over which kappa values are applied
    };
parseargs(varargin, pairs);

%% unpack click inputs for this trial
mydata = mydata(1);
T         = mydata.T;
leftbups  = mydata.leftbups;  leftbups = leftbups(leftbups < T);
rightbups = mydata.rightbups; rightbups = rightbups(rightbups < T);
pokedR    = mydata.pokedR;

%% unpack model parameters
lambda   = param(1);
sigma2_a = param(2);
sigma2_s = param(3);
sigma2_i = param(4);
B        = param(5);
alpha    = param(6);
rho      = param(7);
bias     = param(8);
inatt    = param(9);

%% arguments needed, if not passed in
if isempty(bo), %#ok<NODEF>
    bo = binning(dx, B);
end;

if isempty(P_init), %#ok<NODEF>
    [P_init W_i] = make_P_init3(bo, sigma2_i, inatt);
end;

%% make Markov matrices for each time step
N  = ceil(T/dt);  % number of timesteps
t  = (0:N-1)*dt;  % a time axis

% clicks_L and clicks_R are the same size as leftbups and rightbups, respectively
if adapt_mode == 1,
    [clicks_L clicks_R] = make_adapted_cat_clicks(leftbups, rightbups, alpha, rho);
else
    [clicks_L clicks_R] = make_adapted_clicks(leftbups, rightbups, alpha, rho, 1);
end;

if abs(kappa(1) - 1) > eps,
    clicks_L(clicks_L>=kappa_timerange(1) & clicks_L<=kappa_timerange(2)) = clicks_L(clicks_L>=kappa_timerange(1) & clicks_L<=kappa_timerange(2)) * kappa(1);
end;
if abs(kappa(2) - 1) > eps,
    clicks_R(clicks_R>=kappa_timerange(1) & clicks_R<=kappa_timerange(2)) = clicks_R(clicks_R>=kappa_timerange(1) & clicks_R<=kappa_timerange(2)) * kappa(2);
end;

% net_input and nclicks are the size size as t
[net_input tot_input nclicks] = make_click_inputs35(t, leftbups, rightbups, clicks_L, clicks_R);

% W is a cell array with one element for each time step, containing F and
% all the dFdtheta's used in the backwards, posterior run
if isempty(W0), %#ok<NODEF>
    W0 = make_each_F35(bo, sigma2_a, sigma2_s, lambda, dt, 'net_input', 0, 'tot_input', 0, ...
                       'nclicks', 0, 'total_rate', total_rate, 'var_mode', var_mode);
end;
W = make_each_F35(bo, sigma2_a, sigma2_s, lambda, dt, 'net_input', net_input, 'tot_input', tot_input, ...
                  'nclicks', nclicks, 'total_rate', total_rate, 'W0', W0, 'var_mode', var_mode);


Pf = zeros(nbins(bo), N);
Pf(:,1) = P_init; 

for k = 2:numel(t),
    Pf(:,k) = W(k-1).F*Pf(:,k-1);
end;

Pd = one_sided_posterior(bo, Pf(:,end), bias, pokedR);
likey = sum(Pd); % likelihood of poking the same way the data did
LL = log(likey); % log of the likelihood


if nargout==1
    % We just want the LL not the derivative
    return;
end


%% in order to compute derivatives w.r.t. alpha and rho,
% compute dh/dalpha and dh/drho numerically by finite difference,
% where h is the magnitude of the clicks
if adapt_mode == 1,
    [clicks_L_dalpha clicks_R_dalpha] = make_adapted_cat_clicks(leftbups, rightbups, alpha+ndelta, rho);
else
    [clicks_L_dalpha clicks_R_dalpha] = make_adapted_clicks(leftbups, rightbups, alpha+ndelta, rho, 1);
end;
[net_input_dalpha tot_input_dalpha] = make_click_inputs35(t, leftbups, rightbups, clicks_L_dalpha, clicks_R_dalpha);
net_input_dalpha = (net_input_dalpha - net_input)/ndelta;
tot_input_dalpha = (tot_input_dalpha - tot_input)/ndelta;

if adapt_mode == 1,
    [clicks_L_drho clicks_R_drho] = make_adapted_cat_clicks(leftbups, rightbups, alpha, rho+ndelta);
else
    [clicks_L_drho clicks_R_drho] = make_adapted_clicks(leftbups, rightbups, alpha, rho+ndelta, 1);
end;
[net_input_drho tot_input_drho] = make_click_inputs35(t, leftbups, rightbups, clicks_L_drho, clicks_R_drho);
net_input_drho = (net_input_drho - net_input)/ndelta;
tot_input_drho = (tot_input_drho - tot_input)/ndelta;

%% forwards

if nargout > 3,
    output.Pf = Pf;
    output.x  = bins(bo);
    output.t  = t;
end;


%% backwards
% normalize to use as initial distribution of the backwards, posterior run
Pd = Pd/sum(Pd); 

% Life must be lived forwards, but can only be understood backwards.
[Pb Deltak] = single_trial_backwards2(Pf, Pd, W, net_input);

% now take care of the virtual, 0th timestep that created P_init, and the
% derivative w.r.t. sigma2_i
Deltak = P_init_backwards(Deltak, Pf(:,1), Pb, W_i, bo);

if nargout > 3,
    output.Pb = Pb;
end;

%% compute gradient

% lambda
    dLL(1) = sum(Deltak.l);

% sigma2_a
    % dL/dsigma2_a = dL/dsig2 * dsig2/dsigma2_a
    %              = dL/dsig2 * dt

    dLL(2) = sum(Deltak.sig2) * dt;

% sigma2_s
    if var_mode == 1,
        % dL/dsigma2_s = dL/dsig2 * dsig2/dsigma2_s
        %              = dL/dsig2 * tot_input / total_rate

        dLL(3) = sum(Deltak.sig2.*tot_input')/total_rate;
    else
        % dL/dsigma2_s = dL/dsigma2 * dsigma2/dsigma2_s
        %              = dL/dsigma2 * nclicks / total_rate

        dLL(3) = sum(Deltak.sig2.*nclicks')/total_rate; 
    end;
% sigma2_i
    dLL(4) = Deltak.sigma2_i;
    
% B
    dLL(5) = sum(Deltak.B);

% alpha
    % dL/dalpha = dL/dh * dh/dalpha + dL/dsig2 * dsig2/dalpha
    %           = dL/dh * dh/dalpha + dL/dsig2 * sigma2_s/total_rate * dH/dalpha

    dLL(6) = sum(Deltak.h.*net_input_dalpha')/dt;
    
    if var_mode == 1,
        dLL(6) = dLL(6) + sum(Deltak.sig2.* tot_input_dalpha')*sigma2_s/total_rate; 
    end;

% rho
    % dL/drho = dL/dh * dh/drho + dL/dsig2 * dsig2/drho
    %         = dL/dh * dh/drho + dL/dsig2 * sigma2_s/total_rate * dH/drho

    dLL(7) = sum(Deltak.h.*net_input_drho')/dt;
    
    if var_mode == 1,
        dLL(7) = dLL(7) +sum(Deltak.sig2.*tot_input_drho')*sigma2_s/total_rate;      
    end;

% bias
    % compute deriv numerically,
    LLepsilon = one_sided_posterior(bo,Pf(:,end),bias+ndelta,pokedR);
    dLL(8) = (log(sum(LLepsilon)) - LL)/ndelta;

% inatt
    % compute deriv numerically,
    dLL(9) = (log((likey-inatt/2)/(1-inatt)*(1-(inatt+ndelta))+(inatt+ndelta)/2) - log(likey))/ndelta;
    
    
% single trial done!    