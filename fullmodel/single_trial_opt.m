function [LL ] = single_trial_opt(Wf,param, pokedR,leftbups,rightbups,T, varargin)
% this version scales the sensory noise with the magnitude of the clicks,
% April 2014, JCE.  Modified from June 2012, BWB
% 
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


lambda   = param(1);
sigma2_a = param(2);
sigma2_s = param(3);
sigma2_i = param(4);
B        = param(5);
alpha    = param(6);
rho      = param(7);
bias     = param(8);
inatt    = param(9);
if numel(param)==9
    biased_sigma2_s = param(3);
    biased_input = 1;
    biased_inatt = param(9);

else
    %% unpack model parameters

    biased_sigma2_s = param(10);
    biased_input = param(11);
    biased_inatt = param(12);
end
%% unpack model parameters


%% arguments needed, if not passed in
if isempty(bo), %#ok<NODEF>
    bo = binning(dx, B);
end;

%% make Markov matrices for each time step
N  = ceil(T/dt);  % number of timesteps
t  = (0:N-1)*dt;  % a time axis

% clicks_L and clicks_R are the same size as leftbups and rightbups, respectively
if adapt_mode == 1,
    [clicks_L, clicks_R] = make_adapted_cat_clicks(leftbups, rightbups, alpha, rho);
else
    [clicks_L, clicks_R] = make_adapted_clicks(leftbups, rightbups, alpha, rho, 1);
end;


% net_input and nclicks are the size size as t
[net_input, tot_input, nclicks] = make_click_inputs35(t, leftbups, rightbups, clicks_L, clicks_R);


W = Wf(net_input, tot_input,nclicks);



Pf = zeros(nbins(bo), N);
Pf(:,1) = P_init; 

for k = 2:numel(t),
    Pf(:,k) = W(k-1).F*Pf(:,k-1);
end;

Pd = one_sided_posterior(bo, Pf(:,end), bias, pokedR);
likey = sum(Pd); % likelihood of poking the same way the data did
LL = log(likey); % log of the likelihood


