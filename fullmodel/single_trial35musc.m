function [LL likey output] = single_trial35musc(param, mydata, varargin)
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
    'gain_first' 0; ...
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
biased_sigma2_s = param(10);
biased_input = param(11);
biased_inatt = param(12);


%% arguments needed, if not passed in
if isempty(bo), %#ok<NODEF>
    bo = binning(dx, B);
end;

if isempty(P_init), %#ok<NODEF>
    [P_init W_i] = make_P_init3musc(bo, sigma2_i, inatt,biased_inatt);
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
    W0 = make_each_F35musc(bo, sigma2_a, sigma2_s,biased_sigma2_s, biased_input, lambda, dt, 'net_input', 0, 'tot_input', 0, ...
                       'nclicks', 0, 'total_rate', total_rate, 'var_mode', var_mode,'gain_first',0);
end;
W = make_each_F35musc(bo, sigma2_a, sigma2_s,biased_sigma2_s, biased_input, lambda, dt, 'net_input', net_input, 'tot_input', tot_input, ...
                  'nclicks', nclicks, 'total_rate', total_rate, 'W0', W0, 'var_mode', var_mode,'gain_first',0);


% ~JCE~ 
% I think I understand things up until here.  It gets a little trick after this.


%% in order to compute derivatives w.r.t. alpha and rho,
% compute dh/dalpha and dh/drho numerically by finite difference,
% where h is the magnitude of the clicks
% JCE - why is this computed numerically only for rho and alpha?  


%% forwards
Pf = zeros(nbins(bo), N);
Pf(:,1) = P_init; 

for k = 2:numel(t),
    Pf(:,k) = W(k-1).F*Pf(:,k-1);
end;

Pd = one_sided_posterior(bo, Pf(:,end), bias, pokedR);
% Take only the part of the distribution that led to the rat's choice

likey = sum(Pd); % likelihood of poking the same way the data did
LL = log(likey); % log of the likelihood

if nargout > 2,
    output.Pf = Pf;
    output.x  = bins(bo);
    output.t  = t;
end;


%% backwards
% normalize to use as initial distribution of the backwards, posterior run
% Pd = Pd/sum(Pd); 

% Life must be lived forwards, but can only be understood backwards.

% % Probably need to modify net input here to be like the net_input used in make_each_F35musc.
% net_c_input = (tot_input+net_input)/2;
% net_i_input = (tot_input-net_c_input);

% net_input = net_i_input+net_c_input*biased_input;  % This is used in the make_F_struct but not to calculate variance

% [Pb Deltak] = single_trial_backwards2(Pf, Pd, W, net_input);

% % now take care of the virtual, 0th timestep that created P_init, and the
% % derivative w.r.t. sigma2_i
% Deltak = P_init_backwards(Deltak, Pf(:,1), Pb, W_i, bo);

% if nargout > 3,
%     output.Pb = Pb;
% end;
