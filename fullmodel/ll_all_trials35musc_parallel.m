function [LL likey out] = ll_all_trials35musc_parallel(param, rawdata, varargin)

global history3

% You need to define these here because parfor doesn't handle parseargs well.
dx = 0.25; dt = 0.02; total_rate = 40; var_mode = 1; adapt_mode = 0; kappa = [1 1]; kappa_timerange = [0 1000]; gain_first=0;

pairs = { ...
    'for_fmin'             0  ; ...
    'dx'                0.25  ; ...
    'dt'                0.02  ; ...
    'total_rate'           40 ; ... % total rate of Poisson events per second
    'var_mode'             1  ; ...
    'adapt_mode'           0  ; ...
    'track_history'        0  ; ...
    'kappa'            [1 1]  ; ... % the biased weight of [left right] clicks    
    'kappa_timerange'   [0 1000] ; ... % the time range over which kappa values are applied
    'gain_first'  0;...
    };
parseargs(varargin, pairs);

%% only run if matlabpool is open
% 
% if matlabpool('size')<1,
%     fprintf(1, 'll_all_trials3_parallel: no matlabpool found, not running this function.\n');
%     LL = NaN;
%     dLL = NaN * zeros(size(param));
%     return;
% end;


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

% used to keep track of reflections so that the gradient can be made to
% point in the right direction
reflections = ones(size(param)); 

if sigma2_a < 0, sigma2_a = abs(sigma2_a); reflections(2) = -1; end;
if sigma2_s < 0, sigma2_s = abs(sigma2_s); reflections(3) = -1; end;
if biased_sigma2_s < 0, biased_sigma2_s = abs(biased_sigma2_s); reflections(3) = -1; end;

if sigma2_i < 0, sigma2_i = abs(sigma2_i); reflections(4) = -1; end;
if alpha < 0, alpha = abs(alpha); reflections(6) = -1; end;
if rho <0, rho = abs(rho); reflections(7) = -1; end;
if inatt < 0, inatt = abs(inatt); reflections(9) = -1; end; 
if inatt > 1, inatt = 2-inatt; reflections(9) = -1; end;

if biased_inatt < 0, biased_inatt = abs(biased_inatt); reflections(9) = -1; end; 
if biased_inatt > 1.6, biased_inatt = 3.2-biased_inatt; reflections(9) = -1; end;

param = [lambda sigma2_a sigma2_s sigma2_i B alpha rho bias inatt biased_sigma2_s biased_input biased_inatt];
display(param)
%% make bins and some other vars that can be used many times
bo = binning(dx, B);
% initial distribution
[P_init W_i] = make_P_init3musc(bo, sigma2_i, inatt, biased_inatt);
% the Markov matrix (and its derivatives) used in the absence of clicks
W0 = make_each_F35musc(bo, sigma2_a, 0, 0,0, lambda, ...
                 dt, 'net_input', 0, 'tot_input', 0, 'nclicks', 0);       

% Seems that we could pass sigma2_s and biased_sigma2_s as 0 for this matrix.  
% might make is a tinsy bit faster.



%% now iterate through all trials
ntrials = numel(rawdata);
likey = zeros(ntrials,1);
loglikey = zeros(ntrials,1);
if matlabpool('size')>0 
parfor i = 1:ntrials,
    [loglikey(i) likey(i) out(i)] = single_trial35musc( ...
                param, rawdata(i), ...
                'dx', dx, 'dt', dt, ...
                'total_rate', total_rate, ...
                'var_mode', var_mode, ...
                'adapt_mode', adapt_mode, ...
                'kappa', kappa, 'kappa_timerange', kappa_timerange, ...
                'P_init', P_init, 'W_i', W_i, ...
                'W0', W0,'gain_first',gain_first);
end
else
for i = 1:ntrials,
    [loglikey(i) likey(i) out(i)] = single_trial35musc( ...
                param, rawdata(i), ...
                'dx', dx, 'dt', dt, ...
                'total_rate', total_rate, ...
                'var_mode', var_mode, ...
                'adapt_mode', adapt_mode, ...
                'kappa', kappa, 'kappa_timerange', kappa_timerange, ...
                'P_init', P_init, 'W_i', W_i, ...
                'W0', W0,'gain_first',gain_first);
end
end

%% consolidate log likelihood across all trials
LL  = sum(loglikey);

% if a reflection had taken place, reflect the corresponding derivative

if for_fmin,
    LL  = -LL;
end;

if track_history,
    history3.x    = [history3.x    ; param  ];
    history3.fval = [history3.fval ; LL     ];
%    history3.g    = [history3.g    ; dLL    ];
end;
% and that's it!
