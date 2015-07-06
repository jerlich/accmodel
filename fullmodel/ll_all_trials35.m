function [LL dLL likey loglikey_grad] = ll_all_trials35(param, rawdata, varargin)

global history3

pairs = { ...
    'for_fmin'             0  ; ...
    'dx'                0.25  ; ...
    'dt'                0.02  ; ...
    'total_rate'           40 ; ... % total rate of Poisson events per second
    'var_mode'             1  ; ...
    'adapt_mode'           0  ; ... % if 0, same-side adaptation; if 1, across-side adaptation
    'show_iter'            1  ; ...
    'track_history'        0  ; ...
    'kappa'            [1 1]  ; ... % the biased weight of [left right] clicks
    'kappa_timerange'   [0 1000] ; ... % the time range over which kappa values are applied
    };
parseargs(varargin, pairs);


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

if numel(kappa) == 1,
    kappa = [kappa 1];
end;

% used to keep track of reflections so that the gradient can be made to
% point in the right direction
reflections = ones(size(param)); 

if sigma2_a < 0, sigma2_a = abs(sigma2_a); reflections(2) = -1; end;
if sigma2_s < 0, sigma2_s = abs(sigma2_s); reflections(3) = -1; end;
if sigma2_i < 0, sigma2_i = abs(sigma2_i); reflections(4) = -1; end;
if alpha < 0, alpha = abs(alpha); reflections(6) = -1; end;
if rho <0, rho = abs(rho); reflections(7) = -1; end;
if inatt < 0, inatt = abs(inatt); reflections(9) = -1; end; 
if inatt > 1, inatt = inatt-1; reflections(9) = -1; end;


param = [lambda sigma2_a sigma2_s sigma2_i B alpha rho bias inatt];

%% make bins and some other vars that can be used many times
bo = binning(dx, B);
% initial distribution
[P_init W_i] = make_P_init3(bo, sigma2_i, inatt);
% the Markov matrix (and its derivatives) used in the absence of clicks
W0 = make_each_F35(bo, sigma2_a, sigma2_s, lambda, ...
                 dt, 'net_input', 0, 'tot_input', 0, 'nclicks', 0);             


%% now iterate through all trials
ntrials = numel(rawdata);
likey = zeros(ntrials,1);
loglikey = zeros(ntrials,1);
loglikey_grad = zeros(ntrials, numel(param));

for i = 1:ntrials,
    if show_iter==1 && rem(i,1000)==0,
        fprintf(1, '       ll_all_trials35: running trial %d/%d\n', i, ntrials);
    end;

    [loglikey(i) loglikey_grad(i,:) likey(i)] = single_trial35( ...
                                                                param, rawdata(i), ...
                                                                'dx', dx, 'dt', dt, 'bo', bo, ...
                                                                'total_rate', total_rate, ...
                                                                'var_mode', var_mode, ...
                                                                'adapt_mode', adapt_mode, ...
                                                                'kappa', kappa, 'kappa_timerange', kappa_timerange, ...
                                                                'P_init', P_init, 'W_i', W_i, ...
                                                                'W0', W0);
end;



%% consolidate log likelihood across all trials
LL  = sum(loglikey);
dLL = sum(loglikey_grad, 1);

dLL = dLL.*reflections;

if for_fmin,
    LL  = -LL;
    dLL = -dLL;
end;

if track_history,
    history3.x    = [history3.x    ; param  ];
    history3.fval = [history3.fval ; LL     ];
    history3.g    = [history3.g    ; dLL    ];
end;
% and that's it!