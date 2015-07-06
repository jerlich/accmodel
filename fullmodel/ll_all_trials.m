function [LL dLL likey loglikey_grad output] = ll_all_trials(param, rawdata, varargin)

global history3

pairs = { ...
    'for_fmin'             0  ; ...
    'dx'                0.25  ; ...
    'dt'                0.02  ; ...
    'show_iter'            1  ; ...
    'track_history'        0  ; ...
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

if sigma2_a < 0, sigma2_a = 0; end;
if sigma2_s < 0, sigma2_s = 0; end;
if sigma2_i < 0, sigma2_i = 0; end;
if alpha < 0, alpha = 0; end;
if rho <0, rho = 0; end;
if inatt < 0, inatt = 0; end; if inatt > 1, inatt = 1; end;

param = [lambda sigma2_a sigma2_s sigma2_i B alpha rho bias inatt];

%% make bins and some other vars that can be used many times
bo = binning(dx, B);
% initial distribution
[P_init W_i] = make_P_init(bo, sigma2_i, inatt);
% the Markov matrix (and its derivatives) used in the absence of clicks
W0 = make_each_F(bo, sigma2_a, sigma2_s, lambda, ...
                 dt, 'net_input', 0, 'nclicks', 0);
W0 = W0{1};             


%% now iterate through all trials
ntrials = numel(rawdata);
likey = zeros(ntrials,1);
loglikey = zeros(ntrials,1);
loglikey_grad = zeros(ntrials, numel(param));

if nargout > 4,
    output = struct('Pf', cell(ntrials,1), ...
                    'Pb', cell(ntrials,1), ...
                    'x',  cell(ntrials,1), ...
                    't',  cell(ntrials,1)  ...
                    );
end;

for i = 1:ntrials,
    if show_iter==1 && rem(i,1000)==0,
        fprintf(1, '       ll_all_trials: running trial %d/%d\n', i, ntrials);
    end;
    
    [loglikey(i) loglikey_grad(i,:) likey(i) out] = single_trial( ...
                                                                param, rawdata(i), ...
                                                                'dx', dx, 'dt', dt, 'bo', bo, ...
                                                                'P_init', P_init, 'W_i', W_i, ...
                                                                'W0', W0);
        
    if nargout > 4,
        output(i).Pf = out.Pf;
        output(i).Pb = out.Pb;
        output(i).x  = out.x;
        output(i).t  = out.t;
    end;
end;

%% consolidate log likelihood across all trials
LL  = sum(loglikey);
dLL = sum(loglikey_grad, 1);

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