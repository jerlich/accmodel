function [LL, varargout] = ll_all_trials_opt(param, rawdata,varargin)


dx = 0.25; dt = 0.02; total_rate = 40; var_mode = 1; adapt_mode = 0; kappa = [1 1]; kappa_timerange = [0 1000];
pairs = { ...
    'for_fmin'             0  ; ...
    'dx'                0.25  ; ...
    'dt'                0.02  ; ...
    'total_rate'           40 ; ... % total rate of Poisson events per second
    'var_mode'             1  ; ...
    'adapt_mode'           1  ; ...
    'track_history'        0  ; ...
    'kappa'            [1 1]  ; ... % the biased weight of [left right] clicks
    'kappa_timerange'   [0 1000] ; ... % the time range over which kappa values are applied
    'insert_as'          '';...
    'save_to'            [];...
    'use_parallel'          true;
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
    biasfit=false;
    biased_sigma2_s = param(3);
    biased_input = 1;
    biased_inatt = param(9);
    
elseif numel(param)==12
    biasfit=true;
    %% unpack model parameters
    
    biased_sigma2_s = param(10);
    biased_input = param(11);
    biased_inatt = param(12);
else
    error('Can only process 9 or 12 parameters')
end

param = [lambda sigma2_a sigma2_s sigma2_i B alpha rho bias inatt biased_sigma2_s biased_input biased_inatt];

if check_params(param)
    %
    LL=log(0.01)*numel(rawdata);
    if nargout>1
        varargout{1}=nan(size(rawdata));
    end
    return;
end



%% make bins and some other vars that can be used many times
bo = binning(dx, B);
% initial distribution
[P_init, W_i] = make_P_init3musc(bo, sigma2_i, inatt, biased_inatt);
% the Markov matrix (and its derivatives) used in the absence of clicks
W0 = make_each_F35musc(bo, sigma2_a, 0, 0,0, lambda, ...
    dt, 'net_input', 0, 'tot_input', 0, 'nclicks', 0);

W = @(x,y,z)(make_each_F35musc(bo, sigma2_a, sigma2_s,biased_sigma2_s, biased_input,lambda, dt, 'net_input', x, 'tot_input', y, ...
    'nclicks', z, 'total_rate', total_rate, 'W0', W0, 'var_mode', var_mode));
%% now iterate through all trials
ntrials = numel(rawdata);

slice = 30;
loglikey = zeros(slice, ceil(ntrials/slice));

sto = @(wr,lb,rb,T)(single_trial_opt( W,  param, wr,lb,rb,T, ...
    'dx', dx, 'dt', dt, ...
    'total_rate', total_rate, ...
    'var_mode', var_mode, ...
    'adapt_mode', adapt_mode, ...
    'kappa', kappa, 'kappa_timerange', kappa_timerange, ...
    'P_init', P_init, 'W_i', W_i, ...
    'W0', W0));

if use_parallel
    parfor i = 1:size(loglikey,2)
        startidx = (i-1)*slice+1;
        endidx = startidx + slice - 1 ;
        loglikey(:,i) = do_trials(rawdata(startidx:min(ntrials,endidx)),slice,sto)  ;    
    end
else
    % sometimes it is more efficient to parallelize "higher up" so we run
    % each set of parameters on 1 core.
    for i = 1:size(loglikey,2)
        startidx = (i-1)*slice+1;
        endidx = startidx + slice - 1 ;
        loglikey(:,i) = do_trials(rawdata(startidx:min(ntrials,endidx)),slice,sto) ;       
    end
end


loglikey = loglikey(:);
loglikey = loglikey(1:ntrials);



%% consolidate log likelihood across all trials
LL  = sum(loglikey);

if ~isempty(save_to)
    % Since writing to the DB takes too much time, we are going to write to a file
    % and then process the file in another thread
    if biasfit
        str = '%s, %.8f, %.8f, %.8f, %.8f, %.8f, %.8f, %.8f, %.8f , %.8f, %.8f, %.8f, %.8f, %.8f \n';
        fprintf(save_to,str,insert_as,LL,lambda, sigma2_a, sigma2_s, sigma2_i, B, alpha, rho, bias, inatt ,biased_sigma2_s,biased_input,biased_inatt);
    else
        str = '%s, %.8f, %.8f, %.8f, %.8f, %.8f, %.8f, %.8f, %.8f , %.8f, %.8f \n';
        fprintf(save_to,str,insert_as,LL,lambda, sigma2_a, sigma2_s, sigma2_i, B,alpha, rho, bias, inatt);
  
        
    end
    insert_as=[];
end

if ~isempty(insert_as)
    try
        if biasfit
            sql1 = 'insert into pbups.unifits ';
            sql2 = '(dataset, LL, lambda, sigma2_a, sigma2_s, sigma2_i, alpha, rho, bias, inatt, B,biased_sigma2_s,biased_input,biased_inatt) values ';
            sql3 = '("{S}","{S}","{S}","{S}","{S}","{S}","{S}","{S}","{S}","{S}","{S}","{S}","{S}","{S}")';
            bdata([sql1 sql2 sql3],insert_as,LL,lambda, sigma2_a, sigma2_s, sigma2_i, alpha, rho, bias, inatt, B,biased_sigma2_s,biased_input,biased_inatt);
            
        else
            sql1 = 'insert into pbups.fits ';
            sql2 = '(dataset, LL, lambda, sigma2_a, sigma2_s, sigma2_i, alpha, rho, bias, inatt, B) values ';
            sql3 = '("{S}","{S}","{S}","{S}","{S}","{S}","{S}","{S}","{S}","{S}","{S}")';
            bdata([sql1 sql2 sql3],insert_as,LL,lambda, sigma2_a, sigma2_s, sigma2_i, alpha, rho, bias, inatt, B);
        end
    catch me
        showerror(me);
    end
end


if nargout>1
    varargout{1}=loglikey;
end


function bad=check_params(p)
bad = p(2) < 0 || p(3) < 0 || p(4) < 0 || p(10) < 0 || p(6)< 0 || p(7) < 0 || p(9)<0 || p(9)>2 || p(11)<0 || p(11)>2;

function ll=do_trials(rawdata,slice,sto)

ll = nan(slice,1);

for tx=1:numel(rawdata)
    T = rawdata(tx).T;
    rb = rawdata(tx).rightbups;
    lb = rawdata(tx).leftbups;
    rb = rb(rb<T);
    lb = lb(lb<T);
    ll(tx) = sto(rawdata(tx).pokedR,lb,rb,T);
end



