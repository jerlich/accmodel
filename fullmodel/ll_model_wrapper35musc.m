function [LL likey param_full out] = ll_model_wrapper35musc(param, ratdata, varargin)
% function [LL likey param_full] = ll_model_wrapper35(param, ratdata, varargin)

% BWB, Aug 2011
% BWB, June 2012: sigma2_s scales with magnitude of click input


dx = 0; dt = 0;
pairs = { 
   'track_history'           0  ; ...  
   'do_param'            ones(1,12); ... % 1's for live params, 0's for params to ignore; sum(do_param) should equal numel(rsparam)
   'total_rate'           40    ; ...  % total rate of Poisson events per second
   'default_param'       [0 1 1 0 20 1 0.1 0 0.01 0.1 0.1 0.1]; ... % these are the parameter values defaulted to for entries of 0's in do_param
   'dx'                   0.25  ; ...
   'dt'                   0.02  ; ...
   'var_mode'                1  ; ...  % if 0, sensory variance proportional to nclicks; if 1, sensory variance proportional to total input
   'adapt_mode'              1  ; ...  % if 0, same-side adaptation; if 1, across-side adaptation
   'for_fmin'                1  ; ...  % if 1, then the first 2 outputs are meant to be used for a minimization, so LL and LL_grad are negated to allow maximizing the log likelihood
   'use_parallel'            1  ; ...  % if 1, run parallelized version of code
   'gain_first'              0  ; ...
}; parseargs(varargin, pairs);  

if sum(do_param) ~= numel(param),
    fprintf(2, 'll_model_wrapper35musc: sum(do_param) does not match numel(param)!\n');
    display(do_param);
    display(param);
    return;
end;

ind = find(do_param == 1);
param_full = default_param;
for i = 1:numel(ind),
    param_full(ind(i)) = param(i);
end;

    [LL likey out] = ll_all_trials35musc_parallel(param_full, ratdata, 'dx', dx, 'dt', dt, ...
                                               'total_rate', total_rate, ...
                                               'for_fmin', for_fmin, ...
                                               'var_mode', var_mode, ...
                                               'adapt_mode', adapt_mode, ...
                                               'track_history', track_history, ...
                                               'gain_first',gain_first);



if track_history,
    global history2; % This is where an overall record of all the calls to ll_model_wrapper is kept.

    if isfield(history2, 'x'),
        history2.x    = [history2.x    ; param];
        history2.fval = [history2.fval ; LL];
%        history2.g    = [history2.g    ; LL_grad];
    else
        fprintf(1, '\n\nll_model_wrapper35musc: unable to append current function call to history2\n\n');
    end
end;