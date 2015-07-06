function [LL LL_grad likey param_full] = ll_model_wrapper(param, ratdata, varargin)
% function [LL LL_grad likey param_full] = ll_model_wrapper(param, ratdata, varargin)

% BWB, Aug 2011


dx = 0; dt = 0;
pairs = { 
   'track_history'           0  ; ...  
   'do_param'            ones(1,9); ... % 1's for live params, 0's for params to ignore; sum(do_param) should equal numel(rsparam)
   'default_param'       [0 1 1 0 20 1 0.1 0 0.01]; ... % these are the parameter values defaulted to for entries of 0's in do_param
   'dx'                   0.25  ; ...
   'dt'                   0.02  ; ...
   'for_fmin'                1  ; ...  % if 1, then the first 2 outputs are meant to be used for a minimization, so LL and LL_grad are negated to allow maximizing the log likelihood
   'use_parallel'            1  ; ...  % if 1, run parallelized version of code
}; parseargs(varargin, pairs);  

if sum(do_param) ~= numel(param),
    fprintf(2, 'll_model_wrapper: sum(do_param) does not match numel(param)!\n');
    display(do_param);
    display(param);
    return;
end;

ind = find(do_param == 1);
param_full = default_param;
for i = 1:numel(ind),
    param_full(ind(i)) = param(i);
end;

if use_parallel,
    [LL LL_grad likey] = ll_all_trials3_parallel(param_full, ratdata, 'dx', dx, 'dt', dt, ...
                                               'for_fmin', for_fmin, ...
                                               'track_history', track_history);
else
    [LL LL_grad likey] = ll_all_trials3(param_full, ratdata, 'dx', dx, 'dt', dt, ...
                                       'for_fmin', for_fmin, ...
                                       'track_history', track_history);
end;

LL_grad = LL_grad(do_param==1);

if track_history,
    global history2; % This is where an overall record of all the calls to ll_model_wrapper is kept.

    if isfield(history2, 'x'),
        history2.x    = [history2.x    ; param];
        history2.fval = [history2.fval ; LL];
        history2.g    = [history2.g    ; LL_grad];
    else
        fprintf(1, '\n\nll_model_wrapper: unable to append current function call to history2\n\n');
    end
end;