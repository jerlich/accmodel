% [history f hessian fmincon_exit] = run_fmincon(data, x_init, {'dx', 0.25} ...
%          {'dt', 0.02}, ...
%          {'lower_bound', [-5     0       0       0       4      0.1     0.005   -5      0])
%          {'upper_bound', [+5     200     200     30      32      1.2     0.5     +5      1]; ...
%          {'default_param', [0 1 1 0 20 1 0.1 0 0.01]}, {'do_param', ones(1,9)}, ...
%          {'no_hessian', 0})
%

function [history f hessian fmincon_exit] = run_fmincon(data, x_init, varargin)
    
pairs = { 
    'dx'                   0.25  ; ...
    'dt'                   0.02  ; ...
    'lower_bound'           [-5     0       0       0       4      0.1     0.005   -5      0]; ...
    'upper_bound'           [+5     200     200     30      32      1.2     0.5     +5      1]; ...
    'default_param'       [0 1 1 0 20 1 0.1 0 0.01]; ... % these are the parameter values defaulted to for entries of 0's in do_param 
    'do_param'             ones(1,9); ...
    'no_hessian'           0 ; ...
}; parseargs(varargin, pairs);   

if ischar(data),
   switch data
      case 'get_lower_bound',
         history = lower_bound;
         return;
      case 'get_upper_bound',
         history = upper_bound;
         return;
      otherwise
         error('run_fmincon does not recognize call ''%s''', data);
   end;
end;

% declare global variables to keep track of fmincon's progress
global x_names history history2 history3; %#ok<REDEF>
x_names = {'lambda' 'sigma2_a'  'sigma2_s'  'sigma2_i' 'B'      'alpha'     'rho'   'bias'  'inatt'};

% trim params according to do_param
x_names = x_names(do_param==1);
lower_bound = lower_bound(do_param==1);
upper_bound = upper_bound(do_param==1);
if numel(x_init) > sum(do_param),
    x_init = x_init(do_param==1);
end;

history.x  = []; history.fval  = []; history.g  = []; %#ok<STRNU>
history2.x = []; history2.fval = []; history2.g = []; 
history3.x = []; history3.fval = []; history3.g = []; 


fprintf(1, 'run_fmincon.m: x_init = ');
fprintf_param_in_line(x_init, x_names);
fprintf('\n\n');

tic;
if no_hessian,
   [x_fmincon, f, exitflag, output, ~, grad] = ...
      fmincon(@(x)ll_model_wrapper(x, data, 'dx', dx, 'dt', dt, 'track_history', 1, ...
                                        'for_fmin', 1, 'do_param', do_param, ...
                                        'default_param', default_param), ...
            x_init, [], [], [], [], lower_bound, upper_bound, [], ...
            optimset('Display', 'iter-detailed', ...
                     'DiffMinChange', 0.0001, ...
                     'MaxIter', 200, 'GradObj', 'on', ...
                     'Algorithm', 'interior-point', ...
                     'OutputFcn', 'gradient_search_outfun')...
    );
 hessian = [];
else
   [x_fmincon, f, exitflag, output, ~, grad, hessian] = ...
      fmincon(@(x)ll_model_wrapper(x, data, 'dx', dx, 'dt', dt, 'track_history', 1, ...
                                        'for_fmin', 1, 'do_param', do_param, ...
                                        'default_param', default_param), ...
            x_init, [], [], [], [], lower_bound, upper_bound, [], ...
            optimset('Display', 'iter-detailed', ...
                     'DiffMinChange', 0.0001, ...
                     'MaxIter', 200, 'GradObj', 'on', ...
                     'Algorithm', 'interior-point', ...
                     'OutputFcn', 'gradient_search_outfun')...
    );
end;
toc

fprintf(1, '\n\n\n     run_fmincon: DONE! \n x_fmincon = ');
fprintf_param_in_line(x_fmincon, x_names);
fprintf('\n\n');

if nargout > 3,
    fmincon_exit.exitflag = exitflag;
    fmincon_exit.output = output;
    fmincon_exit.grad = grad;
end;

% ~~~~~~~~~~~~~~~~~~~~~~~~

