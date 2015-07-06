function [H, fval, grad, dfval, dgrad] = run_hessian(data, x, varargin)

pairs = { 
    'dx'                   0.25  ; ...
    'dt'                   0.02  ; ...
    'do_param'             ones(1,9); ...
    'delta'                0.001; ...
}; parseargs(varargin, pairs); 


if numel(x) ~= sum(do_param),
    fprintf(2, 'run_hessian.m: numel(x) does not equal sum(do_param), cannot run!');
    H = [];
    return;
end;

x_names = {'lambda' 'sigma2_a'  'sigma2_s'  'sigma2_i' 'B'      'alpha'     'rho'   'bias'  'inatt'};
x_names = x_names(do_param==1);

fprintf(1, 'run_hessian.m: x = ');
fprintf_param_in_line(x, x_names);
fprintf('\n\n');

tic;
[H, fval, grad, dfval, dgrad] = calc_hessian2(@(x)ll_model_wrapper(x, data, 'dx', dx, 'dt', dt, ...
                                                                   'for_fmin', 1, 'do_param', do_param), ...
                                              x, delta);
toc

fprintf(1, 'run_hessian.m: done!\n');

