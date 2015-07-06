function [plines pfvals pgrads] = run_linescan(data, x_c, varargin)

pairs = { 
    'dx'                   0.25  ; ...
    'dt'                   0.02  ; ...
    'do_param'             ones(1,9); ...
    'extent'               0.1   ; ... 
    'ndots'                20    ; ... % number of points to eval in each dimension
}; parseargs(varargin, pairs);   


x_names = {'lambda' 'sigma2_a'  'sigma2_s'  'sigma2_i' 'B'      'alpha'     'rho'   'bias'  'inatt'};
x_names = x_names(do_param==1);

fprintf(1, '\n\nrun_linescan.m: x_c = ');
fprintf_param_in_line(x_c, x_names);
fprintf('\n\n');


plines = zeros(numel(x_c), ndots);
pfvals = zeros(numel(x_c), ndots);
pgrads = zeros(numel(x_c), ndots, numel(x_c));
for p = 1:numel(x_c),
    fprintf(1, 'running %i/%i parameter\n', p, numel(x_c));
    plines(p,:) = linspace(x_c(p)-extent, x_c(p)+extent, ndots);
    
    for i = 1:ndots,
        x = x_c;
        x(p) = plines(p,i);
        [pfvals(p,i) pgrads(p,i,:)] = ll_model_wrapper(x, data, 'dx', dx, 'dt', dt, ...
            'track_history', 0, 'for_fmin', 1, 'do_param', do_param);
    end;
end;



fprintf(1, '\n\nrun_linescan.m: DONE! \n ');
fprintf('\n\n');