function job_quadfit(varargin)

ratname = varargin{1};
ntrials = str2double(varargin{2});
offset  = str2double(varargin{3});

fprintf(1, 'job_quadfit.m: running\n');

do_param = [1 1 1 0 1 1 1 1 1];
default_param = [0 1 1 0 20 1 0.1 0 0.01];

% load data file for this rat
load(['chrono_' ratname '_rawdata.mat']);
if numel(rawdata) > ntrials,
    decimation = floor(numel(rawdata)/ntrials);
    trials = offset:decimation:numel(rawdata);
    trials = trials(1:ntrials);
else
    trials = 1:numel(rawdata);
end;
mydata = rawdata(trials);
clear rawdata;

% load results from fmincon for best fit peak location
load('best_fits_all_trials.mat');
rati = 1;
while ~strcmp(ratname, all_chrono_rats(rati).ratname) && rati<=numel(all_chrono_rats),
    rati = rati+1;
end;
if rati > numel(all_chrono_rats),
    fprintf('job_quadfit.m: could not find fits for rat %s, aborting.\n', ratname);
    return;
else
    fprintf('\nFound fits for rat %s, number %i\n', ratname, rati);
end;
x_bf = all_chrono_rats(rati).x;
fval_bf = all_chrono_rats(rati).fval;
default_param(1:3) = x_bf(1:3);
default_param(5:9) = x_bf(4:8);

fprintf('job_quadfit.m: evaluating %s around best fit parameters:\n', ratname);
display(x_bf);



hessdata = struct('dims', [], 'dims_names', [], ...
                  'deltas', [], 'xs', [], 'Ls', [], 'grads', [], ...
                  'H', [], 'd', [], 'c', [], 'e', []);

              
%% evaluate hessian at x_bf in slices:
% =============================================================
% first, 1D slice in lambda
fprintf('\n\n\nNow computing quadfit in lambda\n');

dims = 1;
dims_names = {'lambda'};
deltas = 0.2*[-2 -1 1 2];
do_param = zeros(1,9); do_param(dims) = 1;
[xs Ls grads] = eval_func_at_deltas(@(x)ll_model_wrapper(x, mydata, ...
                                        'track_history', 0, 'for_fmin', 1, ...
                                        'do_param', do_param, ...
                                        'default_param', default_param), ...
                                        x_bf(dims), deltas);
xs = [x_bf(dims) xs]; xs = xs - xs(1);
Ls = [fval_bf Ls]; Ls = Ls - min(Ls(:));

H = pinv(xs(:).^2)*Ls(:);
d = 0;
c = 0; 
try
    e = Ls(:) - xs(:).^2*H;
catch
    e = NaN;
end;


for param = fieldnames(hessdata)',
    hessdata(1).(param{1}) = eval(param{1});
end;

fprintf('\n\n\nDone computing quadfit in lambda, H = %g and std = %g\n', H, 1/H);

% =============================================================
%% second, sigmas
fprintf('\n\n\nNow computing quadfit in sigma2_a and sigma2_s\n');

dims = 2:3;
dims_names = {'sigma2_a', 'sigma2_s'};
deltas = [0.2 0.5 0  0 0.2 0.2; ...
          0   0   1 -1 1   -1 ]; 
do_param = zeros(1,9); do_param(dims) = 1;
[xs Ls grads] = eval_func_at_deltas(@(x)ll_model_wrapper(x, mydata, ...
                                        'track_history', 0, 'for_fmin', 1, ...
                                        'do_param', do_param, ...
                                        'default_param', default_param), ...
                                        x_bf(dims), deltas); 
xs(1,:) = xs(1,:) - xs(1,1);
xs(2,:) = xs(2,:) - xs(2,1);
Ls = Ls - min(Ls(:));      

[H d c e] = quadfit2(xs, Ls, 1);

for param = fieldnames(hessdata)',
    hessdata(2).(param{1}) = eval(param{1});
end;

save(sprintf('quadfit_out_%s_off%i.mat', ratname, offset), 'ratname', ...
        'x_bf', 'fval_bf', 'hessdata');

fprintf('\n\n\nDone computing quadfit in sigmas\n');
display(H);
display(inv(H));


% =============================================================
%% third, alpharho
fprintf('\n\n\nNow computing quadfit in alpha and rho\n');

dims = 5:6;
dims_names = {'alpha', 'rho'};
deltas = [0.004 -0.004 0      0      0.004  0.004  -0.004   -0.004; ...
          0      0     0.001 -0.001  0.001 -0.001   0.001   -0.001]; 
do_param = zeros(1,9); do_param(dims+1) = 1;

display(do_param);
display(default_param);
display(deltas);

[xs Ls grads] = eval_func_at_deltas(@(x)ll_model_wrapper(x, mydata, ...
                                        'track_history', 0, 'for_fmin', 1, ...
                                        'do_param', do_param, ...
                                        'default_param', default_param), ...
                                        x_bf(dims), deltas);
display(xs);
display(Ls);

xs(1,:) = xs(1,:) - xs(1,1);
xs(2,:) = xs(2,:) - xs(2,1);
Ls = Ls - min(Ls(:));      

[H d c e] = quadfit2(xs, Ls, 1);

display(H)
display(inv(H))

% now let's find the eigen directions of this distribution and resample
% along those
[V E] = eig(H);
new_deltas = (deltas'*V)';

display(new_deltas);

[xs Ls grads] = eval_func_at_deltas(@(x)ll_model_wrapper(x, mydata, ...
                                        'track_history', 0, 'for_fmin', 1, ...
                                        'do_param', do_param, ...
                                        'default_param', default_param), ...
                                        x_bf(dims), new_deltas);
display(xs);
display(Ls);

xs(1,:) = xs(1,:) - xs(1,1);
xs(2,:) = xs(2,:) - xs(2,1);
Ls = Ls - min(Ls(:));      

[H d c e] = quadfit2(xs, Ls, 1);

display(H);
display(inv(H));

for param = fieldnames(hessdata)',
    hessdata(3).(param{1}) = eval(param{1});
end;

save(sprintf('quadfit_out_%s_off%i.mat', ratname, offset), 'ratname', ...
        'x_bf', 'fval_bf', 'hessdata', 'default_param', 'do_param');

fprintf('\n\n\nDone computing quadfit in alpha and rho\n');



% =============================================================
%% fourth, biasbo
fprintf('\n\n\nNow computing quadfit in bias and blowoff\n');

dims = 7:8;
dims_names = {'bias', 'bo'};
deltas = [0.02 -0.02 0      0      0.02  -0.02  -0.02; ...
          0     0    0.02   0.04   0.02   0.02   0.04]; 
do_param = zeros(1,9); do_param(dims+1) = 1;
[xs Ls grads] = eval_func_at_deltas(@(x)ll_model_wrapper(x, mydata, ...
                                        'track_history', 0, 'for_fmin', 1, ...
                                        'do_param', do_param, ...
                                        'default_param', default_param), ...
                                        x_bf(dims), deltas);
xs(1,:) = xs(1,:) - xs(1,1);
xs(2,:) = xs(2,:) - xs(2,1);
Ls = Ls - min(Ls(:));      

[H d c e] = quadfit2(xs, Ls, 1);
display(H)
display(inv(H))

% now let's find the eigen directions of this distribution and resample
% along those
[V E] = eig(H);
new_deltas = (deltas'*V)';

display(new_deltas);

[xs Ls grads] = eval_func_at_deltas(@(x)ll_model_wrapper(x, mydata, ...
                                        'track_history', 0, 'for_fmin', 1, ...
                                        'do_param', do_param, ...
                                        'default_param', default_param), ...
                                        x_bf(dims), new_deltas);
display(xs);
display(Ls);

xs(1,:) = xs(1,:) - xs(1,1);
xs(2,:) = xs(2,:) - xs(2,1);
Ls = Ls - min(Ls(:));      

[H d c e] = quadfit2(xs, Ls, 1);

display(H);
display(inv(H));

for param = fieldnames(hessdata)',
    hessdata(4).(param{1}) = eval(param{1});
end;

save(sprintf('quadfit_out_%s_off%i.mat', ratname, offset), 'ratname', ...
        'x_bf', 'fval_bf', 'hessdata');

fprintf('\n\n\nDone computing quadfit in bias and blowoff\n');
         
%% do a line scan over B
B = 2:32;
LL_B = zeros(size(B));
LLgrad_B = zeros(numel(B), numel(x_bf));
do_param_linescan = [1 1 1 0 1 1 1 1 1];
fprintf(1, '\njob_quadfit.m: Doing linescan over B\n');
for i = 1:numel(B),
    xx = x_bf;
    xx(4) = B(i);
    [LL_B(i) LLgrad_B(i,:)] = ll_model_wrapper(xx, mydata, 'for_fmin', 1, 'do_param', do_param_linescan);
end;
    
save(sprintf('quadfit_out_%s_off%i.mat', ratname, offset), 'ratname', ...
        'x_bf', 'fval_bf', 'hessdata', 'B', 'LL_B', 'LLgrad_B');
            
fprintf(1, '\njob_quadfit.m: All done!\n');

  