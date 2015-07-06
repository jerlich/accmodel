function job_linescan_sigmas(varargin)

taskid = str2double(getenv('SGE_TASK_ID'));



ratname = 'B083'; % he should have 46000 trials
x_bf = [0.71881    0.0045015       40.242   2.3881e-14       14.987       0.1287     0.049237    -0.067342     0.026089];



load(['chrono_' ratname '_rawdata.mat']);
total_rate = 40;
do_param = ones(1,9);

myfunction = @(x)ll_model_wrapper35(x, rawdata, ...
                                'track_history', 0, 'for_fmin', 1, ...
                                'total_rate', total_rate, ...
                                'do_param', do_param, ...
                                'use_parallel', 0);



bf_sigs = x_bf(2:4);

% set up the sig parameters for the linescan:
% line = -20:20;
line = linspace(-0.2, 0.2, 100);
linescan = repmat(line(:), 1, 3) + repmat(bf_sigs, numel(line), 1);

i = taskid;
LL = myfunction([x_bf(1) linescan(i,:) x_bf(5:end)]);
fprintf('LL(%g, %g, %g) = %g\n', linescan(i,1), linescan(i,2), linescan(i,3), LL);

save(sprintf('linescan_sigmas_point%05i.mat', taskid), 'ratname', 'x_bf', 'taskid', 'line', 'linescan', 'LL');


