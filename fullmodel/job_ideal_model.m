function job_ideal_model(varargin)

ratname = varargin{1};
ntrials = 80000;
offset  = 1;

do_param = [1 1 1 0 1 1 1 1 1];

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
    fprintf('job_ideal_model.m: could not find fits for rat %s, aborting.\n', ratname);
    return;
else
    fprintf('\nFound fits for rat %s, number %i\n', ratname, rati);
end;

fprintf(1, 'run_ideal_model.m: running task %i, %s\n', rati, ratname);


x_bf = all_chrono_rats(rati).x;

param1 = x_bf;
param1([1 2 4 7]) = [0.0001 0.0001 25 0];

fprintf('run_ideal_model.m: evaluating %s at ideal fit parameters:\n', ratname);
display(param1);



[LL LL_grad likey] = ll_model_wrapper(param1, mydata, 'track_history', 0, ...
                                            'for_fmin', 1, 'do_param', do_param, ...
                                            'use_parallel', 0); %#ok<ASGLU>
           
hitfrac = model_predicted_hitfrac(mydata, likey); %#ok<NASGU>

save(sprintf('idealfit_out_%s_off%i_sigmas.mat', ratname, offset), 'ratname', ...
        'x_bf', 'param1', 'LL', 'LL_grad', 'likey', 'hitfrac');

            
fprintf(1, '\run_ideal_model.m: All done!\n');

  