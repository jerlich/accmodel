function job_MC(varargin)

ratname = varargin{1};
ntrials = str2double(varargin{2});
taskid = str2double(getenv('SGE_TASK_ID'));

saveto = 'job_MC_task_%05i.mat';

% for now, ignore bias, alpha, rho, and bo
do_param = [1 1 1 0 1 1 1 0 1];

load(['chrono_' ratname '_rawdata.mat']);

trials = randperm(numel(rawdata));
trials = sort(trials(1:ntrials));
mydata = rawdata(trials);

clear rawdata;


%%
default_params = [0 10 0 0 10 1 0 0 0.02];
lower_bound = run_fmincon('get_lower_bound');
upper_bound = run_fmincon('get_upper_bound');
lower_bound(do_param==0) = default_params(do_param==0);
upper_bound(do_param==0) = default_params(do_param==0);

rand('twister',bitand(round(now*1E10),2^16-1));
param = rand(1,9).*(upper_bound-lower_bound) + lower_bound;

fprintf(1, 'job_MC: trying param = %g %g %g %g %g %g %g %g %g\n', param);


%%
[LL LL_grad likey] = ll_all_trials3(param, mydata, 'show_iter', 0); %#ok<ASGLU>

save(sprintf(saveto, taskid), 'ratname', 'param', 'taskid', 'trials', 'LL', 'LL_grad', 'likey');

fprintf(1, 'job_MC_grid: LL = %g\n', LL);