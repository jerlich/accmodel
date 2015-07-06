function job_grid(varargin)

ratname = varargin{1};
ntrials = str2double(varargin{2});
A_dim = str2double(varargin{3});
B_dim = str2double(varargin{4});
Amin  = str2double(varargin{5});
Astep = str2double(varargin{6});
Amax  = str2double(varargin{7});
Bmin  = str2double(varargin{8});
Bstep = str2double(varargin{9});
Bmax  = str2double(varargin{10});
taskid = str2double(getenv('SGE_TASK_ID'));

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

x_bf = [];
do_param = [1 1 1 0 1 1 1 1 1];

%% construct the grid
A = Amin:Astep:Amax;
B = Bmin:Bstep:Bmax;

AB = [];
for i = 1:numel(A),
    for j = 1:numel(B),
        AB = [AB; A(i) B(j)]; %#ok<AGROW>
    end;
end;

if taskid <= size(AB,1),
    fprintf('job_grid.m: running task %i\n', taskid);

    jobfile = dir(sprintf('gridpoint_%05i.mat', taskid));

    if isempty(jobfile),
        param1 = x_bf;
        param1(A_dim) = AB(taskid, 1);
        param1(B_dim) = AB(taskid, 2);

        display(param1);

        [LL LL_grad likey] = ll_model_wrapper(param1, mydata, 'track_history', 0, ...
                                            'for_fmin', 0, 'do_param', do_param, ...
                                            'use_parallel', 0); %#ok<NASGU,ASGLU>

        save(sprintf('gridpoint_%05i.mat', taskid), 'ratname', 'x_bf', 'param1', 'AB', 'A', 'B', 'taskid', 'trials', 'LL', 'LL_grad', 'likey');

        fprintf('job_grid.m: all done!\n');
    else
        fprintf('job_grid:m task %i has already been run!\n', taskid);
    end;
    
    return;
else
    fprintf('job_grid.m: don''t have that many tasks!\n');
end;
