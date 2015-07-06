function job_linescan(varargin)

ratname = varargin{1};
ntrials = str2double(varargin{2});
offset  = str2double(varargin{3});

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

do_param = [1 1 1 0 1 1 1 1 1];
x = [???];


% linescan in every parameter 

[plines pfvals pgrads] = run_linescan(mydata, x, 'do_param', do_param);


save(sprintf('linescan_%s_off%i.mat', ratname, offset), 'ratname', 'x', 'do_param', ...
                'trials', 'plines', 'pfvals', 'pgrads');