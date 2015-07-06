function job_quadfit2(varargin)

ratname = varargin{1};
ntrials = str2double(varargin{2});
offset  = str2double(varargin{3});

fprintf(1, 'job_quadfit2.m: running\n');

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
best_fits = dir('best_fits_all_trials.mat');
if ~isempty(best_fits),
    load(best_fits(1).name);
    rati = 1;
    while ~strcmp(ratname, all_chrono_rats(rati).ratname) && rati<=numel(all_chrono_rats),
        rati = rati+1;
    end;
    if rati > numel(all_chrono_rats),
        fprintf('job_quadfit2.m: could not find fits for rat %s, aborting.\n', ratname);
        return;
    else
        fprintf('\nFound fits for rat %s, number %i\n', ratname, rati);
    end;
    x_bf = all_chrono_rats(rati).x;
else
    load(sprintf('fmincon_out_%s_off%i.mat', ratname, offset));
    if exist('x_bf', 'var'),
        x_bf = x_bf;
    else
        x_bf = run1.history.x(end,:);
    end;
end;
fprintf('job_quadfit2.m: evaluating %s around best fit parameters:\n', ratname);
display(x_bf);

[H hessdata Bnd] = run_quadfit4(x_bf, mydata, ratname, offset);

fprintf('job_quadfit2.m: DONE!\n', ratname);
