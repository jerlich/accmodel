function job_quadfit35(varargin)

ratname = varargin{1};
ntrials = str2double(varargin{2});
offset  = str2double(varargin{3});

fprintf(1, 'job_quadfit35.m: running\n');

% load data file for this rat
filename = dir(['chrono_' ratname '_rawdata.mat']);
if isempty(filename),
    filename = dir([ratname '.mat']);
end;
if isempty(filename),
    fprintf('job_quadfit35.m: Cannot find data file for %s!  Quitting... \n', ratname);
    return;
end;
load(filename.name);

if ~exist('total_rate', 'var'),
    total_rate = 40;
end;

if numel(rawdata) > ntrials,
    decimation = floor(numel(rawdata)/ntrials);
    trials = offset:decimation:numel(rawdata);
    trials = trials(1:ntrials);
else
    trials = 1:numel(rawdata);
end;
mydata = rawdata(trials);
clear rawdata;

saveto_filename = sprintf('quadfit_out_%s_off%i_%itrials.mat', ratname, offset, numel(mydata));

% load results from fmincon for best fit peak location
best_fits = dir('best_fits_current.mat');
if ~isempty(best_fits),
    load(best_fits(1).name);
    rati = 1;
    while ~strcmp(ratname, all_chrono_rats(rati).ratname) && rati<=numel(all_chrono_rats),
        rati = rati+1;
    end;
    if rati > numel(all_chrono_rats),
        fprintf('job_quadfit35.m: could not find fits for rat %s, aborting.\n', ratname);
        return;
    else
        fprintf('\nFound fits for rat %s, number %i\n', ratname, rati);
    end;
    x_bf = all_chrono_rats(rati).x;
else
    fmincon_filename = sprintf('fmincon_out_%s_off%i_%itrials.mat', ratname, offset, numel(mydata));

    load(fmincon_filename);
    if exist('x_bf', 'var'),
        1;
    elseif exist('run1', 'var'),
        x_bf = run1.history.x(end,:);
    else
        x_bf = history.x(end,:);
    end;
end;
fprintf('job_quadfit35.m: evaluating %s around best fit parameters:\n', ratname);
display(x_bf);

[H hessdata Bnd] = run_quadfit35(x_bf, mydata, ratname, offset, total_rate); %#ok<NASGU,ASGLU>

save(saveto_filename, 'H', 'hessdata', 'Bnd');

fprintf('job_quadfit35.m: %s is DONE!\n', ratname);
