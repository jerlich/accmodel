function job_fmincon35(varargin)

ratname = varargin{1};
ntrials = str2double(varargin{2});
offset  = str2double(varargin{3});

x1 = str2double(varargin{4});
x2 = str2double(varargin{5});
x3 = str2double(varargin{6});
x4 = str2double(varargin{7});
x5 = str2double(varargin{8});
x6 = str2double(varargin{9});
x7 = str2double(varargin{10});
x8 = str2double(varargin{11});
x9 = str2double(varargin{12});
x_init = [x1 x2 x3 x4 x5 x6 x7 x8 x9];

fprintf(1, 'job_fmincon35.m: x_init = %g %g %g %g %g %g %g %g %g\n', x_init);

if nargin<13
    do_param = [1 1 1 0 1 1 1 1 1];
else
    do_param = varargin{13}
end

filename = dir(['../data/chrono_' ratname '_rawdata.mat']);
if isempty(filename),
    filename = dir([ratname '.mat']);
end;
if isempty(filename),
    fprintf('job_fmincon35.m: Cannot find data file for %s!  Quitting... \n', ratname);
    return;
end;
load(['..' filesep 'data' filesep filename.name]);

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

saveto_filename = sprintf('../data/fmincon_out_%s_off%i_%itrials.mat', ratname, offset, numel(mydata));

fprintf('\n\njob_fmincon35.m: mydata has %i trials at offset %i.\n', numel(mydata), offset);
fprintf('    output will be saved as %s\n', saveto_filename);

%% first fmincon run: start from x_init, go where it will
[history f H_fmincon exitflag] = run_fmincon35(mydata, x_init, ...
                                              'do_param', do_param, ...
                                              'total_rate', total_rate); %#ok<*NASGU,*ASGLU>

global history2 history3;  %#ok<*NUSED>


            
run1.history = history;
run1.history2 = history2;
run1.history3 = history3;
run1.H = H_fmincon;
run1.exitflag = exitflag; 

save(saveto_filename, 'ratname', 'x_init', 'do_param', ...
                'trials', 'history', 'history2', 'history3', ...
                'f', 'H_fmincon', 'exitflag', 'run1');

 
%% do a single functional evaluation at best fit parameters and save likely for each trial
myfval = run1.history.fval(end);
x_bf = run1.history.x(end,:);
        
fprintf('\nBest fit parameters: \n');
fprintf(1, 'job_fmincon.m: x_bf = %g %g %g %g %g %g %g %g %g\n', x_bf);
    
[~, ~, likey] = ll_model_wrapper35(x_bf, mydata, 'track_history', 0, 'for_fmin', 1, ...
                                 'do_param', do_param, 'total_rate', total_rate);

save(saveto_filename, 'ratname', 'x_init', 'do_param', ...
                'trials', 'history', 'history2', 'history3', ...
                'f', 'H_fmincon', 'exitflag', 'run1', ...
                'likey', 'x_bf', 'myfval');                