function job_fmincon35musc(varargin)

ratname = varargin{1};
ntrials = +inf;
offset  = 0;

x_init = varargin{2};

%%% unpack model parameters
%lambda   = param(1);
%sigma2_a = param(2);
%sigma2_s = param(3);
%sigma2_i = param(4);
%B        = param(5);
%alpha    = param(6);
%rho      = param(7);
%bias     = param(8);
%inatt    = param(9);
%biased_sigma2_s = param(10);
%biased_input = param(11);
%biased_inatt = param(12);




if numel(varargin)>=3
    do_param = varargin{3};
else
    do_param=ones(size(x_init));
end

if numel(varargin)>=4
    save_suff = varargin{4};;
else
    save_suff = '';
end


if numel(varargin)>=5
    gain_first = varargin{5};;
else
    gain_first=0;
end


fprintf(1, '%s: x_init = %g %g %g %g %g %g %g %g %g %g %g %g\n',mfilename, x_init);


filename = dir(['../data/chrono_' ratname '_rawdata.mat']);
if isempty(filename),
    filename = dir([ratname '.mat']);
end;
if isempty(filename),
    fprintf('%s: Cannot find data file for %s!  Quitting... \n',mfilename, ratname);
    return;
end;
load(['../data/' filename.name]);

if ~exist('total_rate', 'var'),
    total_rate = 40;
end;

mydata = rawdata;

clear rawdata;

saveto_filename = sprintf('../data/fmincon_out_%s_%s.mat', ratname,save_suff);

fprintf('\n\njob_fmincon35musc.m: mydata has %i trials.\n', numel(mydata));
fprintf('    output will be saved as %s\n', saveto_filename);

%% first fmincon run: start from x_init, go where it will
[history f exitflag] = run_fmincon35musc(mydata, x_init, ...
                                              'do_param', do_param, ...
                                              'total_rate', total_rate,...
                                              'gain_first',gain_first); %#ok<*NASGU,*ASGLU>

global history2 history3;  %#ok<*NUSED>


            
run1.history = history;
run1.history2 = history2;
run1.history3 = history3;
run1.exitflag = exitflag; 

% save(saveto_filename, 'ratname', 'x_init', 'do_param', ...
%                 'history', 'history2', 'history3', ...
%                 'f', 'exitflag', 'run1');

 
%% do a single functional evaluation at best fit parameters and save likely for each trial
myfval = run1.history.fval(end);
x_bf = x_init;
x_bf(do_param==1) = run1.history.x(end);


        
 fprintf('\nBest fit parameters: \n');
 fprintf(1, 'job_fmincon35musc.m: x_bf = %g %g %g %g %g %g %g %g %g %g %g %g\n', x_bf);
    
 [~, likey] = ll_model_wrapper35musc(x_bf, mydata, 'track_history', 0, 'for_fmin', 1, ...
                                   'total_rate', total_rate);

 save(saveto_filename, 'ratname', 'x_init', 'do_param', ...
                 'history', 'history2', 'history3', ...
                 'f',  'exitflag', 'run1', ...
                 'likey', 'x_bf', 'myfval');                
