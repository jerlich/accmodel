function job_fmincon(ratname, x_init, do_param, varargin)

ntrials = +inf;

x_init = [x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12];


fprintf(1, 'job_fmincon.m: x_init = %g %g %g %g %g %g %g %g %g %g %g %g\n', x_init);

filename = dir(['chrono_' ratname '_rawdata.mat']);  % make the infusion data so that 
if isempty(filename),e
    filename = dir([ratname '.mat']);
end;
if isempty(filename),
    fprintf('job_fmincon.m: Cannot find data file for %s!  Quitting... \n', ratname);
    return;
end;
load(filename.name);


mydata = rawdata;

clear rawdata;


%% first fmincon run: start from x_init, go where it will
[history f H_fmincon exitflag] = run_fmincon(mydata, x_init, 'do_param', do_param,'default_param',x_init,'no_hessian',1); %#ok<*NASGU,*ASGLU>  % ADD no hessian flag

global history2 history3;  %#ok<*NUSED>


            
run1.history = history;
run1.history2 = history2;
run1.history3 = history3;
run1.H = H_fmincon;
run1.exitflag = exitflag; 

save(sprintf('fmincon_out_%s_off%i.mat', ratname, offset), 'ratname', 'x_init', 'do_param', ...
                'trials', 'history', 'history2', 'history3', ...
                'f', 'H_fmincon', 'exitflag', 'run1');

%% second fmincon run:
% stuff all variance in sigma2_s, use that as new seed for another run
% thats constrained to search within the sigma2_s and sigma2_a dimension
%x_init2 = history.x(end,:);
%
%if x_init2(2) > 0.1,
%    x_init2(2:3) = [0.001 x_init2(2)+x_in
 
%% finally, do a single run at best fit parameters and save likely for each trial
myfval = run1.history.fval(end);
x_bf = run1.history.x(end,:);
if ~isempty(run2),
    if run2.history.fval(end) < myfval,
        myfval = run2.history.fval(end);
        x_bf = run2.history.x(end,:);
    end;
    
    if ~isempty(run3),
        if run3.history.fval(end) < myfval,
            myfval = run3.history.fval(end);
            x_bf(2) = x_init3(2);
            x_bf(3) = run3.history.x;
        end;
    end;
end;
        
fprintf('\nBest fit parameters: \n');
fprintf(1, 'job_fmincon.m: x_bf = %g %g %g %g %g %g %g %g %g\n', x_bf);
    
[~, ~, likey] = ll_model_wrapper(x_bf, mydata, 'track_history', 0, 'for_fmin', 1, ...
                                 'do_param', do_param);

save(sprintf('fmincon_out_%s_off%i.mat', ratname, offset), 'ratname', 'x_init', 'do_param', ...
                'trials', 'history', 'history2', 'history3', ...
                'f', 'H_fmincon', 'exitflag', 'run1', 'run2', 'run3', ...
                'likey', 'x_bf', 'myfval');                

% %% quadfit
% [H hessdata Bnd] = run_quadfit4(x_bf, mydata, ratname, offset);
%              
% save(sprintf('fmincon_out_%s_off%i.mat', ratname, offset), 'ratname', 'x_init', 'do_param', ...
%                 'trials', 'history', 'history2', 'history3', ...
%                 'f', 'H_fmincon', 'exitflag', 'run1', 'run2', 'run3', ...
%                 'likey', 'H', 'hessdata', 'Bnd', ...
%                 'x_bf', 'myfval');                
