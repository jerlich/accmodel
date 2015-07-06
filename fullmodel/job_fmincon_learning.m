function job_fmincon_learning(varargin)

ratname = varargin{1};
fort  = str2double(varargin{2});

x1 = str2double(varargin{3});
x2 = str2double(varargin{4});
x3 = str2double(varargin{5});
x4 = str2double(varargin{6});
x5 = str2double(varargin{7});
x6 = str2double(varargin{8});
x7 = str2double(varargin{9});
x8 = str2double(varargin{10});
x9 = str2double(varargin{11});
x_init = [x1 x2 x3 x4 x5 x6 x7 x8 x9];

fprintf(1, 'job_fmincon_learning.m: x_init = %g %g %g %g %g %g %g %g %g\n', x_init);


do_param = [1 1 1 0 1 1 1 1 1];


load(sprintf('chrono_%s_fortnight_%i_rawdata.mat', ratname, fort));

mydata = rawdata; 
clear rawdata;

fprintf(1, 'rat %s from %s to %s has %i trials\n\n\n', ratname, twoweeks{1}, twoweeks{2}, numel(mydata)); %#ok<USENS>

%% first fmincon run: start from x_init, go where it will
[history f H_fmincon exitflag] = run_fmincon(mydata, x_init, 'do_param', do_param); %#ok<*NASGU,*ASGLU>

global history2 history3;  %#ok<*NUSED>

save(sprintf('fmincon_out_%s_off%i.mat', ratname, fort), 'ratname', 'x_init', 'do_param', ...
                'history', 'history2', 'history3', ...
                'f', 'H_fmincon', 'exitflag');
            
run1.history = history;
run1.history2 = history2;
run1.history3 = history3;
run1.H = H_fmincon;
run1.exitflag = exitflag; 


%% second fmincon run:
% stuff all variance in sigma2_s, use that as new seed for another run
% that's constrained to search within the sigma2_s and sigma2_a dimension
x_init2 = history.x(end,:);

if x_init2(2) > 0.1,
    x_init2(2:3) = [0.001 x_init2(2)+x_init2(3)-0.001];
    new_defaults = [x_init2(1:3) 0 x_init2(4:end)];
    new_do_param = [1 1 1 0 1 1 1 1 1];
    
    fprintf(1, '\n\nStarting second fmincon run: \n');
    
    [history f H_fmincon exitflag] = run_fmincon(mydata, x_init2, 'do_param', new_do_param, ...
                                                 'default_param', new_defaults); %#ok<*NASGU,*ASGLU>

    save(sprintf('fmincon_out_%s_off%i.mat', ratname, fort), 'ratname', 'x_init', 'do_param', ...
                    'history', 'history2', 'history3', ...
                    'f', 'H_fmincon', 'exitflag', 'run1');
                
    run2.history = history;
    run2.history2 = history2;
    run2.history3 = history3;
    run2.H = H_fmincon;
    run2.exitflag = exitflag;
else
    run2 = [];
    fprintf(1, '\n\nNo need to do second run.\n');
end;

%% third fmincon run:
% just because the second run is sometimes stupid and acts up, we'll do a
% line scan

x_init3 = history.x(end,:);

if x_init3(2) > 0.1,
    x_init3(2:3) = [0.001 x_init3(2)+x_init3(3)-0.001];
    new_defaults = [x_init3(1:3) 0 x_init3(4:end)];
    new_do_param = [0 0 1 0 0 0 0 0 0];
    
    fprintf(1, '\n\nStarting third, restricted fmincon run: \n');
    
    [history f H_fmincon exitflag] = run_fmincon(mydata, x_init3(3), 'do_param', new_do_param, ...
                                                 'default_param', new_defaults); %#ok<*NASGU,*ASGLU>

    save(sprintf('fmincon_out_%s_off%i.mat', ratname, fort), 'ratname', 'x_init', 'do_param', ...
                    'history', 'history2', 'history3', ...
                    'f', 'H_fmincon', 'exitflag', 'run1', 'run2');
                
    run3.history = history;
    run3.history2 = history2;
    run3.history3 = history3;
    run3.H = H_fmincon;
    run3.exitflag = exitflag; 
else
    run3 = [];
    fprintf(1, '\n\nNo Need to do a third run.\n');
end;
 
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
fprintf(1, 'job_fmincon_learning.m: x_bf = %g %g %g %g %g %g %g %g %g\n', x_bf);
    
[~, ~, likey] = ll_model_wrapper(x_bf, mydata, 'track_history', 0, 'for_fmin', 1, ...
                                 'do_param', do_param);

save(sprintf('fmincon_out_%s_off%i.mat', ratname, fort), 'ratname', 'x_init', 'do_param', ...
                'history', 'history2', 'history3', ...
                'f', 'H_fmincon', 'exitflag', 'run1', 'run2', 'run3', ...
                'likey', 'x_bf');                

%% quadfit
[hessdata Bnd] = run_quadfit(x_bf, mydata, ratname, fort);
             
save(sprintf('fmincon_out_%s_off%i.mat', ratname, fort), 'ratname', 'x_init', 'do_param', ...
                'history', 'history2', 'history3', ...
                'f', 'H_fmincon', 'exitflag', 'run1', 'run2', 'run3', ...
                'likey', 'x_bf', 'hessdata', 'Bnd');                
