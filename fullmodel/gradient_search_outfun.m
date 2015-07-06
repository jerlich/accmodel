function stop = gradient_search_outfun(x, optimValues, state)  

global x_names history history2 history3;

stop = false;

switch state
    case 'iter'
        for i=1:numel(x),
            fprintf(1, '-- %s = %.6f  ', x_names{i}, x(i));
        end;
        fprintf(1, '\n');

        history.fval = [history.fval; optimValues.fval];
        history.x = [history.x; x];
        history.g = [history.g; optimValues.gradient'];
        
        jobid = getenv('JOB_ID');
        if ~isempty(jobid),
           try
            save(sprintf('job%s_backup.mat', jobid), 'history', 'history2', 'history3');
           end;
        end;

    case 'interrupt'
        % sometimes fmincon can get into a real funk!
        if length(history.fval) > 30 && std(history.fval(end-30:end))/mean(history.fval(end-30:end))<1e-5,
            stop = true;
        end;
    otherwise
end;
