function job_all_hessians(varargin)

fprintf(1, 'job_all_hessian.m: computing all hessians in subdirectories of %s \n', pwd);

fldrs = dir('B*');

all_hessians = struct('ratname', [], ...
                      'x', [], ...
                      'H', [], ...
                      'ci', []);

for i = 1:numel(fldrs)
    if fldrs(i).isdir,
        fprintf(1, 'going into directory %s \n', fldrs(i).name);
        
        cd(fldrs(i).name);
        
        datafile = dir('fmincon_out*');
        rawdatafile = dir('chrono*.mat');
        if ~isempty(datafile) && ~isempty(rawdatafile),
            load(rawdatafile(1).name);

            for runs = 1:numel(datafile),
                load(datafile(runs).name);
                
                H = run_hessian(rawdata(trials), history.x(end,:), 'do_param', do_param);
                ci = sqrt(diag(inv(H)));
                
                fprintf('for rat %s: \n', ratname);
                display(H);
                display(ci);
                
                n = numel(all_hessians)+1;
                all_hessians(n).ratname = ratname;
                all_hessians(n).x = history.x(end,:);
                all_hessians(n).H = H;
                all_hessians(n).ci = ci;
            end;
        else
            fprintf(1, '     this directory doesn"t have an fmincon output file!\n');
        end;
        
        save(sprintf('all_hessians.mat'), 'all_hessians');
        
        cd('..');
    end;
end;

save(sprintf('all_hessians.mat'), 'all_hessians');
