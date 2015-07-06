function fit_musc_rats

%%
addpath ~/ratter/Analysis/Pbups

%% Get the data

ratnames=bdata('select distinct(ratname) from pbups.pbupssumm')

for rx=1:numel(ratnames)
	this_rat=ratnames{rx}; 
	if ~exist(['chrono_' this_rat '_rawdata.mat'],'file')
    % Use all control sessions, not just the fof ones.

	    first_sess=bdata('select min(sessid) from ratinfo.infusions where ratname="{S}"',this_rat);
	    last_sess=bdata('select max(sessid) from ratinfo.infusions where ratname="{S}"',this_rat);
	    inf_sess=bdata('select sessid from ratinfo.infusions where ratname="{S}"',this_rat);
	    [sessid,n_trials]=bdata('select sessid, n_done_trials from sessions where left_correct>0.70 and right_correct>0.70 and protocol="Pbups" and ratname="{S}" and sessid>"{S}"',this_rat,first_sess);
	    [sesstodo,si]=setdiff(sessid,inf_sess);
	    tot_trials=sum(n_trials(si));
	    fprintf('Packaging %s, %d sessions, %d trials....',this_rat,numel(sesstodo),numel(tot_trials));
	    rawdata=package_pbups_data(this_rat,sesstodo,'gamma_range',[0.001 10]);
		fprintf('done\n');

		tot_trials=numel(rawdata);
	    if tot_trials>30000
	    	rp=randperm(tot_trials);
	    	rawdata=rawdata(rp(1:30000));
	    end

	    save(['chrono_' this_rat '_rawdata.mat'],'rawdata')

	end


	cntrl_f=dir(['fmincon_out_'  this_rat '*']);
	if numel(cntrl_f)==0
		fprintf('Fitting %s\n',this_rat);
	   
	    job_fmincon35(this_rat,'+inf','0', '0.0384','8.2200','53.8019','1','16.0429','0.3455','0.0423','0.0817','0.0835')
	    
	    % params = [lambda, sigma2_a, sigma2_s, sigma2_i, B, alpha, rho, bias, inatt]
	else
		fprintf('Skipping %s, already done\n',this_rat);
	end

end


