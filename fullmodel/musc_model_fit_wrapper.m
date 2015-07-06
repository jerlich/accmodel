function musc_model_fit_wrapper

do_cntrl=0;
%%
addpath ~/ratter/Analysis/Pbups

%% Get the data


if ~exist('chrono_cntrl_rawdata.mat','file')
    % Use all control sessions, not just the fof ones.
    [csessid]=bdata('select distinct(cntrl_sessid) from pbups.pbupssumm');
    [rawdata] = package_pbups_data('cntrl', [csessid],'gamma_range',[0.001 10],'hitfrac_thresh',0.7);
    save chrono_cntrl_rawdata rawdata
end

if ~exist('chrono_fof_rawdata.mat','file')
    [sessid,infside]=bdata('select sessid,sign(dose) from pbups.pbupssumm where expgroup="fof"');
    iside(infside==-1)='l';
    iside(infside==1)='r';
    % for consistency let's try the first 250.  If we want we can try with all trials... will give more data.
    [rawdata] = package_pbups_data('cntrl', sessid, 'ipsi',iside,'gamma_range',[0.001 10],'upto',250,'hitfrac_thresh',0,'n_done_trials_thresh',0);    
    orig_rawdata=rawdata;
    rawdata = flipit(rawdata);
    save chrono_fof_rawdata rawdata orig_rawdata
end


%% Run the control model
cntrl_f=dir('fmincon_out_cntrl*');
if numel(cntrl_f)==0
    job_fmincon35('cntrl','+inf','0', '0.0384','8.2200','53.8019','1','16.0429','0.3455','0.0423','0.0817','0.0835')
    
    % params = [lambda, sigma2_a, sigma2_s, sigma2_i, B, alpha, rho, bias, inatt]
end

%% Run the fof model
fof_f=dir('fmincon_out_fof*');
cfit=load('fmincon_out_cntrl_off0_47580trials');
x_init=[cfit.x_bf cfit.x_bf(3) 1 cfit.x_bf(9)];
%	x_names = {'lambda' 'sigma2_a'  'sigma2_s'  'sigma2_i' 'B'      'alpha'     'rho'   'bias'  'inatt' 'biased_sigma2_s' 'biased_input' 'biased_inatt'};

if ~exist('fmincon_out_fof_lapse.mat','file')
    
    do_param=zeros(1,12);
    do_param(12)=1;
    job_fmincon35musc('fof',x_init,do_param,'lapse');
    if do_cntrl
        job_fmincon35musc('cntrl',x_init,do_param,'lapse');
    end
    
end


if ~exist('fmincon_out_fof_gain.mat','file')
    
    do_param=zeros(1,12);
    do_param(11)=1;
    job_fmincon35musc('fof',x_init,do_param,'gain');
    if do_cntrl
        job_fmincon35musc('cntrl',x_init,do_param,'gain');
    end
end


%if ~exist('fmincon_out_fof_gain_first.mat','file')
if 0
    do_param=zeros(1,12);
    do_param(11)=1;
    job_fmincon35musc('fof',x_init,do_param,'gain_first',1);
end


if ~exist('fmincon_out_fof_noise.mat','file')
    
    do_param=zeros(1,12);
    do_param(10)=1;
    job_fmincon35musc('fof',x_init,do_param,'noise');
    if do_cntrl
        job_fmincon35musc('cntrl',x_init,do_param,'noise');
    end
end



if ~exist('fmincon_out_fof_bias.mat','file')
    do_param=zeros(1,12);
    do_param(8)=1;
    job_fmincon35musc('fof',x_init,do_param,'bias');
    if do_cntrl
        job_fmincon35musc('cntrl',x_init,do_param,'bias');
    end
end
%%

names={'bias' 'noise' 'gain_first' 'gain' 'lapse'};
clear ll param;
for nx=1:numel(names)
    F(nx)=load(['fmincon_out_fof_' names{nx} '.mat']);
    ll(nx)=F(nx).run1.history.fval(end);
    param(nx)=F(nx).run1.history.x(end);
    fprintf('%s = %f, -LL=%.1f\n',names{nx},param(nx),ll(nx));
end















function D=flipit(D)

% The model code assumes delta = R - L vs PokedR
% We want the model to do delta = contra - ipsi vs PokedC
% So for left infusions things are correct.
% For right infusions things are backwards

for tx=1:numel(D)
    
    
    
    if D(tx).ipsi=='r'
        % Then we have to flip it.
        
        oR=D(tx).rightbups;
        oL=D(tx).leftbups;
        
        D(tx).leftbups=oR;
        D(tx).rightbups=oL;
        D(tx).Delta=-D(tx).Delta;
        D(tx).gamma=-D(tx).gamma;
        D(tx).pokedR=~D(tx).pokedR;
    end
end
