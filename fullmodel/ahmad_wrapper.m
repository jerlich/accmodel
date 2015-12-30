function ahmad_wrapper

do_cntrl=0;
%%

%% Get the data




cntrl_f=dir('fmincon_out_visin*');
if numel(cntrl_f)==0
    job_fmincon35('visin','+inf','1', '0','10','50','0','50','0.3','0.05','0','0')         
    % params = [lambda, sigma2_a, sigma2_s, sigma2_i, B, alpha, rho, bias, inatt]
end


cntrl_f=dir('fmincon_out_visout*');
if numel(cntrl_f)==0
    job_fmincon35('visout','+inf','1', '0','10','50','0','50','0.3','0.05','0','0')         
    % params = [lambda, sigma2_a, sigma2_s, sigma2_i, B, alpha, rho, bias, inatt]
end


cntrl_f=dir('fmincon_out_audin*');
if numel(cntrl_f)==0
    job_fmincon35('audin','+inf','1', '0','10','50','0','50','0.3','0.05','0','0')         
    % params = [lambda, sigma2_a, sigma2_s, sigma2_i, B, alpha, rho, bias, inatt]
end


cntrl_f=dir('fmincon_out_audout*');
if numel(cntrl_f)==0
    job_fmincon35('audout','+inf','1', '0','10','50','0','50','0.3','0.05','0','0')         
    % params = [lambda, sigma2_a, sigma2_s, sigma2_i, B, alpha, rho, bias, inatt]
end





cntrl_f=dir('fmincon_out_wilke*');
if numel(cntrl_f)==0
    job_fmincon35('wilke','100','1', '0.0384','8.2200','53.8019','1','16.0429','0.3455','0.0423','0.0817','0.0835')
    % job_fmincon35musc('wilke','100','1', '0.0384','8.2200','53.8019','1','16.0429','0.3455','0.0423','0.0817','0.0835','1','54','.08')
        
    % params = [lambda, sigma2_a, sigma2_s, sigma2_i, B, alpha, rho, bias, inatt]
end



