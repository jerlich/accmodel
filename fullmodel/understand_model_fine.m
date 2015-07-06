function understand_model_fine

%{
     'lambda' 
     'sigma2_a' 
     'sigma2_s' 
     'sigma2_i' 
     'B' 
     'alpha' 
     'rho' 
     'bias' 
     'inatt' 
     'biased_sigma2_s' 
     'biased_input' 
     'biased_inatt'
%}
ax(1)=subplot(4,5,1);
%% data


cfit=load('fmincon_out_cntrl_off0_47580trials');
x_init=[cfit.x_bf cfit.x_bf(3) 1 cfit.x_bf(9)];
load chrono_fof_rawdata.mat

[nR,nL,wr]=foo(rawdata);

% isosum_psycho(nL,nR,wr,'nsumbins',3,'plot_it',true,'ax',ax(1));

%% lapse
if ~exist('lapse_test_fine.mat','file')

figure(11); clf
ax=axes; set(ax,'NextPlot','add');
clear likey out mwr t_param
xind=0;
for x=0.9:0.005:1.1
    params=x_init;
    params(12)=x;
xind=xind+1;
    [~, likey{xind},~,out] = ll_model_wrapper35musc(params, rawdata, 'track_history', 0, 'for_fmin', 1, ...
                                   'total_rate', 40);
                               
     mwr{xind}=collapse_dist(out,params(8));                              
    % plot(ax,nR-nL,mwr{xind},'.','Color',[x x x]);
    % drawnow;
     t_param{xind}=params;
     LL(xind)=sum(log(likey{xind}));                               
end
set(ax,'Color','b')

save lapse_test_fine nR nL wr likey t_param mwr LL

end


%% gain
if ~exist('gain_test_fine.mat','file')

figure(12); clf
ax=axes; set(ax,'NextPlot','add');
clear likey out mwr t_param
xind=0;
for x=[0:0.002:0.05 0.06:0.01:0.2] 
    params=x_init;
    params(11)=x;
    xind=xind+1;
    [~, likey{xind},~,out] = ll_model_wrapper35musc(params, rawdata, 'track_history', 0, 'for_fmin', 1, ...
                                   'total_rate', 40);
                               
     mwr{xind}=collapse_dist(out,params(8));                              
    % plot(ax,nR-nL,mwr{xind},'.','Color',[x x x]);
    % drawnow;
     t_param{xind}=params;
     LL(xind)=sum(log(likey{xind}));                               
                               
end
set(ax,'Color','b')

save gain_test_fine nR nL wr likey t_param mwr

end

%% noise
if ~exist('noise_test_fine.mat','file')
figure(13); clf
ax=axes; set(ax,'NextPlot','add');
xind=0;
clear likey out mwr t_param
for x=[100:100:700 800:10:900 905:5:980 985:1009 1010:0.5:1012 1013:1025 1030:5:1100 1110:10:1200 1300:100:2000]
    params=x_init;
    params(10)=x;
xind=xind+1;
    [~, likey{xind},~,out] = ll_model_wrapper35musc(params, rawdata, 'track_history', 0, 'for_fmin', 1, ...
                                   'total_rate', 40);
                               
     mwr{xind}=collapse_dist(out,params(8));                              
    % plot(ax,nR-nL,mwr{xind},'.','Color',[x x x]);
    % drawnow;
     t_param{xind}=params;
     LL(xind)=sum(log(likey{xind}));                               
                               
end
set(ax,'Color','b')

save noise_test_fine nR nL wr likey t_param mwr
end

%% bias
if ~exist('bias_test_fine.mat','file')
figure(13); clf
ax=axes; set(ax,'NextPlot','add');
xind=0;
clear likey out mwr t_param
for x=[0:0.25:2 2.1:0.1:4 4.5 5]
    params=x_init;
    params(8)=x;
    xind=xind+1;
    [~, likey{xind},~,out] = ll_model_wrapper35musc(params, rawdata, 'track_history', 0, 'for_fmin', 1, ...
                                   'total_rate', 40);
                               
     mwr{xind}=collapse_dist(out,params(8));                              
    % plot(ax,nR-nL,mwr{xind},'.','Color',[x x x]);
    % drawnow;
     t_param{xind}=params;
     LL(xind)=sum(log(likey{xind}));                               
                               
end
set(ax,'Color','b')

save bias_test_fine nR nL wr likey t_param mwr
end

function [R,L,wr]=foo(D)
R=nan(size(D));
L=R; wr=R;

for tx=1:numel(D)
    R(tx)=numel(D(tx).rightbups);
    L(tx)=numel(D(tx).leftbups);
    wr(tx)=D(tx).pokedR;
end

function O=collapse_dist(P,bias)
O=nan(size(P));
   for tx=1:numel(P)
       endP=P(tx).Pf(:,end);
       O(tx)=(sum(endP(P(tx).x>bias)));
   end
