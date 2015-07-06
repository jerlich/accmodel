function [hessdata Bnd] = run_quadfit3(x_bf, mydata, ratname, offset)
% function [hessdata Bnd] = run_quadfit3(x_bf, mydata, [ratname], [offset])

format short g;

if nargin<3,
    ratname = [];
    offset = 1;
end;


default_param = [0 1 1 0 20 1 0.1 0 0.01];
default_param(1:3) = x_bf(1:3);
default_param(5:9) = x_bf(4:8);
do_param = [1 1 1 0 0 1 1 1 1];

x_peak = x_bf([1:3 5:8]); x_peak = x_peak(:); % ignore B
N = 7; % number of dimensions
delta_scales = [0.1 0.5 0.5 0.05 0.005 0.02 0.004];
at_ori = 1;


              
%% evaluate this 7-dimensional hessian:

% first, let's decide what steps we're going to take
deltas = eye(N, N);
for i = 1:N-1,
    for j = i+1:N,
        step = zeros(N,1);
        step([i j]) = 1;
        deltas = [deltas step]; %#ok<AGROW>
    end;
end;
deltas = deltas .* repmat(delta_scales', 1, size(deltas,2));
              
        
fprintf('\n\nNow computing quadfits, first evaluating in cardinal coordinates\n');

[xs Ls grads] = eval_func_at_deltas(@(x)ll_model_wrapper(x, mydata, ...
                                          'track_history', 0, 'for_fmin', 1, ...
                                          'do_param', do_param, ...
                                          'default_param', default_param), ...
                                    x_peak, deltas);

if at_ori,                                
    % normalize so that x_peak is at the origin                                
    xs = xs - repmat(xs(:,1), 1, size(xs,2));
    Ls = Ls - (Ls(1));
end;

[H d c e] = quadfitN(xs, Ls, at_ori);


fprintf('\n\nDone first pass in cardinal coordinates:\n');

display(H);

hessdata.deltas = deltas;
hessdata.H1 = H;
hessdata.d1 = d;
hessdata.c1 = c;
hessdata.e1 = e;
hessdata.xs1 = xs;
hessdata.Ls1 = Ls;
   
%% do a second pass by computing the eigendirections

[V E] = eig(H);

% new_deltas = (deltas'*V)';
new_deltas = [V/10 V/100 V/1000];

display(new_deltas);

[xs2 Ls2 grads2] = eval_func_at_deltas(@(x)ll_model_wrapper(x, mydata, ...
                                          'track_history', 0, 'for_fmin', 1, ...
                                          'do_param', do_param, ...
                                          'default_param', default_param), ...
                                    x_peak, new_deltas);

if at_ori,                                
    % normalize so that x_peak is at the origin                                
    xs2 = xs2 - repmat(xs2(:,1), 1, size(xs2,2));
    Ls2 = Ls2 - (Ls2(1));
end;

[H d c e] = quadfitN([xs xs2], [Ls Ls2], at_ori);

fprintf('\n\nDone second pass in eigen coordinates:\n');

display(H);

hessdata.H = H;
hessdata.H2 = H;
hessdata.d2 = d;
hessdata.c2 = c;
hessdata.e2 = e;
hessdata.xs2 = xs2;
hessdata.Ls2 = Ls2;

save(sprintf('quadfit_out_%s_off%i.mat', ratname, offset), 'ratname', ...
        'x_bf', 'hessdata');
    
    
    
fprintf('\n\nFINAL HESSIANS:\n');

display(inv(hessdata.H1));
display(inv(hessdata.H2));

%% do a line scan over B
% B = 2:32;
% LL_B = zeros(size(B));
% LLgrad_B = zeros(numel(B), numel(x_bf));
% do_param_linescan = [1 1 1 0 1 1 1 1 1];
% fprintf(1, '\njob_quadfit.m: Doing linescan over B\n');
% for i = 1:numel(B),
%     xx = x_bf;
%     xx(4) = B(i);
%     [LL_B(i) LLgrad_B(i,:)] = ll_model_wrapper(xx, mydata, 'for_fmin', 1, 'do_param', do_param_linescan);
% end;
% 
% Bnd.B = B;
% Bnd.LL_B = LL_B;
% Bnd.LLgrad_B = LLgrad_B;
%     
% save(sprintf('quadfit_out_%s_off%i.mat', ratname, offset), 'ratname', ...
%         'x_bf', 'hessdata', 'Bnd');