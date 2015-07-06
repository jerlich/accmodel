% function [hessdata Bnd] = run_quadfit2(x_bf, mydata, ratname, offset)

% Another version in the long line of efforts to compute the Hessian around
% the best-fit point.
%
% We're going back to fitting a quadratic to the gradient of the LL, but
% this time, it's going to be a 3-step process:
%
% 1) estimate the hessian H using steps in the cardinal axes
% 2) using this H, get the eigendirections of the landscape; refit H
% 3) with the revised H, compute steps along each eigendirection that would
%    result in a change in LL of 0.05
% 4) refit H again with new steps along eigendirections
%
% BWB, Jan 2012

format short g;


default_param = [0 1 1 0 15 1 0.1 0 0.01];
default_param(1:3) = x_bf(1:3);
default_param(5:9) = x_bf(4:8);

% hessdata = struct('dims', [], 'dims_names', [], ...
%                   'deltas', [], 'xs', [], 'Ls', [], 'grads', [], ...
%                   'H', [], 'd', [], 'c', [], 'e', []);



x_peak = x_bf([1:3 5:8]); % ignore B for this computation: it's non-gaussian
do_param = [1 1 1 0 0 1 1 1 1];

myfunction = @(x)ll_model_wrapper(x, mydata, ...
                                'track_history', 0, 'for_fmin', 1, ...
                                'do_param', do_param, ...
                                'default_param', default_param);

%% step 1
% estimate H along steps in cardinal axes
stepsize = 0.001;
cardinal_deltas = stepsize*eye(numel(x_peak));

[xs Ls grads] = eval_func_at_deltas(myfunction, x_peak, cardinal_deltas);

[H_cardinal, d, c, e] = quadratic_fit(xs, Ls, grads);



%% step 2
% use H_cardinal to estimate eigendirections of LL

display(H_cardinal);
display(inv(H_cardinal));

[V E] = eig(H_cardinal);

% eigen_deltas = stepsize*V;
% [xs Ls grads] = eval_func_at_deltas(myfunction, x_peak, eigen_deltas);
% [H_eigen, d, c, e] = quadratic_fit(xs, Ls, grads);

%% step 3
% what steps would we take to change LL by 0.05?

% [V E] = eig(H_eigen);

new_deltas = zeros(size(cardinal_deltas,1),1);
H = H_cardinal;
peak = + diag(0.5*xs(:,1)'*H*xs(:,1))' + d'*xs(:,1) + c;
for i = 1:numel(x_peak),
    % I'm sure there's a prettier way of doing this, but I haven't thought
    % of it yet
    mystep = 0;
    mypeak = peak;
    myV = V(:,i);
    if myV(2) < 0, myV = -myV; end;
    while abs(mypeak - peak) < 0.10,
        mystep = mystep + 0.0001;
        newx = x_peak' + myV*mystep;
        mypeak = diag(0.5*newx'*H*newx)' + d'*newx + c;
        
    end;
    fprintf('dimension %i: step = %f, mypeak = %f\n', i, mystep, mypeak);
    new_deltas(:,i) = myV*mystep;
end;

%%
new_eigen_deltas = [new_deltas 10*new_deltas];
[xs2 Ls2 grads2] = eval_func_at_deltas(myfunction, x_peak, new_eigen_deltas);
[H_new_eigen, d, c, e] = quadratic_fit([xs xs2(:,2:end)], [Ls Ls2(:,2:end)], [grads grads2(:,2:end)]);


%%
hessdata.H_cardinal = H_cardinal;
hessdata.H_eigen = H_eigen;
hessdata.H = H_new_eigen;

fprintf('\n\n\nDone computing quadfit:\n');
display(hessdata.H);
         
%% do a line scan over B

Bnd = 0;

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