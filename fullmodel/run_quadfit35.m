function [H hessdata Bnd] = run_quadfit35(x_bf, mydata, ratname, offset, total_rate)
% function [hessdata Bnd] = run_quadfit35(x_bf, mydata, [ratname], [offset], [total_rate])

if nargin < 3,
    ratname = [];
    offset = 1;
end;

if nargin < 5,
    total_rate = 40;
end;

default_param = x_bf;
do_param = [1 1 1 1 0 0 0 0 1];

x_peak = x_bf(do_param==1); x_peak = x_peak(:); % ignore B
N = sum(do_param); % number of dimensions
cardinal_steps = [0.1 0.01 3 0.1 .1 0.01 0.002 0.03 0.1]; % has as many dimensions as x_peak
cardinal_steps = cardinal_steps(do_param==1)
at_ori = 1;


if isempty(gcp('nocreate'))
    use_parallel=0;
else
    use_parallel=1;
end


myfunction = @(x)ll_model_wrapper35(x, mydata, ...
    'track_history', 0, 'for_fmin', 1, ...
    'total_rate', total_rate, ...
    'do_param', do_param, ...
    'default_param', default_param,...
    'use_parallel',use_parallel);


%%
hessdata = struct('deltas', [], 'xs', [], 'Ls', [], 'grads', [], ...
  'H', [], 'd', [], 'c', [], 'e', [], ...
  'eV', [], 'eE', [], ...
  'bestV', [], 'bestE', [], 'new_basis', []);



%%
deltas = make_deltas(N);
deltas = [deltas 2*deltas];
deltas = deltas .* repmat(cardinal_steps', 1, size(deltas,2));

fprintf('\n\n');
fprintf('run_quadfit35.m: Now computing quadfits, first evaluating in cardinal coordinates, n = %i dimensions\n', N);

[xs Ls grads] = eval_func_at_deltas(myfunction, x_peak, deltas);

xs = xs - repmat(xs(:,1), 1, size(xs,2));
Ls = Ls - Ls(1);

% fit a quadratic function to points Ls evaluated at xs
% and estimate the eigenvectors and eigenvalues of this Hessian matrix
[H d c e] = quadfitN(xs, Ls, at_ori);
[eV eE] = eig(H); 
eE = diag(eE);

bestV = [];
bestE = [];
new_basis{1} = eye(N);

fprintf('\n\n');
fprintf('run_quadfit35.m: Eigenvalues of estimate of H based on cardinal steps only:\n');
display(sort(eE));
fprintf('H = \n');
display(H);
fprintf('inv(H) = \n');
display(inv(H));

for param = fieldnames(hessdata)',
    hessdata(1).(param{1}) = eval(param{1});
end;

%% iterate
n = N-1;
while n>=1 && (max(abs(eE))/min(abs(eE))>100 || any(eE<0)),
    [eigen ind] = max(eE); % find the largest eigenvalue and its index
    newV = eV(:, ind);
    for k = numel(new_basis):-1:2,
        newV = new_basis{k}*newV;
    end;
    
    bestV = [bestV newV]; % and add it to the list of well-estimated eigenvectors
    bestE = [bestE; eigen];
    new_basis{numel(new_basis)+1} = eV(:, 1:n+1~=ind); % use the rest as the new basis vector
    
    deltas = make_deltas(n);
    deltas = [0.5*deltas 1.0*deltas];
    
    % estimate the appropriate step size in every eigendirection in our new
    % basis vector space
    eE = eE(1:n+1~=ind);
    step = ones(n,1);
    for i = 1:n,
        if eE(i) < 10, p = 0.1;
        else           p = eE(i); end;
            step(i) = sqrt(2*0.1/p); 
        end;
        deltas = deltas .* repmat(step, 1, size(deltas,2));

        NDdeltas = deltas;
        for k = numel(new_basis):-1:2,
            NDdeltas = new_basis{k}*NDdeltas;
        end;

        [xs Ls grads] = eval_func_at_deltas(myfunction, x_peak, NDdeltas);
        xs = xs - repmat(xs(:,1), 1, size(xs,2));
        Ls = Ls - Ls(1);

        [H d c e] = quadfitN([zeros(size(deltas,1),1) deltas], Ls, at_ori);
        [eV eE] = eig(H); 
        eE = diag(eE);


        nh = numel(hessdata);
        for param = fieldnames(hessdata)',
            hessdata(nh+1).(param{1}) = eval(param{1});
        end;

        fprintf('\run_quadfit35.m: finished evaluations at %i dimensions, with %i function evals...\n', n, size(xs,2));
        fprintf('eV = \n');
        display(eV);
        fprintf('eE = \n');
        display(eE);


        n = n - 1;
    end;

    bestE = [bestE; eE(:)];

    newV = eV;
    for k = numel(new_basis):-1:2,
        newV = new_basis{k}*newV;
    end;
    bestV = [bestV newV];

% if we somehow still wind up with negative eigenvalues...
if any(bestE < 0)
    fprintf('run_quadfit35.m: We still estimate at least one negative eigenvalue!\n');
    display(sort(bestE));
    fprintf('since this is (hopefully) a negative eigenvalue of small magnitude, \nand we know this is a peak, we"re going to totally hack it and pretend it"s positive\n');
    bestE = abs(bestE);
end;

fprintf('\n\n');
fprintf('run_quadfit35.m: Eigenvalues of estimate of H based on iterative approach:\n');
display(sort(bestE));
fprintf('with eigenvectors:\n');
display(bestV);

H = bestV*diag(bestE)*bestV';

fprintf('\n covariance matrix = inv(H) = \n');
display(inv(H));


%% save data to file
save(sprintf('quadfit_out_%s_off%i.mat', ratname, offset), 'ratname', ...
    'x_bf', 'hessdata', 'H', 'bestV', 'bestE');

%% do a line scan over B
B = 2:32;
LL_B = zeros(size(B));
LLgrad_B = zeros(numel(B), numel(x_bf));
do_param_linescan = [1 1 1 1 1 1 1 1 1];
fprintf(1, '\run_quadfit35.m: Doing linescan over B\n');
for i = 1:numel(B),
    xx = x_bf;
    xx(5) = B(i);
    [LL_B(i) LLgrad_B(i,:)] = ll_model_wrapper35(xx, mydata, 'for_fmin', 1, ...
        'total_rate', total_rate, ...
        'do_param', do_param_linescan);
end;

Bnd.B = B;
Bnd.LL_B = LL_B;
Bnd.LLgrad_B = LLgrad_B;

save(sprintf('quadfit_out_%s_off%i.mat', ratname, offset), 'ratname', ...
    'x_bf', 'hessdata', 'H', 'bestV', 'bestE', 'Bnd');

function d = make_deltas(N)

    deltas = eye(N, N);
    for i = 1:N-1,
        for j = i+1:N,
            step = zeros(N,1);
            step([i j]) = 1;
        deltas = [deltas step]; %#ok<AGROW>
    end;
end;

d = deltas;
return;