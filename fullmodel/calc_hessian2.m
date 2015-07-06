function [H, fval, grad, dfval, dgrad] = calc_hessian2(func, x, delta)
% function [H, fval, grad, dfval, dgrad] = calc_hessian2(func, x, delta)
%
% numerical evalutes the hessian of func at x
% assumes func's first two outputs are its functional value and its
% gradient

% compute grad at x
[fval grad] = func(x);

H = zeros(numel(x), numel(x));

dfval = zeros(size(x));
dgrad  = zeros(numel(x), numel(x));

if numel(delta)==1,
    delta = delta*ones(size(x));
end;

% compute grad at a step of delta in each direction of x
% dgrad(i,j) = gradient w.r.t. param j when param i has been incremented by delta
for i=1:numel(x),
   dx = zeros(size(x)); 
   dx(i) = delta(i);
   [dfval(i) dgrad(i,:)] = func(x+dx); 
   fprintf(1, 'Computed %d/%d gradients\n', i, numel(x));
end;

% calculate the hessian from the gradients
for i=1:numel(x),
    for j=i:numel(x),
        H(i,j) = (dgrad(j,i) - grad(i))/delta(j);
        H(j,i) = H(i,j);
    end;
end;

