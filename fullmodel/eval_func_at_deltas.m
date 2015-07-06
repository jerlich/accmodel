function [xs, Ls, grads] = eval_func_at_deltas(func, x, deltas)
% function [xs, Ls, grads] = eval_func_at_deltas(func, x, deltas)
%
% evaluates func at x and at x+deltas
% 
% INPUTS:
%   func    function which takes x as an input
%           it is assumed that the first argout is the functional value
%           and the second argout is the gradient
%
%   x       the input to func of size 1 by n (or n by 1)
%
%   deltas  the deltas around x at which we want to evalute func
%           if numel(deltas)==1, then take a step of that magnitude in each
%           direction while keeping all other parameters constant;
%           if numel(deltas)==numel(x), then take steps of variable
%           magnitudes in each direction of x;
%           else, deltas must be a matrix of n by m, where n is the number
%           of elements in x and m is the number of different steps to
%           take, where each column in deltas correspond to a different
%           step to be taken
%
% OUTPUTS:
%   xs      the x values tried as functional inputs
%           a matrix of n by m
%
%   Ls      the function values of func evaluted at xs
%           a vector of 1 by m
%
%   grads   the gradients of func at each point xs
%           a matrix of n by m
%
%
%
% BWB, Oct 2011


% make x a column vector
x = x(:);

% display(x);
% display(deltas); 

if numel(deltas)==1,
    deltas = deltas*ones(size(x));
end;

if numel(deltas)==numel(size(x));
    deltas = diag(deltas);
end;

% set up xs at which to evaluate func
xs = repmat(x, 1, size(deltas,2)+1);
xs(:, 2:end) = xs(:,2:end) + deltas;

Ls = ones(1, size(xs,2));
if nargout > 2,
    grads = ones(size(xs));
else
    grads = [];
end;

for i = 1:size(xs,2),
    display([i xs(:,i)']);
    [Ls(i) temp] = func(xs(:,i));
    if nargout > 2,
        grads(:,i) = temp;
    end;
end;

