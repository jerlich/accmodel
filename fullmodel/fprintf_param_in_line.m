function fprintf_param_in_line(x, x_names, a)

if nargin < 3, a = 1; end;

for i = 1:numel(x)
    fprintf(a, '-- %s = %8.6f  ', x_names{i}, x(i));
end;