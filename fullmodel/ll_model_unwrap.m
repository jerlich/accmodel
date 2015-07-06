function param = ll_model_unwrap(rsparam, varargin)
%
%
% This function complements ll_model_wrapper and returns sensible param values
%
% Rescaling is done such that the range 0 to 1 in rsparam is mapped to 
% lower_bound to upper_bound for each parameter.

% x_names = {'lambda' 'sigma2_a'  'sigma2_s'  'sigma2_i' 'B'      'alpha'     'rho'   'bias'  'inatt'};

pairs = { 
   'wrap'        0; ...
   'lower_bound'           [-5  0    0    0    2     0.1     0.005     -5     0  ]; ...
   'upper_bound'           [+5  400  400  20   32    1.5     1.005     +5     1  ]; ...
}; parseargs(varargin, pairs);  
   
if ischar(rsparam),
   switch rsparam
      case 'get_lower_bound',
         param = lower_bound;
         return;
      case 'get_upper_bound',
         param = upper_bound;
         return;
      otherwise
         error('ll_model_unwrap does not recognize call ''%s''', rsparam);
   end;
end;


if wrap == 0,    
    param = rsparam.*(upper_bound-lower_bound) + lower_bound;
else
    param = (rsparam-lower_bound)./(upper_bound-lower_bound);
end;