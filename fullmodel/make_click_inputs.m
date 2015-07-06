function [net_input nclicks] = make_click_inputs(t, leftbups, rightbups, clicks_L, clicks_R, NL, NR)
% function [net_input nclicks] = make_click_inputs(t, leftbups, rightbups, clicks_L, clicks_R, NL, NR)
% returns two vectors the same length as t, representing the input in each
% timestep
%
% net_input is the effective input (in click units), clicks_R - clicks_L
% for each dt
%
% nclicks is the total number of clicks (left + right) in each dt; 

here_L = qfind(t, leftbups);
here_R = qfind(t, rightbups);
if nargout > 1,
    nclicks = zeros(length(t),1);
    for i = 1:numel(here_L),
        nclicks(here_L(i)) = nclicks(here_L(i)) + NL(i);
    end;
    for i = 1:numel(here_R),
        nclicks(here_R(i)) = nclicks(here_R(i)) + NR(i);
    end;
end;


net_input = zeros(length(t),1);
for i = 1:numel(leftbups),
    net_input(here_L(i)) = net_input(here_L(i)) - clicks_L(i);
end;
for i = 1:numel(rightbups),
    net_input(here_R(i)) = net_input(here_R(i)) + clicks_R(i);
end;