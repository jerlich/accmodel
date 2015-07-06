function r2=r2tjur(dep,ind)

% r2 = r2kjur(dep,ind)
%
% dep   a float vector of model predictions (-value predict 0 and +values
%       predict 1)
% ind   a binomal vector of predictions
%
% Based on 
% Tjur, T. (2009). Coefficients of determination in logistic regression models?A new proposal: The coefficient of discrimination. The American Statistician.
% 
% Discovered here:
% http://www.statisticalhorizons.com/r2logistic
%


r2 = mean(ind(dep>0))-mean(ind(dep<0));