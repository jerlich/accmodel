function [F dFds dFdsig2 dFdB dFdl dFdh] = make_F2(bo, sigma2, lambda, h, dt, varargin)
% version 2 takes lambda, h, and dt as inputs directly, instead of gamma
% and phi; correspondingly, dFdl and dFdh are returned instead of dFdg and
% dFdphi

x  = bins(bo);
dB = dxdB(bo);

% Declare F, the forwards markov matrix, and its derivatives:
F       = zeros(numel(x), numel(x));  % This will be F_ij
dFds    = zeros(numel(x), numel(x));  % This will be dF_ij / ds, where s is the center of mass of a dollop of distribution that came from bin j
dFdsig2 = zeros(numel(x), numel(x));  % dF_ij/dsigma2
dFdB    = zeros(numel(x), numel(x));  % dF_ij/dB
F(1,1)  = 1; F (end,end)=1;

% define the slices of a gaussian with sigma2 variance
[deltas ps] = gaussbins(binsize(bo), sigma2);

% construct F and its derivatives, stepping through from the bins j we came
% from:
for j = 2:numel(x)-1,
    % mu is the center of mass of where the probability starting at bin j will go
    if abs(lambda) < eps
        mu = x(j) + h*dt;
    else
        mu = exp(lambda*dt)*(x(j) + h/lambda) - h/lambda;
    end;
    
    ss = mu + deltas; % these are the centers of mass of the slices we're looking at
    [hps, lps] = ceil_and_floor(bo, ss);
    % we're going to look over all of the slices of the gaussian
    for k = 1:numel(deltas),
%         s = mu + deltas(k); % this is the center of mass of the slice we're looking at
        s = ss(k);
        if s <= x(1),       % if s is more negative than the first bin, put it all in the first bin
            F(1,   j) = F(1, j) + ps(k);
            % there's nothing to add to the derivatives
        elseif s >= x(end), % similarly, if s is more positive than the last bin
            F(end, j) = F(end, j) + ps(k);
        else
            % find the bin id of the lowest bin whose position is greater
            % than or equal to s (= hp), and the bin id of the highest bin
            % whose position is less than or equal to s (= lp)
%             [hp, lp] = ceil_and_floor(bo, s);
            hp = hps(k); lp = lps(k);
            if ~(x(lp)-eps<=s && s<=x(hp)+eps), error('woops!'); end;
            dd = x(hp) - x(lp);
            if dd == 0, % if s landed _exactly_ on a bin
                F(lp, j) = F(lp, j) + ps(k);
                
                dFds(lp,j)   = dFds(lp,j)   - ps(k)/(x(lp+1)-x(lp)); % if s grows further, this bin will lose probability
                dFds(lp+1,j) = dFds(lp+1,j) + ps(k)/(x(lp+1)-x(lp)); % and the next bin up will gain a corresponding amount
                
                % deltas(k)/(2*sigma2) is the derivative of s with respect to sigma2.  (dF/ds)*(ds/dsigma2) = dF/dsigma2.
                dFdsig2(lp,  j) = dFdsig2(lp,  j) - (ps(k)/(x(lp+1)-x(lp)))*deltas(k)/(2*sigma2);  
                dFdsig2(lp+1,j) = dFdsig2(lp+1,j) + (ps(k)/(x(lp+1)-x(lp)))*deltas(k)/(2*sigma2);
            else
                % if s did not land exactly on one bin, its probability
                % gets distributed to the two immediate neighbors
                F(hp,j) = F(hp,j) + ps(k)*(s-x(lp))/dd;      % amount of probability assigned to hp   
                F(lp,j) = F(lp,j) + ps(k)*(x(hp)-s)/dd;      % amount of probability assigned to lp   

                dFds(hp,j) = dFds(hp,j) + ps(k)/dd;     % now we just differentiate the above expressions
                dFds(lp,j) = dFds(lp,j) - ps(k)/dd;

                % again, deltas(k)/(2*sigma2) is the derivative of s with respect to sigma2.  (dF/ds)*(ds/dsigma2) = dF/dsigma2.
                dFdsig2(hp,j) = dFdsig2(hp,j) + (ps(k)/dd)*deltas(k)/(2*sigma2);   
                dFdsig2(lp,j) = dFdsig2(lp,j) - (ps(k)/dd)*deltas(k)/(2*sigma2);

                if lp==1, % we're at the bound; compute derivative with respect to x(lp)
                   dFdB(lp,j) = dFdB(lp,j) + dB(1)*ps(k)*(x(hp)-s)/dd.^2;
                   dFdB(hp,j) = dFdB(hp,j) + dB(1)*ps(k)*((s - x(lp)) - dd)/dd^2;
                end;
                if hp==numel(x),
                   dFdB(lp,j) = dFdB(lp,j) + dB(end)*ps(k)*(dd - (x(hp)-s))/dd.^2;
                   dFdB(hp,j) = dFdB(hp,j) - dB(end)*ps(k)*(s - x(lp))/dd^2;               
                end;
            end;
        end;
    end;
end;

if sigma2 < eps,
    dFdsig2 = 0;
end;

% The derivative of F with respect to lambda and c both follow from the
% derivative w.r.t. s
if abs(lambda)<eps,
    dFdl = 0;
    dFdh = dFds*dt;
else
    dFdl = dFds.*(ones(numel(x),1)*x + h/lambda)*dt*exp(lambda*dt)  ...
           - dFds*h/lambda^2*(exp(lambda*dt)-1);
    dFdh = dFds*(exp(lambda*dt)-1)/lambda;
end;