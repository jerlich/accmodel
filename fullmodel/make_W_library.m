function [W h_bit] = make_W_library(bo, bo_exp, lambda, sigma2_a, sigma2_s, dt, varargin)

pairs = { ...
    'N'            24  ; ... % clicks sizes will be quantized into 2^N bins
    };
parseargs(varargin, pairs);

W = struct('h', [], 'vh', [], 'F', [], ...
           'dFds', [], 'dFdsig2',  [], 'dFdB', [], ...
           'dFdl', [], 'dFdh', []);

% W(1): the contracting Markov matrix
%       takes whatever is in bins past B (in the extra, virtual bins) and
%       collapses them onto the sticky bounds
[F dFds dFdsig2 dFdB dFdl dFdh] = make_contracting_F(bo, bo_exp); %#ok<*ASGLU,*NASGU>
h  = 0;
vh = 0;
for d = fieldnames(W)',
    W(1).(d{1}) = eval(d{1});
end;

% W(2): to use in the absence of clicks
%       has sigma2_a and lambda
[F dFds dFdsig2 dFdB dFdl dFdh] = make_F2(bo_exp, sigma2_a*dt, lambda, 0, dt);
F = W(1).F*F; % contracts everything past B to the sticky bins
              % this operation has no effect on the derivatives (I hope)
h  = 0;
vh = sigma2_a*dt;
for d = fieldnames(W)',
    W(2).(d{1}) = eval(d{1});
end;

% W(3:42): the matrices to add variance for 1-40 clicks in this dt
%          has sigma2_a, sigma2_s, instability, and a contraction
for nclicks = 1:40;
    sigma2_tot = tot_var(sigma2_a, sigma2_s, nclicks, dt);
    [F dFds dFdsig2 dFdB dFdl dFdh] = make_F2(bo_exp, sigma2_tot, lambda, 0, dt);
    F = W(1).F*F;
    h  = 0;
    vh = sigma2_tot;
    for d = fieldnames(W)',
        W(nclicks+2).(d{1}) = eval(d{1});
    end;
end;

% W(43:43+N): the binarized matrices for rightward clicks
% W(43+N+1:43+N+1+N): the binarized matrices for leftward clicks
offset = numel(W);
height  = round(max(bins(bo_exp)));
for i = 1:N,
    % W up, for right clicks
    [F dFds dFdsig2 dFdB dFdl dFdh] = make_F2(bo_exp, 0, 0, height/dt, dt);
    h  = height;
    vh = 0;
    for d = fieldnames(W)',
        W(i+offset).(d{1}) = eval(d{1});
    end;
 

    % W down, for left clicks
    [F dFds dFdsig2 dFdB dFdl dFdh] = make_F2(bo_exp, 0, 0, -height/dt, dt);
    h  = -height;
    vh = 0;
    for d = fieldnames(W)',
        W(i+N+offset).(d{1}) = eval(d{1});
    end;

    height  = height/2;
end;    
h_bit  = height*2;

% this last one is used for the sole purpose of computing the derivative in
% alpha and rho
offset = numel(W);
[F dFds dFdsig2 dFdB dFdl dFdh] = make_F2(bo_exp, 0, 0, 0, dt);
h  = 0;
vh = 0;
for d = fieldnames(W)',
    W(1+offset).(d{1}) = eval(d{1});
end;
