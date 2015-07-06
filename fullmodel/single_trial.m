function [LL dLL likey output] = single_trial(param, mydata, varargin)

pairs = { ...
    'dx'          0.25  ; ...
    'dt'          0.02  ; ...
    'bo'            []  ; ...
    'P_init'        []  ; ...
    'W_i'           []  ; ...
    'W0'            []  ; ...
    'ndelta'      1e-4  ; ...
    };
parseargs(varargin, pairs);

%% unpack click inputs for this trial
mydata = mydata(1);
T         = mydata.T;
leftbups  = mydata.leftbups;  leftbups = leftbups(leftbups < T);
rightbups = mydata.rightbups; rightbups = rightbups(rightbups < T);
pokedR    = mydata.pokedR;

%% unpack model parameters
lambda   = param(1);
sigma2_a = param(2);
sigma2_s = param(3);
sigma2_i = param(4);
B        = param(5);
alpha    = param(6);
rho      = param(7);
bias     = param(8);
inatt    = param(9);

%% arguments needed, if not passed in
if isempty(bo),
    bo = binning(dx, B);
end;

if isempty(P_init), %#ok<NODEF>
    [P_init W_i] = make_P_init(bo, sigma2_i, inatt);
end;

%% make Markov matrices for each time step
N  = ceil(T/dt);  % number of timesteps
t  = (0:N-1)*dt;  % a time axis

% clicks_L and clicks_R are the same size as leftbups and rightbups, respectively
[clicks_L clicks_R] = make_facdep_clicks(leftbups, rightbups, alpha, rho);
% net_input and nclicks are the size size as t
[net_input nclicks] = convert_clicks_to_inputs(t, leftbups, rightbups, clicks_L, clicks_R, 'is_facdep', abs(alpha-1)>eps);

% W is a cell array with one element for each time step, containing F and
% all the dFdtheta's used in the backwards, posterior run
if isempty(W0), %#ok<NODEF>
    W0 = make_each_F(bo, sigma2_a, sigma2_s, lambda, dt, 'net_input', 0, 'nclicks', 0);
    W0 = W0{1};
end;
W = make_each_F(bo, sigma2_a, sigma2_s, lambda, dt, 'net_input', net_input, 'nclicks', nclicks, 'W0', W0);


%% in order to compute derivatives w.r.t. alpha and rho,
% compute dh/dalpha and dh/drho numerically by finite difference,
% where h is the magnitude of the clicks
[clicks_L_dalpha clicks_R_dalpha] = make_facdep_clicks(leftbups, rightbups, alpha+ndelta, rho);
net_input_dalpha = convert_clicks_to_inputs(t, leftbups, rightbups, clicks_L_dalpha, clicks_R_dalpha, 'is_facdep', abs(alpha+ndelta-1)>eps);
net_input_dalpha = (net_input_dalpha - net_input)/ndelta;
[clicks_L_drho clicks_R_drho] = make_facdep_clicks(leftbups, rightbups, alpha, rho+ndelta);
net_input_drho = convert_clicks_to_inputs(t, leftbups, rightbups, clicks_L_drho, clicks_R_drho, 'is_facdep', abs(alpha-1)>eps);
net_input_drho = (net_input_drho - net_input)/ndelta;

%% forwards
Pf = zeros(nbins(bo), N);
Pf(:,1) = P_init; 

for k = 2:numel(t),
    Pf(:,k) = W{k-1}.F*Pf(:,k-1);
end;

Pd = one_sided_posterior(bo, Pf(:,end), bias, pokedR);
likey = sum(Pd); % likelihood of poking the same way the data did
LL = log(likey); % log of the likelihood

if nargout > 3,
    output.Pf = Pf;
    output.x  = bins(bo);
    output.t  = t;
end;


%% backwards
% normalize to use as initial distribution of the backwards, posterior run
Pd = Pd/sum(Pd); 

% Life must be lived forwards, but can only be understood backwards.
[Pb Deltak] = single_trial_backwards(Pf, Pd, bo, W);

% now take care of the virtual, 0th timestep that created P_init, and the
% derivative w.r.t. sigma2_i
Deltak = P_init_backwards(Deltak, Pf(:,1), Pb, W_i, bo);

if nargout > 3,
    output.Pb = Pb;
end;

%% compute gradient

% lambda
    % dL/dlambda = dL/dphi * dphi/dlambda + dL/dg * dgamma/dlambda
    %            = dL/dphi * (-c/lambda^2) + dL/dg * (gamma*dt)

    c = net_input'/dt; % the effective input per timestep
    g = exp(lambda*dt)*dt;

    dLL(1) = sum(Deltak.phi.*c) / (-lambda^2) + sum(Deltak.gamma)*g;

% sigma2_a
    % dL/dsigma2_a = dL/dsigma2 * dsigma2/dsigma2_a
    %              = dL/dsigma2 * dt

    dLL(2) = sum(Deltak.sigma2) * dt;

% sigma2_s
    % dL/dsigma2_s = dL/dsigma2 * dsigma2/dsigma2_s
    %              = dL/dsigma2 * nclicks / 41

    dLL(3) = sum(Deltak.sigma2.*nclicks')/41;

% sigma2_i
    dLL(4) = Deltak.sigma2_i;
    
% B
    dLL(5) = sum(Deltak.B);

% alpha
    % dL/dalpha = dL/dphi * dphi/dh * dh/dalpha;  

    dLL(6) = sum(Deltak.phi.*net_input_dalpha')/lambda/dt; 

% rho
    % dL/drho = dL/dphi * dphi/dh * dh/drho;  

    dLL(7) = sum(Deltak.phi.*net_input_drho')/lambda/dt;      

% bias
    % compute deriv numerically,
    LLepsilon = one_sided_posterior(bo,Pf(:,end),bias+ndelta,pokedR);
    dLL(8) = (log(sum(LLepsilon)) - LL)/ndelta;

% inatt
    % compute deriv numerically,
    dLL(9) = (log((likey-inatt/2)/(1-inatt)*(1-(inatt+ndelta))+(inatt+ndelta)/2) - log(likey))/ndelta;
    
    
% single trial done!    