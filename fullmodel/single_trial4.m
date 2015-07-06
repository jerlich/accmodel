function [LL dLL likey output] = single_trial4(param, mydata, varargin)
% new in version 4: has cross-side adaptation, and sensory noise scales with
% magnitude of clicks
% June 2012, BWB


pairs = { ...
    'dx'          0.25  ; ...
    'dt'          0.02  ; ...
    'bo'            []  ; ...
    'P_init'        []  ; ...
    'W_i'           []  ; ...
    'W0'            []  ; ...
    'ndelta'      1e-4  ; ... % delta used for numerical finite-difference estimates of derivatives
    'cross_side_suppression' 0 ; 
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
phi      = param(6);
tau_phi  = param(7);
bias     = param(8);
inatt    = param(9);
psi      = param(10);
tau_psi  = param(11);

%% arguments needed, if not passed in
if isempty(bo), %#ok<NODEF>
    bo = binning(dx, B);
end;

if isempty(P_init), %#ok<NODEF>
    [P_init W_i] = make_P_init3(bo, sigma2_i, inatt);
end;

%% make Markov matrices for each time step
N  = ceil(T/dt);  % number of timesteps
t  = (0:N-1)*dt;  % a time axis

% clicks_L and clicks_R are the same size as leftbups and rightbups, respectively
[clicks_L clicks_R NL NR] = make_adapted_clicks(leftbups, rightbups, phi, tau_phi, psi, tau_psi, ...
                                          'cross_side_suppression', cross_side_suppression);
% net_input and nclicks are the same size as t
[net_input tot_input] = make_click_inputs35(t, leftbups, rightbups, clicks_L, clicks_R, NL, NR);

% W is a cell array with one element for each time step, containing F and
% all the dFdtheta's used in the backwards, posterior run
if isempty(W0), %#ok<NODEF>
    W0 = make_each_F35(bo, sigma2_a, sigma2_s, lambda, dt, 'net_input', 0, 'tot_input', 0);
end;
W = make_each_F35(bo, sigma2_a, sigma2_s, lambda, dt, 'net_input', net_input, 'tot_input', tot_input, 'W0', W0);


%% in order to compute derivatives w.r.t. phi and tau_phi,
% compute dh/dphi and dh/dtauphi numerically by finite difference,
% where h is the magnitude of the clicks

% dh/dphi
[clicks_L_dphi clicks_R_dphi NL NR] = make_adapted_clicks(leftbups, rightbups, phi+ndelta, tau_phi, psi, tau_psi, ...
                                          'cross_side_suppression', cross_side_suppression);
[net_input_dphi tot_input_dphi] = make_click_inputs35(t, leftbups, rightbups, clicks_L_dphi, clicks_R_dphi, NL, NR);
net_input_dphi = (net_input_dphi - net_input)/ndelta;
tot_input_dphi = (tot_input_dphi - tot_input)/ndelta;

% dh/dtau_phi
[clicks_L_dtauphi clicks_R_dtauphi NL NR] = make_adapted_clicks(leftbups, rightbups, phi, tau_phi+ndelta, psi, tau_psi, ...
                                          'cross_side_suppression', cross_side_suppression);
[net_input_dtauphi tot_input_dtauphi] = make_click_inputs35(t, leftbups, rightbups, clicks_L_dtauphi, clicks_R_dtauphi, NL, NR);
net_input_dtauphi = (net_input_dtauphi - net_input)/ndelta;
tot_input_dtauphi = (tot_input_dtauphi - tot_input)/ndelta;

% dh/dpsi
[clicks_L_dpsi clicks_R_dpsi NL NR] = make_adapted_clicks(leftbups, rightbups, phi, tau_phi, psi+ndelta, tau_psi, ...
                                          'cross_side_suppression', cross_side_suppression);
[net_input_dpsi tot_input_dpsi] = make_click_inputs35(t, leftbups, rightbups, clicks_L_dpsi, clicks_R_dpsi, NL, NR);
net_input_dpsi = (net_input_dpsi - net_input)/ndelta;
tot_input_dpsi = (tot_input_dpsi - tot_input)/ndelta;

% dh/dtau_psi
[clicks_L_dtaupsi clicks_R_dtaupsi NL NR] = make_adapted_clicks(leftbups, rightbups, phi, tau_phi, psi, tau_psi+ndelta, ...
                                          'cross_side_suppression', cross_side_suppression);
[net_input_dtaupsi tot_input_dtaupsi] = make_click_inputs35(t, leftbups, rightbups, clicks_L_dtaupsi, clicks_R_dtaupsi, NL, NR);
net_input_dtaupsi = (net_input_dtaupsi - net_input)/ndelta;
tot_input_dtaupsi = (tot_input_dtaupsi - tot_input)/ndelta;

%% forwards
Pf = zeros(nbins(bo), N);
Pf(:,1) = P_init; 

for k = 2:numel(t),
    Pf(:,k) = W(k-1).F*Pf(:,k-1);
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
[Pb Deltak] = single_trial_backwards2(Pf, Pd, W, net_input);

% now take care of the virtual, 0th timestep that created P_init, and the
% derivative w.r.t. sigma2_i
Deltak = P_init_backwards(Deltak, Pf(:,1), Pb, W_i, bo);

if nargout > 3,
    output.Pb = Pb;
end;

%% compute gradient

% lambda
    dLL(1) = sum(Deltak.l);

% sigma2_a
    % dL/dsigma2_a = dL/dsigma2 * dsigma2/dsigma2_a
    %              = dL/dsigma2 * dt

    dLL(2) = sum(Deltak.sig2) * dt;

% sigma2_s
    % dL/dsigma2_s = dL/dsigma2 * dsigma2/dsigma2_s
    %              = dL/dsigma2 * tot_input / 40

    dLL(3) = sum(Deltak.sig2.*tot_input')/40;

% sigma2_i
    dLL(4) = Deltak.sigma2_i;
    
% B
    dLL(5) = sum(Deltak.B);

% phi
    % dL/dphi = dL/dh * dh/dphi + dL/dsig2 * dsig2/dphi
    %         = dL/dh * dh/dphi + dL/dsig2 * sigma2_s/40 * dH/dphi
    % where h is net_input, and H is tot_input

    dLL(6) = sum(Deltak.h.*net_input_dphi')/dt + ...
             sum(Deltak.sig2.* tot_input_dphi')*sigma2_s/40; 

% tau_phi
    % dL/dtauphi = dL/dh * dh/dtauphi + dL/dsig2 * dsig2/dtauphi
    %            = dL/dh * dh/dtauphi + dL/dsig2 * sigma2_s/40 * dH/dtauphi

    dLL(7) = sum(Deltak.h.*net_input_dtauphi')/dt+ ...
             sum(Deltak.sig2.* tot_input_dtauphi')*sigma2_s/40;      

% bias
    % compute deriv numerically,
    LLepsilon = one_sided_posterior(bo,Pf(:,end),bias+ndelta,pokedR);
    dLL(8) = (log(sum(LLepsilon)) - LL)/ndelta;

% inatt
    % compute deriv numerically,
    dLL(9) = (log((likey-inatt/2)/(1-inatt)*(1-(inatt+ndelta))+(inatt+ndelta)/2) - log(likey))/ndelta;
    
% psi
    % dL/dphi = dL/dh * dh/dpsi + dL/dsig2 * dsig2/dpsi
    %         = dL/dh * dh/dpsi + dL/dsig2 * sigma2_s/40 * dH/dpsi
    % where h is net_input, and H is tot_input
    dLL(10) = sum(Deltak.h.*net_input_dpsi')/dt + ...
              sum(Deltak.sig2.* tot_input_dpsi')*sigma2_s/40;

% tau_psi
    % dL/dtaupsi = dL/dh * dh/dtaupsi + dL/dsig2 * dsig2/dtaupsi
    %            = dL/dh * dh/dtaupsi + dL/dsig2 * sigma2_s/40 * dH/dtaupsi
    
    dLL(11) = sum(Deltak.h.*net_input_dtaupsi')/dt+ ...
              sum(Deltak.sig2.* tot_input_dtaupsi')*sigma2_s/40;   
    
% single trial done!    