function [LL dLL likey output output_exp] = single_trial2(param, mydata, varargin)
% new in version 2: uses binarized Markov matrices and virtual timesteps!
%
% BWB, Aug 2011

pairs = { ...
    'dx'          0.25  ; ...
    'dt'          0.02  ; ...
    'bo'            []  ; ... % bin object
    'bo_exp'        []  ; ... % expanded bin object, with virtual bins
    'P_init'        []  ; ...
    'W_i'           []  ; ...
    'W'             []  ; ...
    'h_bit'         []  ; ...
    'Nclickbins'    24  ; ... % clicks sizes will be quantized into 2^24 bins
    'ndelta'      1e-5  ; ...
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
if isempty(bo), %#ok<NODEF>
    bo = binning(dx, B);
    bo_exp  = binning(dx, 2*B);
end;

if isempty(P_init), %#ok<NODEF>
    [P_init W_i] = make_P_init(bo_exp, sigma2_i, inatt);
end;

% the library of Markov matrices to handle timesteps with clicks
if isempty(W), %#ok<NODEF>
    [W h_bit] = make_W_library(bo, bo_exp, lambda, sigma2_a, sigma2_s, dt, 'N', Nclickbins);
end;

%% approximate Markov process with a series of virtual timesteps and contractions
N  = ceil(T/dt);  % number of timesteps
t  = (0:N-1)*dt;  % a time axis

% clicks_L and clicks_R are the same size as leftbups and rightbups, respectively
[clicks_L clicks_R] = make_facdep_clicks(leftbups, rightbups, alpha, rho);
% net_input and nclicks are the size size as t
[net_input nclicks] = convert_clicks_to_inputs(t, leftbups, rightbups, clicks_L, clicks_R, 'is_facdep', abs(alpha-1)>eps);


% in order to compute derivatives w.r.t. alpha and rho:
% compute dh/dalpha and dh/drho numerically by finite difference,
% where h is the magnitude of the clicks

[clicks_L_dalpha clicks_R_dalpha] = make_facdep_clicks(leftbups, rightbups, alpha+ndelta, rho);
net_input_dalpha = convert_clicks_to_inputs(t, leftbups, rightbups, clicks_L_dalpha, clicks_R_dalpha, 'is_facdep', abs(alpha+ndelta-1)>eps);
net_input_dalpha = (net_input_dalpha - net_input)/ndelta;
[clicks_L_drho clicks_R_drho] = make_facdep_clicks(leftbups, rightbups, alpha, rho+ndelta);
net_input_drho = convert_clicks_to_inputs(t, leftbups, rightbups, clicks_L_drho, clicks_R_drho, 'is_facdep', abs(alpha-1)>eps);
net_input_drho = (net_input_drho - net_input)/ndelta;


% figure out which matrix to use when and save the indices in vector m
m  = [];
mc = []; % same size as m, contains nclicks for each virtual step
dhda = []; % dh/dalpha
dhdr = []; % dh/drho
is_real = []; % same size as m, marks the real timesteps
for k = 1:numel(net_input),
    % at each real timestep:
    if nclicks(k) == 0,
        m = [m 2]; % W(2): timesteps with no clicks
        mc = [mc 0];
        dhda = [dhda 0]; dhdr = [dhdr 0];
        is_real = [is_real 1];
    else
        h  = net_input(k); 

        if abs(h) > eps,
            % first reproduce the desired click height
            str = dec2bin(round(abs(h)/h_bit));
            h_approx = 0;

            hct = 42; % HACK ALERT: this is the index before the first click height matrix
            offset = Nclickbins - numel(str);
            for i = numel(str):-1:1,
                if str(i)=='1',
                    if h > 0,
                        m = [m hct+offset+i];
                    else
                        m = [m hct+offset+Nclickbins+i];
                    end;
                    h_approx = h_approx + W(m(end)).h;            
                    mc = [mc 0];
                    dhda = [dhda 0]; dhdr = [dhdr 0];
                    is_real = [is_real 0];
                end;
            end;
            % fprintf(1, 'h=%g, approx=%g, (h-approx)/h=%g\n', h, h_approx, (h-h_approx)/h);
        end;
        
%         m = [m numel(W)];
%         mc = [mc 0];
        dhda = [dhda net_input_dalpha(k)]; dhdr = [dhdr net_input_drho(k)];
%         is_real = [is_real 0];
        
        % next use the correct matrix that reproduces the desired variance
        % and instability
        m = [m 2+nclicks(k)];
        mc = [mc nclicks(k)];
%         dhda = [dhda 0]; dhdr = [dhdr 0];
        is_real = [is_real 1];
    end;
end;

mc = [mc 0]; % HACK ALERT: padded because there's an extra step at the end
dhda = [dhda 0]; dhdr = [dhdr 0];


%% forwards
Pf = zeros(nbins(bo_exp), numel(m)+1);
Pf(:,1) = P_init; 

% iterate through all real and virtual timesteps
for k = 2:numel(m)+1,
    Pf(:,k) = W(m(k-1)).F*Pf(:,k-1);
end;

Pd = one_sided_posterior(bo_exp, Pf(:,end), bias, pokedR);
likey = sum(Pd); % likelihood of poking the same way the data did
LL = log(likey); % log of the likelihood

if nargout > 4,
    output_exp.Pf = Pf;
    output_exp.x  = bins(bo_exp);
    output_exp.t  = t;
    output_exp.W  = W;
end;



dLL = zeros(size(param));

%% backwards
% normalize to use as initial distribution of the backwards, posterior run
Pd = Pd/sum(Pd); 

Pb = zeros(size(Pf));
Pb(:,end) = Pd;
n = nbins(bo_exp);

Deltak.lambda = zeros(1, size(Pf,2));
Deltak.h      = zeros(1, size(Pf,2));
Deltak.sigma2 = zeros(1, size(Pf,2));
Deltak.B      = zeros(1, size(Pf,2));

% Life must be lived forwards, but can only be understood backwards.
for k = numel(m)+1:-1:2,
    B = W(m(k-1)).F'.*(Pf(:,k-1)*ones(1,n))./(ones(n,1)*Pf(:,k)');  % The backwards transition matrix, from the notes
    B(isnan(B))=0;
    Pb(:,k-1) = B*Pb(:,k); % Iterate one step backward
    
    JmF = (ones(n,1)*Pf(:,k-1)')./(Pf(:,k)*ones(1, n)).*(Pb(:,k)*ones(1, n));   
    JmF(isnan(JmF)) = 0;
    JmF(W(m(k-1)).F==0) = 0;

    Deltak.lambda(k-1) = sum(sum(JmF.*W(m(k-1)).dFdl));
    Deltak.h(k-1)      = sum(sum(JmF.*W(m(k-1)).dFdh));
    Deltak.sigma2(k-1) = sum(sum(JmF.*W(m(k-1)).dFdsig2));
    Deltak.B(k-1)      = sum(sum(JmF.*W(m(k-1)).dFdB));
end;

% now take care of the virtual, 0th timestep that created P_init, and the
% derivative w.r.t. sigma2_i
JmF = (ones(n,1)*W_i.P0')./(P_init*ones(1,n)).*(Pb(:,1)*ones(1,n));   
JmF(isnan(JmF)) = 0;
JmF(W_i.F==0) = 0;

Deltak.sigma2_i = sum(sum(JmF.*W_i.dFdsig2));
Deltak.B        = [sum(sum(JmF.*W_i.dFdB)) Deltak.B];



Pf_real = zeros(nbins(bo), N); % the real slim Pf, without all the virtual bins and virtual timesteps
Pb_real = zeros(nbins(bo), N);
real_dt = find(is_real==1);
for k = 1:N,
    Pf_real(:,k) = collapse_P(Pf(:,real_dt(k)), bo, bo_exp);
    Pb_real(:,k) = collapse_P(Pb(:,real_dt(k)), bo, bo_exp);
end;

if nargout > 3,
    output.Pf = Pf_real;
    output.Pb = Pb_real;
    output.x  = bins(bo);
    output.t  = t;
end;

if nargout > 4,
    output_exp.Pb = Pb;
end;

%% compute gradient

% lambda
    % dL/dlambda = dL/dphi * dphi/dlambda + dL/dg * dgamma/dlambda
    %            = dL/dphi * (-c/lambda^2) + dL/dg * (gamma*dt)

    dLL(1) = sum(Deltak.lambda);

% sigma2_a
    % dL/dsigma2_a = dL/dsigma2 * dsigma2/dsigma2_a
    %              = dL/dsigma2 * dt

    dLL(2) = sum(Deltak.sigma2) * dt;

% sigma2_s
    % dL/dsigma2_s = dL/dsigma2 * dsigma2/dsigma2_s
    %              = dL/dsigma2 * nclicks / 41

    dLL(3) = sum(Deltak.sigma2.*mc)/41;

% sigma2_i
    dLL(4) = Deltak.sigma2_i;
    
% B
    dLL(5) = sum(Deltak.B);

% alpha
    % dL/dalpha = dL/dh * dh/dalpha;  

    dLL(6) = sum(Deltak.h.*dhda)/dt;

% rho
    % dL/drho = dL/dh * dh/drho

    dLL(7) = sum(Deltak.h.*dhdr)/dt;      

% bias
    % compute deriv numerically,
    LLepsilon = one_sided_posterior(bo_exp,Pf(:,end),bias+ndelta,pokedR);
    dLL(8) = (log(sum(LLepsilon)) - LL)/ndelta;

% inatt
    % compute deriv numerically,
    dLL(9) = (log((likey-inatt/2)/(1-inatt)*(1-(inatt+ndelta))+(inatt+ndelta)/2) - log(likey))/ndelta;
    
%     
% single trial done!    