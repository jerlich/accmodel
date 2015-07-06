function [Pb Deltak] = single_trial_backwards2(Pf, Pd, W, net_input)

% the backwards, posterior probability distribution
Pb = zeros(size(Pf));

N = size(Pf,2); % number of time steps

Deltak.l    = zeros(1,N);
Deltak.h    = zeros(1,N);
Deltak.sig2 = zeros(1,N);
Deltak.B    = zeros(1,N);

Pb(:,end) = Pd; % last timestep of backwards distribution is Pd
for k = N:-1:2,
    F = W(k-1).F;
    
    B = backwards_matrix(F, Pf(:,k-1), Pf(:,k));
% The code for backwards_matrix is missing. :(


%     B = F'.*(Pf(:,k-1)*ones(1,n))./(ones(n,1)*Pf(:,k)');  % The backwards transition matrix, from the notes


    B(isnan(B))=0;
    Pb(:,k-1) = B*Pb(:,k); % Iterate one step backward
    
    % J = B'.*(Pb(:,k)*ones(1, numel(x))); %#ok<NASGU>  This line not used, here for annotation purposes only
    % J2 = (F'.*(Pf(:,k-1)*ones(1,numel(x)))./(ones(numel(x),1)*Pf(:,k)'))'.*(Pb(:,k)*ones(1, numel(x)));
    % J3 = (F .*(ones(numel(x),1)*Pf(:,k-1)')./(Pf(:,k)*ones(1, numel(x)))).*(Pb(:,k)*ones(1, numel(x)));

    % Now compute JmF, "J minus F", defined as JmF_ij = J_ij/F_ij.  Note that when we
    % compute J_ij * d(ln F_ij)/dtheta = J_ij/F_ij * dF_ij/dtheta = JmF_ij * dF_ij/dtheta 
    JmF = JmF_matrix(Pf(:,k-1),Pf(:,k),Pb(:,k));
    % JmF = (ones(n,1)*Pf(:,k-1)')./(Pf(:,k)*ones(1, n)).*(Pb(:,k)*ones(1, n));   
    JmF(isnan(JmF) | isinf(JmF)) = 0;
    JmF(F==0) = 0;
    
    Deltak.l(k-1)    = sum(sum(JmF.*W(k-1).dFdl));
    if abs(net_input(k-1)) > 0,
        Deltak.h(k-1) = sum(sum(JmF.*W(k-1).dFdh));
    else % if there's no net input at this time step, no need wasting time 
        Deltak.h(k-1) = 0;
    end;
    Deltak.sig2(k-1) = sum(sum(JmF.*W(k-1).dFdsig2));
    Deltak.B(k-1)    = sum(sum(JmF.*W(k-1).dFdB));
end;