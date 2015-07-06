function [Pb Deltak] = single_trial_backwards(Pf, Pd, bo, W)

% the backwards, posterior probability distribution
Pb = zeros(size(Pf));

n = nbins(bo);
N = size(Pf,2); % number of time steps

Deltak.gamma  = zeros(1,N);
Deltak.phi    = zeros(1,N);
Deltak.sigma2 = zeros(1,N);
Deltak.B      = zeros(1,N);

Pb(:,end) = Pd; % last timestep of backwards distribution is Pd
for k = N:-1:2,
    F = W{k-1}.F;
    
    B1 = backwards_matrix(F, Pf(:,k-1), Pf(:,k));
    B = F'.*(Pf(:,k-1)*ones(1,n))./(ones(n,1)*Pf(:,k)');  % The backwards transition matrix, from the notes
    B(isnan(B))=0;
    Pb(:,k-1) = B*Pb(:,k); % Iterate one step backward
    
    % J = B'.*(Pb(:,k)*ones(1, numel(x))); %#ok<NASGU>  This line not used, here for annotation purposes only
    % J2 = (F'.*(Pf(:,k-1)*ones(1,numel(x)))./(ones(numel(x),1)*Pf(:,k)'))'.*(Pb(:,k)*ones(1, numel(x)));
    % J3 = (F .*(ones(numel(x),1)*Pf(:,k-1)')./(Pf(:,k)*ones(1, numel(x)))).*(Pb(:,k)*ones(1, numel(x)));

    % Now compute JmF, "J minus F", defined as JmF_ij = J_ij/F_ij.  Note that when we
    % compute J_ij * d(ln F_ij)/dtheta = J_ij/F_ij * dF_ij/dtheta = JmF_ij * dF_ij/dtheta 
    JmF = (ones(n,1)*Pf(:,k-1)')./(Pf(:,k)*ones(1, n)).*(Pb(:,k)*ones(1, n));   
    JmF(isnan(JmF)) = 0;
    JmF(F==0) = 0;

    Deltak.gamma(k-1)  = sum(sum(JmF.*W{k-1}.dFdg));
    Deltak.phi(k-1)    = sum(sum(JmF.*W{k-1}.dFdphi));
    Deltak.sigma2(k-1) = sum(sum(JmF.*W{k-1}.dFdsig2));
    Deltak.B(k-1)      = sum(sum(JmF.*W{k-1}.dFdB));
end;