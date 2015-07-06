function Deltak = P_init_backwards(Deltak, P_init, Pb, W, bo)

n = nbins(bo);

F       = W.F;
dFdsig2 = W.dFdsig2;
dFdB    = W.dFdB;
P0      = W.P0;

% the backwards transition matrix
% B = F'.*(P0*ones(1,n))./(ones(n,1)*P_init');
% B(:,P_init==0) = 0;
% Pb0 = B*Pb(:,1);

JmF = (ones(n,1)*P0')./(P_init*ones(1,n)).*(Pb(:,1)*ones(1,n));   
JmF(isnan(JmF)) = 0;
JmF(F==0) = 0;

Deltak.sigma2_i = sum(sum(JmF.*dFdsig2));
Deltak.B        = [sum(sum(JmF.*dFdB)) Deltak.B];