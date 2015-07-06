function [p,D]=bootsigmoid(A,B,varargin)
% [p,D]=bootsigmoid(A,B,varargin)
% Takes two sets of N x 2 binomial data A and B, and fits sigmoids to them
% both.  It then uses the mean and covariance of the fits to calculate the
% distance (using the projection of the fits onto fisher's linear
% discriminant) of the fits.  
% Then we permute the rows of A and B to generate permuted data sets and
% perform the same fits and estimates of "distance". Finally, the distance
% between A and B fits are compared to the distribution of distances
% generated by permuting A and B.


BOOTS=10000;

overridedefaults(who,varargin);

n_A=size(A,1);
n_B=size(B,1);


[betaA,~,~,covA,~]=nlinfit(A(:,1),A(:,2),@sig4,[0 1 0 10]);
[betaB,~,~,covB,~]=nlinfit(B(:,1),B(:,2),@sig4,[0 1 0 10]);

vAB=flda(betaA,betaB,covA,covB,n_A,n_B);

dAB=abs(betaA*vAB-betaB*vAB);

permAB=nan(BOOTS,1);

M=[A;B];

parfor bx=1:BOOTS
    
    rperm=randperm(n_A+n_B);
    rA=M(rperm(1:n_A),:);
    rB=M(rperm((n_A+1):end),:);
    
    [rbetaA,~,~,rcovA,~]=nlinfit(rA(:,1),rA(:,2),@sig4,[0 1 0 10]);
    [rbetaB,~,~,rcovB,~]=nlinfit(rB(:,1),rB(:,2),@sig4,[0 1 0 10]);
    
    
    rvAB=flda(rbetaA,rbetaB,rcovA,rcovB,n_A,n_B);
    
    permAB(bx)=abs(rbetaA*rvAB-rbetaB*rvAB);
    
end

p=get_p(dAB,permAB);

D.ld=vAB;
D.betaA=betaA;
D.betaB=betaB;
D.permAB=permAB;
D.dAB=dAB;