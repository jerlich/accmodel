function [tval,pval,cival]=ttest2_w(X,Y)
if isempty(X) || isempty(Y)
    tval=nan; pval=1; cival=[-inf inf];
else
[b,pval,cival,t4]=ttest2(X,Y);

tval=t4.tstat;
end
