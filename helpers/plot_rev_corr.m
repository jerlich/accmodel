function plot_rev_corr(ax,data,model,varargin)


dR = data(2).rev_corr_r;
dL = data(2).rev_corr_l;
mR = model(2).rev_corr_r;
mL = model(2).rev_corr_l;
t = data(2).tax;
M = {  dR dL mR mL};
clrs = { [1 0 0] [64 128 0]/255 [252 161 164]/255 [128 254 0]/255};
	
set(ax,'NextPlot','add');

for px=1:2
    lb = nanmean(M{px},1) - nanstderr(M{px},1);
    ub = nanmean(M{px},1) + nanstderr(M{px},1);
    
    shadeplot(t,lb,ub,{clrs{px},ax,0.5});
end

for px=3:4
   plot(ax,t,nanmean(M{px},1),'Color',clrs{px},'LineWidth',1); 
end

ylabel(ax,sprintf('Excess\nClick Rate (Hz)'));
xlabel(ax,sprintf('Time from\nend of stimulus (s)'));

xlim(ax,[-.65 0]);
ylim(ax,[-3.5 3.5]);
xhairs(ax,':k')



function y=lite(y)

y = (y+[1 1 1])/2;

