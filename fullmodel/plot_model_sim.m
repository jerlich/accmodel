function plot_model_sim
%%
clear
load lapse_test_fine.mat
name='Lapse';
fnum=1; pslot=12;
plotit(fnum, name, mwr, nR, nL, t_param, likey,pslot,wr)


%% gain
clear
load gain_test_fine.mat
fnum=2; pslot=11;
name='Gain';
plotit(fnum, name, mwr, nR, nL, t_param, likey,pslot,wr)

%% Noise

clear
load noise_test_fine.mat
% load noise_test.mat
fnum=3; pslot=10;
name='Noise';
plotit(fnum, name, mwr, nR, nL, t_param, likey,pslot,wr)


clear
load bias_test_fine.mat
% load noise_test.mat
fnum=4; pslot=8;
name='Bias';
plotit(fnum, name, mwr(4:25), nR, nL, t_param(4:25), likey(4:25),pslot,wr)

function plotit(fnum, name, mwr, nR, nL, t_param, likey,pslot,wr)

figure(fnum); clf
ax=axes('Position',[0.1 0.1 0.3 0.5]) ;set(ax,'NextPlot','add');
ylim(ax,[0 1]);
for px=1:numel(mwr)
    X(px)=t_param{px}(pslot);
    LL(px)=sum(log(likey{px}));
end

LL(isinf(LL))=nan;
[m,mi]=min(-LL);
[worst]=max(-LL);

fprintf('%s = %.3g, -LL=%.1f, fit/trial %.4g\n',name,X(mi),-LL(mi),mean(likey{mi}));
tstr=sprintf('Best %s = %.4g\n-LL=%.1f\nfit/trial %.4g',name,X(mi),-LL(mi),exp(LL(mi)/numel(likey{mi})));


for px=1:numel(mwr)
    
    [x,y,e]=binned(nR'-nL', mwr{px}','n_bins',15);
    clr=((-LL(px)-m)/(worst-m)/2)^.5;
    if ~isnan(clr)
        plot(ax,x,y,'Marker','none','Color',[clr 0 1-clr]);
    end
end

[x,y,e]=binned(nR'-nL', wr,'n_bins',7);
eh=errorplot(ax,x,y,e,'Marker','o','Color',[0.5 0.5 0.5]);
set(eh(2),'MarkerSize',3);

th=text(-32,0.95,tstr,'HorizontalAlignment','left','VerticalAlignment','top')

xlabel(ax,'#Contra - #Ipsi Clicks');
ylabel(ax,'P(Went Contra)');
set(ax,'Box','off','TickDir','out','TickLength',[0.025 0.1],'XLim',[-37 37],'YLim',[0 1])


ax2=axes('Position',[0.55 0.1 0.3 0.5]); set(ax2,'NextPlot','add');
set(ax2,'Box','off','TickDir','out','TickLength',[0.025 0.1])
buf=0.98;

ybuf=abs(m-worst)/2;

min_y=m-ybuf
max_y=max(-LL)+ybuf;
min_x=min(X);
max_x=max(X);


plot(ax2,X,-LL,'k.');
PP=spline(X,-LL);
PV=ppval(PP,min_x:0.001:max_x);

plot(min_x:0.001:max_x,PV,'k-');
set(ax2,'XLim',[min_x max_x],'YLim',[min_y max_y]);

xlabel(ax2,[name ' parameter']);
ylabel(ax2,'Neg. Log Likelihood')

title(name);
set(fnum,'PaperPosition',[1 1 6 4])
saveas(fnum,sprintf('%s_fine.eps',lower(name)),'epsc2');

