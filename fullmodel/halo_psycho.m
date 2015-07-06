
load chuck_data
% This was already done: 
%    rawdata=flipit(ratdata)
%    allrawdata=rawdata;
%   C=allrawdata(stim==0);

stim=cell2mat({rawdata.stim});
rawdata=allrawdata(stim==1);
dc=cell2mat({rawdata.Delta});
wr=cell2mat({rawdata.pokedR});
cdc=cell2mat({C.Delta});
cwr=cell2mat({C.pokedR});


%%


figure(1);clf
ax=axes;

[x,y,e]=binned(dc,wr,'n_bins',6);
errorplot(ax,x,y,e,'Color',[1 0.5 0])
[cx,cy,ce]=binned(cdc,cwr,'n_bins',6);
errorplot(gca,cx,cy,ce,'Color','k')

xhairs(ax,':k',0,0.5)
