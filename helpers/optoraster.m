function optoraster(cellid,varargin)

S=[];
C=[];
if ispc
   save_path='c:\';
else
   save_path='~/';
end
save_flag=true;

overridedefaults(who,varargin);

if isempty(C)
   C=get_celldata(cellid);
end

if isempty(S)
   S=get_sessdata(C.sessid);
end




for tx=1:numel(C.cellid)

   sx=find(C.sessid(tx)==S.sessid);
   
peh=S.peh{sx};
pd=S.pd{sx};
%%
stim=~isnan(extract_waves(peh,'stimulator_wave1'));
right=pd.sides=='r';
ref=extract_event(peh,'nicstim(1,1)');

cnd=10*right+stim;
%%
ax=[];
   figure(100+tx);clf
[ax(tx,:),D{tx}]=exampleraster(ref,C.ts{tx},'cnd',cnd,'pre',1.5,'post',4,'corner',[0.2 0.2],'total_height',0.5,'psth_height',0.15,'renderer','painters','errorbars',0,'krn',0.15,'clrs',{'b' 'b' 'r' 'r'},'x_label','Time from Sound Onset (s)');


% Make the PSTH for stim dashed
ch=get(ax(tx,end),'Children');
set(ch([1 3]),'LineStyle','--');

% Make the 0 line dashed

ch=get(ax(tx,2),'Children');
set(ch([1]),'LineStyle','-.','LineWidth',2);
ch=get(ax(tx,4),'Children');
set(ch([1]),'LineStyle','-.','LineWidth',2);



set(100+tx,'PaperPosition',[0.25 0.25 5 5]);

if save_flag
saveas(100+tx,sprintf('%s%s%d_opto.eps',save_path, filesep, C.cellid(tx)),'epsc2');
end

end


