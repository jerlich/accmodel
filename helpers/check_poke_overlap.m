


function [pcout,pgocue,stime,cout]=check_poke_overlap(peh,varargin)
%[precout,pregocue]=check_poke_overlap(peh)
% returns 2 binary vectors that indicated whether there was an spoke before
% cpoke1(end,end) or before cout
%

pairs={'NIC_state', 'cpoke1'};
parseargs(varargin, pairs);

pcout=zeros(numel(peh),1)+nan;
pgocue=pcout;
stime=pcout;
cout=pcout;

cin=extract_event(peh,[NIC_state '(end,1)']);
gos=extract_event(peh,[NIC_state '(end,end)']);


for tx=1:numel(peh)
	if ~isnan(cin(tx)) 
		
		rp=peh(tx).pokes.R;
		lp=peh(tx).pokes.L;
		cp=peh(tx).pokes.C(:,2);
		
		LRp=sort([rp(:,1); lp(:,1)]);
		coi=find(cp>gos(tx),1,'first');
		sti=find(LRp>cin(tx),1,'first');
		if ~isempty(coi) && ~isempty(sti)
			stime(tx)=LRp(sti);
			cout(tx)=cp(coi);
			pcout(tx)=stime(tx)<cout(tx);
			pgocue(tx)=stime(tx)<gos(tx);
		end
		
	end
end