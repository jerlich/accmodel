
function savef(fig_h,savename,varargin)
% savef(fig_h,savename,varargin)
%
% A wrapper for saveas that formats figures for print.
%
% Usage:
% savef(gcf,'foo.pdf','PaperSize',[2 1]);
% Will make a pdf that is 2" wide and 1" high.
% 
% PaperSize=[2 2];
% PaperPosition=[0 0 PaperSize];
%


if nargin==0
   fig_h=gcf;
end

if nargin==1
   savename=['saved_fig_' datestr(now,'ddhhmm') '.pdf'];
   saveext='pdf';
else
   pind=find(savename=='.',1,'last');
   saveext=savename(pind+1:end);
end

PaperSize=[2 2];
PaperPosition=[];

overridedefaults(who,varargin);

if isempty(PaperPosition)
PaperPosition=[0 0 PaperSize(1) PaperSize(2)];
end

if strcmpi(saveext,'eps')
    saveext='epsc2';
end

set(gcf,'PaperSize',PaperSize,'PaperPosition',PaperPosition);

saveas(fig_h,savename,saveext);



   




   
