function [S, extra_args]=get_sessdata(varargin)
% [S, extra_args]=get_sessdata(varargin)
% [S, extra_args]=get_sessdata(sessid)
% [S, extra_args]=get_sessdata(ratname,experimenter, date)
% [S, extra_args]=get_sessdata(ratname,experimenter, daterange)
%
% A frontend to get data from the sessions table that does some nice input
% parsing.   Useful to use in other functions to avoid having to parse
% inputs.  If you pass it all the args from a parent function the leftover
% args are returned as extra_args.  For a good example of this see
% psychoplot_delori.m (in ExperPort/Analysis/SameDifferent)
%
% S by default contains fields: sessid, sessiondate,protocol_data and peh
%
% sessid can be a single sessid or a vector of sessids
% date should be of the form "YYYY-MM-DD" or a relative date like -5
% daterange should be a numeric vector in relative form like -10:-1 or a
% cell array of date string of the from "YYYY-MM-DD"
%
% pairs={'do_tracking' false;...
% 	   'fetch_peh' true...
% 	   };
%
%

if nargin==0 || isempty(varargin{1})
    S.sessid=[];
    S.pd={};
    S.peh={};
    S.ratname={};
    S.sessiondate={};
    S.protocol={};
    return;
end
    

if iscell(varargin{1})
	varargin=varargin{1};
	nargs=numel(varargin);
else
	nargs=nargin;
end

%% parse inputs
use_sessid=0;
if isnumeric(varargin{1})
	% Case 1, we've got a vector of sessids
	sessid=varargin{1};
    sessid=sessid(~isnan(sessid));
	use_sessid=1;
	[ratname, experimenter]=bdata('select ratname, experimenter from sessions where sessid="{S}"',sessid(1));
	varargin=varargin(2:end);
elseif nargs>=3
	% Case 2, we've got a ratname and experimenter
	ratname=varargin{1};
	experimenter=varargin{2};
	datein=varargin{3};
	if isnumeric(datein)
		%Case 2a, we've got relative dates (e.g. -10:0)
		for dx=1:numel(datein)
			dates{dx}=to_string_date(datein(dx));
		end
	elseif ischar(datein)
		%Case 2b, we've got a single date (e.g. '2009-05-01')
		dates{1}=datein;
	else
		%Case 2c, we've got a cell array of dates
		dates=datein;
	end
	% In future we might allow in extra parameters
	varargin=varargin(4:end);
else
	S=[];
	warning('Failed to parse inputs.');
	extra_args=varargin;
	return
end

extra_args=varargin;

pairs={'do_tracking' false;...
	   'fetch_peh' true...
	   };
parseargs(extra_args, pairs,[],1);


%% get data from sql

if ~use_sessid
	% If we are not in Case 1 (see above)
	% then transform the cell array of strings into a long comma separated
	% string.
	datestr='';
	for dx=1:numel(dates)
		datestr=[datestr ',"'  dates{dx}  '"'];
	end
	% Use the datestr for a select ... where sessiondate in (datestr) type sql command to get all the relevant sessions. 
	[sessid]=bdata(['select sessid from sessions where ratname regexp "{S}" and experimenter="{S}" and sessiondate in (' datestr(2:end) ') order by sessiondate'],ratname, experimenter);
end



	% We have a list of sessids.  Transform that into a comman seperated string
	sessstr='';
	for sx=1:numel(sessid)
		sessstr=[sessstr, ',' num2str(sessid(sx))];
	end
	if fetch_peh
	[S.sessiondate,S.pd,S.sessid, S.protocol,S.peh,S.ratname]=bdata(['select sessiondate, protocol_data,s.sessid, protocol, peh,s.ratname from sessions s,parsed_events p where s.sessid in (' sessstr(2:end) ') and s.sessid=p.sessid order by sessiondate']);
	else
	[S.sessiondate,S.pd,S.sessid, S.protocol,S.ratname]=bdata(['select sessiondate, protocol_data,s.sessid, protocol,s.ratname from sessions s where s.sessid in (' sessstr(2:end) ') order by sessiondate']);
	end	
	
if do_tracking
    S.a=cell(numel(S.sessid),1);
    
    [T.sessid, T.ts, T.theta]=bdata(['select sessid, ts, theta from tracking where sessid in (' sessstr(2:end) ')']);
    
    for sx=1:numel(S.sessid)
        
       tx=find(T.sessid==S.sessid(sx));
	   if ~isempty(tx)
       a.ts=T.ts{tx};
       a.theta=T.theta{tx};
       
       S.a{sx}=a(:);
	   end
    end
     


end


    


