
function varargout=sdata(sqlstr, varargin)
%function bdata
% connection_id=bdata
% out=bdata(sqlstr)
% out=bdata('connect',host,user,pass)
% connects to sql server via mym and maintains a connection id.
%
% See also MYM
%
%
% This bdata is compatible with mym 1.36 which returns a single struct output with each requested column as a fieldname.
% earlier versions of mym (<=1.0.9) will not work with this.
%

persistent CON_ID
persistent IP_ADDR
writehost= writehost;
readhost = readhost;


if ~isempty(CON_ID) && usejava('jvm')
    cur_addr=java.net.InetAddress.getLocalHost.getHostAddress;
    if ~isequal(cur_addr, IP_ADDR)
        CON_ID=[];
        fprintf(2,'IP address changed, discarding previous connection\n');
    end
end

try
    
    if nargin==0
        varargout{1}=CON_ID;
        return;
    end
    
    % first check that the server is there
    if ~check_connection(readhost)
        warning('sdata:noserver','Server is not accessible')
        CON_ID=[];
        varargout{1}=-1;
        return;
    end
    
    action=lower(strtok(sqlstr,' '));
    
    switch action,
        
        %% connect
        case 'connect'
            if nargin<4
                muser='user';
                mpass='pass';
            else
                readhost=varargin{1};
                muser=varargin{2};
                mpass=varargin{3};
            end
            
            try
                
                % mym supports multiple simultaneous connections.
                % Passing -1 asks for the next available connection id.
                CON_ID(1)=mym(-1, 'open',writehost,muser, mpass);
                CON_ID(2)=mym(-1, 'open',readhost,muser, mpass);

                mym(CON_ID(1),'use bdata');
                mym(CON_ID(2),'use bdata');

                display('connected')
                varargout{1}=CON_ID;
                if usejava('jvm')
                    IP_ADDR = java.net.InetAddress.getLocalHost.getHostAddress;
                else
                    IP_ADDR=1;
                end
            catch
                varargout{1}=-1;
                showerror(lasterror);
                varargout{2}=lasterror;
            end
            
            
            %% close
        case 'close'
            if isempty(CON_ID)
                warning('sdata:nohandle','No handle to sdata server')
            else
                varargout{1}=mym(CON_ID(1), 'close');
                varargout{1}=mym(CON_ID(2), 'close');

            end
            CON_ID=[];
            IP_ADDR=[];
            %% status
        case 'status'
            varargout{1}=mym(CON_ID(1), 'status');
            %% sql
        case {'select','insert','show','explain','describe','call','update'}  % this is an sql statment.
            % by only allowing select and insert it means that properly
            % inserted data cannot be corrupted.
            
            if isempty(CON_ID)
                not_connected=1;
                fprintf(1,'First time connecting with bdata server since matlab start.\n');
            else
                not_connected=mym(CON_ID(1), 'status');
                if not_connected
                    warning('bdata:lostconnection','May have lost connection with mysql server');
                end
            end
            if not_connected
                [cid]=sdata('connect');
                
                % This prevents the code from looping endlessly if the server
                % isn't up.
                if cid==-1
                    warning('bdata:noserver','Failed to connect to bdata server');
                    return;
                end
                
            end
            
            
            %       mym can  be used like sprintf where place holders like "{S}"
            %       are placed in the sql statement and those are filled by a comma
            %       seperated list of variables.
            
            varlist='';
            if nargin>1
                vs=varargin;
%                 for vx=1:numel(vs)
%                     if (isscalar(vs{vx}) && isnan(vs{vx})) || isempty(vs{vx})
%                         vs{vx}='NULL';
%                     end
%                 end
                varlist=',vs{1}';
                for vx=2:numel(vs)
                    varlist=[varlist ', vs{' num2str(vx) '} '];
                end
            end
            
            %       When mym is used with no outputs it prints out a table.  However,
            %       we need to contruct a varargout string for the case where outputs
            %       are requested.
            %
            if nargout>0
                outstr='S=';
            else
                outstr='';
            end
            
            if isequal(lower(action),'insert') || isequal(lower(action),'update') || isequal(lower(action),'call')
                evalstr=[outstr 'mym(CON_ID(1),''' sqlstr ''''  varlist ');'];
            else
                evalstr=[outstr 'mym(CON_ID(2),''' sqlstr ''''  varlist ');'];
            end
            eval(evalstr)
            
            if nargout>0
                fn=fieldnames(S);
                for ox=1:nargout
                    varargout{ox}=S.(fn{ox});
                end
            end
            
            
        otherwise
            warning('Mysql:bdata',['The interface does not support the action: ' action])
            
    end
    
catch me
        showerror(me)

    fprintf(2,'ERROR in bdata: %s\n',evalstr);
end
