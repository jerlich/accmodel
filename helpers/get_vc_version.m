function [cvver,varargout]=get_vc_version(fulln)
% [ver,varargout]=get_vc_version(fullpathandfilename)
% returns the current version of an mfile
% This can be used when saving data to .mat files or the database, so that
% you know what code was used to create that data.
% Output
% ver   0 if file is locally modified
%      -1 if file is not found or not in versioning system
%       otherwise the version
%
% If no inputs are passed the function will get the fullname of the caller


if nargin==0
    S=dbstack;
    fulln=which(S(2).name);
end

olddir=pwd;
cvver=-1;
hre='File not found';
try
    
        lastslash=find(fulln==filesep,1,'last');
        pathstr=fulln(1:lastslash);
        filestr=fulln((lastslash+1):end);
        cd(pathstr);
        % Try cvs first
        cvsdir=dir('CVS');
        if ~isempty(cvsdir)
            [hre,cvver]=parse_cvsoutput(filestr);
        else
            
            [hre,cvver]=parse_svnoutput(filestr);
        end
        
    
    
catch le
    showerror(le);
    hre=le.message;
    cvver=-1;
    
end

cd(olddir);

if nargout==2
    varargout{1}=hre;
end

function [s,v]=parse_cvsoutput(filename)
[s,r]=system(['cvs stat ' filename]);

[s1,s2]=regexp(r,'Status:');
eol=find(r(s2:end)==10,1,'first');
s=strtrim(r(s2+1:(s2+eol)));



switch s
    case 'Locally Modified'
        v=0;
    case 'Unknown'
        v=-1;
    case 'Up-to-date'
        [v1,v2]=regexp(r,'Working revision:'); %#ok<*ASGLU>
        eol=find(r(v2:end)==10,1,'first');
        v=strtrim(r(v2+1:(v2+eol-1)));
        
end


function [s,v]=parse_svnoutput(filename)
[s,r]=system(['svn info ' filename]);
[s1,s2]=regexp(r,'Revision:');
if isempty(s1)
    s=r;
    v=-1;
else
    [v1,v2]=regexp(r,'Revision:'); %#ok<*ASGLU>
    eol=find(r(v2:end)==10,1,'first');
    v=strtrim(r(v2+1:(v2+eol-1)));
    
    [ss,rr]=system(['svn stat ' filename]);
    if ~isempty(rr)
        v=0;
        s='File is locally modified';
    end
end
    
    
    
    
    
