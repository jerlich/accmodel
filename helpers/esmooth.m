%esmooth  [S] = esmooth(S, sd)   smooth with care over edge effects
%
% This function smooths each row of S using a Gaussian with
% standard deviation sd, where sd is measured in samples. The
% function takes proper care of edge effects.
%
% If S is a column vector, the vector is smoothed and then also
% returned as a column vector.
%

% Carlos Brody June 1998

function [NewS] = esmooth(S, sd)

if iscolumn(S), was_column = 1; S = S';
else was_column = 0; end;

colsS = cols(S);


if colsS < 11*sd,  % Why 9x - think the 31,3 bug should be caught here
    NewS = space_smooth(S, sd);
    if was_column, NewS = NewS'; end;
    return;
end;

N = 2.^nextpow2(colsS);
x = -N/2 : N/2-1;
smoother = exp(-x.*x/(2.*sd.*sd));
smoother = smoother./sum(smoother);
fsmoother = fft(fftshift(smoother));
csum = cumsum(smoother);
smzero = N/2 + 1;

tsd = floor(3*sd);

NewS = zeros(size(S));

if tsd > 0,
    ML = zeros(tsd, 2*tsd);
    for r=1:tsd,
        
        ML(r,:) = smoother(smzero-r+1:smzero-r+2*tsd);
        ML(r,:) = ML(r,:)/sum(ML(r,:));
    end;
end;

padding = zeros(1, 2.^nextpow2(colsS)-colsS);

parfor i=1:rows(S),
    
    
    
    news = real(ifft(fft([S(i,:) padding]).*fsmoother));
    
    news = news(1:colsS); s = S(i,:);
    
    news(1:tsd)            = (ML*s(1:2*tsd)')';
    news(end:-1:end-tsd+1) = (ML*s(end:-1:end-2*tsd+1)')';
    NewS(i,:) = news;
end;

if was_column, NewS = NewS'; end;
return;



% -------------------------

function [NewS] = space_smooth(S, sd)

x = -(ceil(6*sd)+1):(ceil(6*sd)+1);
smoother = exp(-x.*x/(2.*sd.*sd));
smoother = smoother./sum(smoother);
csum = cumsum(smoother);
smzero = ceil(6*sd) + 2;
smlen = length(smoother);

len = cols(S);
NewS = zeros(size(S));

tsd = floor(4*sd);


sz = size(S,2);

for i=1:rows(S),
    
    s = S(i,:);
    
    news = zeros(1,sz);
    
    if tsd > 1,
        for k=1:numel(s),
            news(k) = sum(s(max(1,k-tsd):min(len,k+tsd)) .* ...
                smoother((smzero-tsd) + max(1-(k-tsd),0) ...
                : smzero+tsd + min(len-(k+tsd),0)));
            
            news(k) = news(k) ./ ...
                (csum(smzero+tsd + min(len-(k+tsd),0)) - ...
                csum((smzero-tsd) + max(1-(k-tsd),0) - 1));
        end;
    end;
    
    NewS(i,:) = news;
end;

return;
