function z=fill_in_nan(z,t)

if nargin==1 
    t=1;
end

for tx=1:t
[ni,nj]=find(isnan(z));

for nx=1:numel(ni)
    z(ni(nx),nj(nx))=avgnear(z,ni(nx),nj(nx));
end
end


function y=avgnear(z,ni,nj)

[rows,cols]=size(z);
sump=0;
nump=0;
for x = [-1 0 1]
    for y= [-1 0 1]
        newx=ni+x;
        newy=nj+y;
        if newx>0 && newy>0 && newx<=rows && newy<=cols && ~isnan(z(newx,newy))
            sump=sump+z(newx,newy);
            nump=nump+1;
        end
    end
end

y=nanmean(sump)/nump;

            