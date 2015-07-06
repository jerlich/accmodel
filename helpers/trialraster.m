function [ys, xs] = trialraster(ref, trialidx, trialvals, pre, post,prewin,postwin)
% [y, x] = trialraster(ref, trialidx, trialsvals, pre, post)
%
% returns a raster of some trial data relative to ref trials




if length(trialidx) ~= length(trialvals),
    y = [];
    x = [];
      warning('trial index must be the same length as trial values');
    return;
end

max_pre=max(pre);
max_post=max(post);

xs = prewin:postwin;
ys = nan(numel(ref), numel(xs));

for i = 1:numel(ref),
    starting = ref(i) - pre(i);
    ending   = ref(i) + post(i);
    
    if isnan(starting) || isnan(ending)
        times=[];
        s=[];
    else
    [trials s] = qbetween2(trialidx, starting, ending);
    end
    if isempty(trials),
        ys(i,:) = NaN * ones(size(xs));
    else
        trials = trials - ref(i);
        vals  = trialvals(s(1):s(2));
		try
	        ys(i,max_pre-pre(i)+1:post(i)+max_pre) = vals(:)'; 
		
		catch
			ys(i,:) = NaN * ones(size(xs));
            fprintf(2,'Probably something wrong in trialraster\n');
		end;
    end
end