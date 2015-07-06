

% function [y x W] = warpfilter(ref, ts, kernel, varargin)
%
% produces smoothed and warped single-trial peth's in a window [-pre post] from each
% element of the reference (ref) matrix.  If there is only one column of
% ref than no warping takes place.
% all times are assumed to be in seconds
%
%
% Inputs:
%	ref					matrix of times of reference events, it should be
%                       size [Trials, refs].  Each element of a row should
%                       be greater than the previous element of that row.
%	ts					time stamps of spikes
%	kernel				kernel used for filtering
%
% Varargin:
%	kernel_bin_size		size of a kernel bin (sec)
%	pre					sec before each ref event in window
%	post				sec after each ref event in window
%
% Outputs:
%	x					a vector of time steps, [-pre:kernel_bin_size:post]
%	y					a M x N matrix, where M is the length of ref and N
%						is the length of x.
%						each row is the smoothed spike raster for one trial
%   W                   The times of the warped spikes.  If you take W and
%                       use ref(:,1) and W as the inputs to a function that takes
%                       ref and spike times (like rasterplot) you should
%                       produce a figure that matches y.
%
% if there are nans in ref, that row of y is nans
%
% J. Erlich 2011/8/1

function [y x W]=warpfilter(refs,ts,kernel,varargin)


pairs = {'kernel_bin_size'			5e-4	; ...
    'pre'						2		; ...
    'post'						3		; ...
    }; parseargs(varargin, pairs);



[ntrials, nrefs]=size(refs);

kernel = kernel/sum(kernel)/kernel_bin_size; % normalize
offset = ceil(length(kernel)/2); % The offset compensates for the way that filter deals with edges.
buffered_pre=pre+offset*kernel_bin_size;
x = -buffered_pre:kernel_bin_size:post;
y = zeros(ntrials, numel(x)-1);
ts=ts(:)';   % make ts a row vector;

% compute the median time of each of the events relative to the first
% event.
mutau=zeros(1,nrefs);
for rx=2:nrefs
    mutau(rx)=nanmedian(refs(:,rx)-refs(:,1));
end

W=[];  % Initialize the warped timestamps vector
for tx=1:ntrials
    if any(isnan(refs(tx,:)))
        % if any refs on this trial are nan, we can't warp.
        y(tx,:)=y(tx,:)+nan;
    else
        start = refs(tx,1) - buffered_pre;
        fin   = refs(tx,nrefs) + post;
        spks = qbetween(ts, start, fin); % spike times relative to ref

        if ~isempty(spks), % if there are no spikes on this trial we don't do anything
            
            % get the spikes that are before the first ref event and align
            % them with no warping.
            pre_spks=spks(spks<refs(tx,1));
            warped=pre_spks-refs(tx,1);
            
            for rx=2:nrefs
                % Get the spikes between ref(n-1) and ref(n) relative to
                % ref(n-1)
                tspks=qbetween(spks,refs(tx,rx-1),refs(tx,rx))-refs(tx,rx-1);
                % Normalize the times so that they go between mutau(n-1)
                % and mutau(n).
                tspks=tspks/(refs(tx,rx)-refs(tx,rx-1))*(mutau(rx)-mutau(rx-1))+mutau(rx-1);
                warped=[warped tspks];
            end
            % Get the spikes after the last ref event and align them with
            % no warping.
            post_spks=spks(spks>refs(tx,end))-refs(tx,end)+mutau(end);
            warped=[warped, post_spks];
            W=[W warped+refs(tx,1);];
            
            % Bin the warped spikes.
            ty = histc(warped,x);
            y(tx, :) = ty(1:end-1);
        end
    end
end

y=[y zeros(ntrials, offset)];  % pad with extra zeros

% Apply the smoothing kernel using the builtin matlab filter function.
y = filter(kernel, 1, y, [], 2);
y = y(:, 2*offset:end-1); % trim extra columns
x=x(offset+1:end-1);
