function inds = find_closest(range, vector)
% flatten everything into row vectors
origshape = size(vector);
sz_r = size(range);
numTargets = sz_r(1) * sz_r(2);
numInputs = origshape(1) * origshape(2);

range = reshape(range, [1, numTargets]);
vector = reshape(vector, [1, numInputs]);

% find the difference between the target values and the actual values
% range is replicated across columns,
% vector is replicated across rows
diff = abs(repmat(range', 1, numInputs) - repmat(vector, numTargets, 1));

[~, inds] = min(diff);
inds = reshape(inds, origshape);
end