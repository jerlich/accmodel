function [L R] = make_adapted_cat_clicks(leftbups, rightbups, phi, tau_phi, varargin)

pairs = {
    'cross_side_suppression'  0; ... % clicks on different sides closer than this duration (in sec) cancel and annihilate
}; parseargs(varargin, pairs);

if abs(phi - 1) > eps,
    lefts  = [leftbups(:)  -ones(numel(leftbups),1)];
    rights = [rightbups(:) +ones(numel(rightbups),1)];
    allbups = sortrows([lefts; rights])'; % one bup in each col, second row has side bup was on
    ici = diff(allbups(1,:));
    
    adapted = ones(1, size(allbups,2));
    
    
    
    for i = 2:size(allbups,2),
        if ici(i-1) <= cross_side_suppression,
%             adapted(i-1) = 0.5*adapted(i-1);
%             adapted(i) = 0.5*adapted(i);
            adapted(i-1) = 0;
            adapted(i) = 0;
        else
            last = tau_phi * log(1 - adapted(i-1)*phi);
            adapted(i) = 1 - exp((-ici(i-1) + last)/tau_phi);
        end;
    end;
end;

adapted = real(adapted);

L = adapted(allbups(2,:)==-1);
R = adapted(allbups(2,:)==+1);

