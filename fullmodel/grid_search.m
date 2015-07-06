function [D]=grid_search(x_bf,rawdata,varargin)
% D = grid_search(x_bf, rawdata, [optional inputs]);
% Inputs:
% x_bf      The best_fit_params from run_fmincon35
% rawdata   The data that was used to generate the best fit parameters
%
% Optional Inputs:
% do_param  A 1x9 logical vector which chooses which cross sections to
%           sample.  At least 2 elements must be true.  The code will
%           generate n choose 2 cross sections where n is the # of true
%           elements in the vector.  
%           lambda sig2_a sig2_s sig2_i bounds phi tau_phi bias lapse
%           Default: [1   1 1    0   0  0   0 0 1]
% samp_mech The sampling method.  Default norm.  
%           'norm': samples using norm_sd around the best fits
%           'uni' : samples from a uniform distribution with range
%           uni_range around the best fits.
%           'grid': samples in a grid from the best fits using step_size
% step_size [.2 .2 5   .5   2  .1 .1 1 .05];
% uni_range [10 30 300  20  20  2  1 10 1];
% norm_sd   [3   5 100  20  20  1  1 5  1];
% n_bins    At max (2*n_bins+1)^2 points will be in each cross section.
%           Code generates points and then discards parameters that are not
%           meaningful. e.g. negative variances.
% 
% Outputs:
% D             A struct with 2 fields
% D.p_search    A n x 9 matrix with all the parameter sets sampled
% D.LL          An n x 1 vector with the log likelihood of the model at 
%               each of the searched parameters sets 

%%
% lambda sig2_a sig2_s sig2_i bounds phi tau_phi bias lapse
if numel(x_bf)==9
do_param =  [1   1 1    0   0  0   0 0 1]==1;
step_size = [.2 .2 5   .5   2  .1 .1 1 .05];
uni_range = [10 30 300  20  20  2  1 10 1];
norm_sd =   [3   5 100  20  20  1  1 5  1];
elseif numel(x_bf)==12
    do_param =  [1   1 1    0   0  0   0 1 1 1 1 1]==1;
    step_size = [.2 .2 5   .5   2  .1 .1 1 .05  5 0.1 .05];
    uni_range = [3  10 50  20  20  2  1 10 1  300 1 1];
    norm_sd =   [3   5 100  20  20  1  1 5  1 100 .2  1];
end


n_bins = 10;
samp_mech = 'norm';

overridedefaults(who,varargin);
if isnumeric(do_param)
    do_param = do_param==1;
end

p_search = repmat(x_bf,(2*n_bins+1)^2*nchoosek(sum(do_param),2),1);

p_todo=find(do_param);

p_pairs=nchoosek(p_todo,2);
ind = 1;
switch samp_mech,
case 'grid',
    % A straight up grid search
    for ppx=1:size(p_pairs,1)
        p1 = p_pairs(ppx,1);
        p2 = p_pairs(ppx,2);
        ss1 = step_size(p1);
        ss2 = step_size(p2);
        mu1 = x_bf(p1); mu2 = x_bf(p2);
        vec1 = (mu1-n_bins*ss1):ss1:(mu1+n_bins*ss1); 
        vec2 = (mu2-n_bins*ss2):ss2:(mu2+n_bins*ss2); 
        for v1x = 1:numel(vec1)
            for v2x = 1:numel(vec2)
                p_search(ind,p1)=vec1(v1x);
                p_search(ind,p2)=vec2(v2x);
                ind=ind+1;
            end
        end   
    end


case 'uni',
    % Sample randomly in a given range
    for ppx=1:size(p_pairs,1)
        p1 = p_pairs(ppx,1);
        p2 = p_pairs(ppx,2);
        mu1 = x_bf(p1); mu2 = x_bf(p2);
        tot_points = (2*n_bins+1)^2;
        vec1 = ((-0.5+rand(1,tot_points))*2*uni_range(p1))+mu1; 
        vec2 = ((-0.5+rand(1,tot_points))*2*uni_range(p2))+mu2; 

        for np = 1:numel(vec1)
                p_search(ind,p1)=vec1(np);
                p_search(ind,p2)=vec2(np);
                ind=ind+1;
            
        end   
    end

case 'norm',
    % Sample mostly near the best fit params and sparesly farther away.
    for ppx=1:size(p_pairs,1)
        p1 = p_pairs(ppx,1);
        p2 = p_pairs(ppx,2);
        mu1 = x_bf(p1); mu2 = x_bf(p2);
        tot_points = (2*n_bins+1)^2;
        
        vec1 = ((randn(1,tot_points))*norm_sd(p1))+mu1; 
        vec2 = ((randn(1,tot_points))*norm_sd(p2))+mu2; 

        for np = 1:numel(vec1)
                p_search(ind,p1)=vec1(np);
                p_search(ind,p2)=vec2(np);
                ind=ind+1;
        end   
    end
end


p_search = unique(p_search,'rows');

p_search=validate_psearch(p_search);

ll = @(x)(ll_all_trials_opt(x,rawdata)); 

LL=nan(size(p_search,1),1);

for px=1:size(p_search,1)
    fprintf('%.2f ',p_search(px,:))
    fprintf('searching %d/%d \n',px,size(p_search,1));

    LL(px) = ll(p_search(px,:));
end

D.LL=LL;
D.p_search=p_search;
D.do_param=do_param;


function p=validate_psearch(p)
    % Make sure that the parameter values are meaningful.
if size(p,2)==9
    bad = p(:,2)<0 | p(:,3)<0 | p(:,4)<0 | p(:,5)<0 | p(:,6)<0 | p(:,7)<0 | p(:,9)<0 | p(:,9)>1;
    % lambda sig2_a sig2_s sig2_i bounds phi tau_phi bias lapse
elseif size(p,2)==12
    
    bad = p(:,2)<0 | p(:,3)<0 | p(:,4)<0 | p(:,5)<0 | p(:,6)<0 | p(:,7)<0 | p(:,9)<0 | p(:,9)>1 | p(:,12)<0 | p(:,12)>2 | p(:,11)<0 | p(:,10)<0;
end
    
    p=p(~bad,:);
