% script to process all best fit information synthesized from their various
% outputs

% expects all_chrono_rats, a struct array, to be populated by something
% like extract_fmincon_quadfit_out.m
% and has (at a miniumu) these fields for each rat:
%       .ratname
%       .fval
%       .x
%       .H
%       .Bnd
%           .Bnd.B
%           .Bnd.LL_B
%           .Bnd.LLgrad_B
%       .likey

data_dir = '~/ratter/PoissonClicks/BingModeling/data/';
filename = 'best_fits_all_trials.mat'; 

load([data_dir filename]);

%% unpack all_chrono_rats
ratnames = [];
for rati = 1:numel(all_chrono_rats),
    ratnames         = [ratnames; all_chrono_rats(rati).ratname];
    allfval(rati)    = all_chrono_rats(rati).fval;
    allfits(rati,:)  = all_chrono_rats(rati).x;
%     allcov(rati,:,:) = inv(all_chrono_rats(rati).H);
%     allci(rati,:)    = sqrt(diag(inv(all_chrono_rats(rati).H)));
%     alllikey{rati}   = all_chrono_rats(rati).likey;
    if ~isempty(all_chrono_rats(rati).Bnd),
        allB(rati)       = all_chrono_rats(rati).Bnd;
    end;
end;

%% process the B linescan
% what we're trying to do is parse out the B dimension sensibly, since it
% does not always have a minimum in LL
% epsilon = 2;
epsilon = 0.1;
figure;
for rati = 1:numel(allB),
    B        = allB(rati).B;
    LL_B     = allB(rati).LL_B;
    LLgrad_B = allB(rati).LLgrad_B;

    if isempty(B),
        fprintf('rat number %i, %s doesn'' have a B linescan.\n', rati, ratnames(rati,:));
    else
        x = allfits(rati,:);
        B_bf = x(4);
        fval = allfval(rati);

        if LL_B(1) > LL_B(end),
            mbf = min(LL_B);
            dip = abs(LL_B - mbf) < epsilon;
            
        else
            % what's the range over which the LL_B is within epsilon of fval?
            dip = abs(LL_B - fval) < epsilon;
            LL_B = -LL_B;
        end;
        
        LL_B = LL_B - min(LL_B);
        
        fprintf('rat number %i, %s\n', rati, ratnames(rati,:));
        display(dip);

        alldips(rati,:) = dip;
        
        subplot(4,5,rati); plot(B, LL_B, '.-'); title(ratnames(rati,:));
        set(gca, 'XLim', [0 32], 'YLim', [0 50]);
        
     
    end;
end;

%% save
data_dir = '~/ratter/PoissonClicks/BingModeling/data/';
save([data_dir 'best_fits_all_trials.mat'], 'all_chrono_rats', 'ratnames', 'allfits', ...
    'allH', 'allci', 'alllikey', 'allfval', 'alldips');