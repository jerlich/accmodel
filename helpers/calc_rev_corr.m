function G=calc_rev_corr(rawdata,varargin)

%% Inputs
%
% rawdata         the output of package_pbups_data/or scrach code.  If rawdata contains a
%                 field probR, the function assumes that that is a model
%                 prediction and generates the rev_corr in a weighted
%                 fashion.  E.g. if the model predicts that the subject
%                 will go left 72% of the time, then that trial is added to
%                 the left trials with a weight of 0.72 and to right trials
%                 with a weight of 0.28.  If probR does not exist then the
%                 function assumes that is a real data without weighting.
%
% Output 
% G             
%
% Optional Imputs (and their default values)
%
% inc = [];
%       A vector of logicals the same length as rawdata.  If this is not
%       empty it overrides time_range and gamma_range
% bin = 0.005;
%       The bin size.
% time_range = [0.4 1];
%       Trials in this time range will be included
% gamma_range = [0 4];
%       Trials in this gamma range will be included
% plot_it = 1;
%       If this is true, then the rev corrs are plotted, forward on the
%       left and backward on the right
% krn_sd = 0.04;
%       The standard deviation of the smoothing kernel
% use_esmooth = true;
%       If this is true then the esmooth function is used.  Otherwise the
%       built-in smoothing in spike_filter (which uses filter)
%       esmooth is better with edge conditions, but filter is faster
% as_data = false;
%       If this is set to true then even if probR is a field, it is ignored
%       and the rev corrs are generated using pokedR.  This allows you to
%       have one rawdata structure with both pokedR and probR.


%% Defaults

inc = [];
bin = 0.01;
time_range = [1 4];
gamma_range = [0 3];
plot_it = 1;
krn_sd =0.04;
use_esmooth = true;
as_data = false;

G=overridedefaults(who,varargin);

if use_esmooth
    krn=1;
    smooth_krn=krn_sd/bin;
    
else
    krn = normpdf(-.5:bin:.5,0,krn_sd);
end
left_clicks=[];
right_clicks=[];



gamma = cell2mat({rawdata.gamma});
T = cell2mat({rawdata.T});

if isempty(inc)
    inc = abs(gamma)>= gamma_range(1) & abs(gamma)<=gamma_range(2) & T>=time_range(1) & T<=time_range(2);
end

rawdata=rawdata(inc);
upto = 1*max(cell2mat({rawdata.T})); % play with which data durations 
wr = cell2mat({rawdata.pokedR});
% concatenate the clicks=
for tx=1:numel(rawdata)
    left_clicks = [left_clicks; rawdata(tx).leftbups(:)+10*(tx-1)];
    right_clicks = [right_clicks; rawdata(tx).rightbups(:)+10*(tx-1)];
    start_time(tx)=10*(tx-1);
    end_time(tx)=10*(tx-1)+rawdata(tx).T;
end

gamma = [rawdata.gamma];  % this rawdata has been filtered by gamma_range and time_range.
tt = unique(gamma);

G(1).gamma=gamma;

for ex=1:2
    if ex==1
        [fly, flx] = spike_filter(start_time,left_clicks,krn,'pre',0,'post',upto,'kernel_bin_size',bin);
        [fry, ~] = spike_filter(start_time,right_clicks,krn,'pre',0,'post',upto,'kernel_bin_size',bin);
        fy = fry - fly;
        fy = maskraster(flx,fy,-inf,end_time-start_time);
        
    else
        [fly, flx] = spike_filter(end_time,left_clicks,krn,'pre',upto,'post',0,'kernel_bin_size',bin);
        [fry, ~] = spike_filter(end_time,right_clicks,krn,'pre',upto,'post',0,'kernel_bin_size',bin);
        fy = fry - fly;
        fy = maskraster(flx,fy,-end_time+start_time,+inf); 
    end
    if use_esmooth
       fy = do_smooth(fy,smooth_krn);
    end
    
    for ttx = 1:numel(tt)
        mu(ttx,:) = nanmedian(fy(gamma==tt(ttx),:)); % multiply by zeros to get rid of the mu if want to check each gamma alone
    end
    
    for tx=1:size(fy,1)
        nfy(tx,:) = fy(tx,:) - mu(gamma(tx)==tt,:);
    end
    
    if isfield(rawdata,'probR') && ~as_data
        probR = cell2mat({rawdata.probR});
        G(ex).rev_corr_r = bsxfun(@times,nfy,probR')/mean(probR);
        G(ex).rev_corr_l = bsxfun(@times,nfy,1-probR')/mean(1-probR);
    else
        G(ex).rev_corr_r = nfy(wr>0.5,:);
        G(ex).rev_corr_l = nfy(wr<0.5,:);
    end 
    G(ex).tax = flx;

end

if plot_it
    do_plot(G)
end

function do_plot(G)
figure(148); clf
ax(1) = axes('Position',[0.15 0.3 0.3 0.4]);
ax(2) = axes('Position',[0.55 0.3 0.3 0.4]);
xlab={'Time from Stim Onset' , 'Time from Stim Offset'};
for x=1:2
    mur = nanmean(G(x).rev_corr_r);
    ser = nanstderr(G(x).rev_corr_r);
    mul = nanmean(G(x).rev_corr_l);
    sel = nanstderr(G(x).rev_corr_l);
    
    he(2*x-1,:)=errorplot(ax(x),G(x).tax,mur,ser,'Color','r');
    he(2*x,:)  =errorplot(ax(x),G(x).tax,mul,sel,'Color','g');
    
    xlabel(ax(x),xlab{x});
    if x==1
        xlim(ax(x),[0 G(x).tax(end)*.85])
        ylabel(ax(x),'Excess Click Rate (Hz)');
    else
        xlim(ax(x),[G(x).tax(1)*.85 0])
    end
    ylim(ax(x),[-4 4])
    
end



% end do_plot

%%

function sy = do_smooth(y,krn)

sy = nan+y;

for rx=1:size(y,1)
    ty=y(rx,:);
    goody = ~isnan(ty);
    nty = ty(goody);
    sty = esmooth(nty,krn);
    sy(rx,goody) = sty;
end




