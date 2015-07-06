% run this on the cluster, in the directory containing sub-directories with
% rats' data

total_rate = 40;

for rati = 1:numel(all_chrono_rats),
    ratname = all_chrono_rats(rati).ratname;
    x_bf = all_chrono_rats(rati).x;
    
    load([ratname '/chrono_' ratname '_rawdata.mat']);
    
    oneclick_trials = find(abs(avgdata.Delta)==1);
    
    fprintf('running one click trials for %s... \n', ratname);
    
    [~, ~, likey] = ll_model_wrapper35(x_bf, rawdata(oneclick_trials), ...
                      'track_history', 0, 'for_fmin', 1, ...
                      'do_param', ones(1,9), 'total_rate', total_rate);
                  
    all_chrono_rats(rati).oneclick_trials = oneclick_trials;
    all_chrono_rats(rati).oneclick_likey = likey;
    
    
    save('best_fits_v35_g2c.mat', 'all_chrono_rats');
end;
    