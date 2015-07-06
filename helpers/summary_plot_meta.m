function fh=summary_plot_meta(cellid, varargin)

pairs={'config'       '';...    
    
    }; parseargs(varargin, pairs);


switch config
    case 'PBups' 
        % aligned to stimulus onset
%         fh=summary_plot2(cellid,'align_on','PBUPS','by_bup_diff', 2,'pre',0.2,'post',0.5,'mask_post_cpoke1',1,'correct_only',0,...
%             'plot_hv',1,'plot_isi',0,'plot_wave',0,'num_plots',2,'this_plot_num',1);
        fh=summary_plot2(cellid,'align_on','PBUPS','pre',0.2,'post',0.8,'mask_post_cpoke1',1,'by_choice',1,'correct_only',0,...
            'plot_hv',1,'plot_isi',0,'plot_wave',0,'num_plots',2,'this_plot_num',1);
        
        % aligned to end of stimulus
%         summary_plot2(cellid,'align_on','PBUPS','by_bup_diff', 2,'pre',2.0,'post',2.0,'plot_hv',1,'plot_isi',1,'correct_only',1,...
%             'plot_wave',1,'fh',fh,'num_plots',2,'this_plot_num',2);
%         summary_plot2(cellid,'pre',2.0,'post',2.0,'plot_hv',1,'plot_isi',1,'plot_wave',1,'fh',fh,'num_plots',2,'this_plot_num',2);
        summary_plot2(cellid,'pre',2.0,'post',2.0,'by_choice',1,'plot_hv',1,'plot_isi',1,'plot_wave',1,'fh',fh,'num_plots',2,'this_plot_num',2);
    otherwise
        fh = summary_plot(cellid);
end