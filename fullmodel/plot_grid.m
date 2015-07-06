function E=plot_grid(fname,varargin)
% plot_grid(fname, [optional inputs]);
%
% Inputs:
% fname     The filename containing the D struct from grid_search
%
% Optional Inputs:
% do_fit    If true, the data will be fit with a quadratic.
%
% Outputs:
% E         A struct with the fits of a 2D quadratic for each of the cross
%           sections

do_fit=false;
thresh = 0.8;
norm_plot=false;
n_x_bins=18;
n_y_bins=18;
plot_type='heat';

overridedefaults(who,varargin)

if iscell(fname)
  x_bf=fname{1};
  D = fname{2};
else
load(fname)


end

ps=D.p_search;
good=D.p_search;


for px=1:size(good,2)
    good(:,px)=D.p_search(:,px)==x_bf(px);
end

do_param=find(var(ps)>eps);
p_names={'lambda', 'sigma2_a', 'sigma2_s','sigma2_i','Bounds','phi','tau_phi','bias', 'lapse' ,'biased_sigma2_s','biased_input','biased_inatt'};
p_labs={'lambda', 'sigma2_a', 'sigma2_s','sigma2_i','Bounds','phi','tau_phi','bias', 'lapse' ,'Unbalanced Input Noise','Unbalanced Input Gain','Post-Categorization Bias'};

p_pairs=nchoosek(do_param,2);

for px = 1:size(p_pairs,1)
    figure(100+px); clf
    p1=p_pairs(px,1);
    p2=p_pairs(px,2);
    
    both_good=ones(size(good,1),1)==1;
    
    for ppx=1:size(good,2)
        if ppx~=p1 && ppx~=p2
            both_good=both_good & ps(:,ppx)==x_bf(ppx);
        end
    end
    
    
    f1=ps(both_good,p1);
    f2=ps(both_good,p2);
    Y=D.LL(both_good);
    zY = exp(Y-max(Y(:)));
    
    if norm_plot
        zY = sqrt(exp(Y-max(Y(:))));
    else
        zY=Y;
    end
    
    switch plot_type
        case 'scatter'
            h=scatter(f1,f2,50*ones(size(f1)),zY,'Marker','o');
            set(gca,'Color',[0.5 0.5 0.5])
            
        case 'heat'
            [x,y,z]=binned2(f1,f2,Y,'n_x_bins',n_x_bins,'n_y_bins',n_y_bins);
            z=fill_in_nan(z);
            if norm_plot
                z=(exp(z-max(z(:))));
            end
            ph=pcolor(x,y,z);
            set(ph,'LineStyle','none');
            if norm_plot
            set(gca,'Clim',[-0.03 1]);
            else
            % set(gca,'Clim',[max(z(:))-20 max(z(:))]);
            end
                
    end
    
    
    if do_fit
        nY=Y-min(Y); nY=nY/max(nY);
        
        bady=nY<thresh;
        f1=f1(~bady);
        f2=f2(~bady);
        Y=Y(~bady);
        
        
        [H, d, c, e,fitf] = quadfitN([f1 f2]', Y', 0);
        
        rx=rand(10000,1)*range(f1)+min(f1);
        ry=rand(10000,1)*range(f2)+min(f2);
        fit=fitf([rx ry]');
        
        
        
        [bx,by,bz]=binned2(rx,ry,fit','n_x_bins',30,'n_y_bins',30);
        
        hold on;
        if norm_plot
            bz=exp(bz-max(bz(:)));
        end
        mesh(bx,by,bz);
        [~,fstr]=fileparts(fname);
        %saveas(gcf,[fstr '_' p_names{p1} '_v_' p_names{p2} '.fig'])
        E(px).f1=f1;
        E(px).f2=f2;
        E(px).LL=Y;
        E(px).fit = {H, d, c, e,fit};
    end
    
    
    
    %     if norm_plot
    %         z=exp(z-max(z(:)));
    %     end
    %
    %    surf(x,y,z);
    % imagesc(x,y,z)
    if norm_plot
        colormap(hot);
    else
        colormap(jet);
    end
    cbh=colorbar;
    if norm_plot
        ylabel(cbh,'Normalized Likelihood');
    else
        ylabel(cbh,'Log Likelihood');
    end
    %  axis xy
    hx=xlabel(p_labs{p1});
    hy=ylabel(p_labs{p2});
    set([hx hy],'Interpreter','none')
    E(px).p1=p1; %#ok<*AGROW>
    E(px).p2=p2;
    E(px).allf1=ps(both_good,p1);
    E(px).allf2=ps(both_good,p2);
    E(px).allY=D.LL(both_good);
    E(px).thresh = thresh;
    
    E(px).param1 = p_names{p1};
    E(px).param2 = p_names{p2};
    set(100+px,'PaperSize',[5 4])
    set(100+px,'PaperPosition',[.5 .5 4 3])
    if norm_plot
        saveas(100+px,sprintf('~/Desktop/grid/%s_%s_norm.pdf',p_names{p1},p_names{p2}));
    else
        saveas(100+px,sprintf('~/Desktop/grid/%s_%s.pdf',p_names{p1},p_names{p2}));
    end
    close(100+px);
end


if do_fit
    for nx = 1:numel(E)
      H = E(nx).fit{1};
      SDM = sqrt(abs(inv(H)));
      sd1 = SDM(1,1);
      sd2 = SDM(2,2);
      fprintf('%s s.d.=%f, %s s.d.=%f\n',E(nx).param1,sd1,E(nx).param2,sd2)
    end
end
    
    
    

end