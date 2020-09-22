function plot_prepost(output)
    color_lab = {[1,0,0],[0,0,0.5],[1,0.4,0],[1,1,0],...
        [0.4,0.4,1],[0.4,1,0.4],[1,0,1],[0,1,1]};
    location_lab = repmat({'northwest'},[length(output.parameters),1]);
    % in a couple special cases put the legend in the bottom left corner
    location_lab{2} = 'southwest';
    location_lab{8} = 'southwest';
    
    pre_geno_mean = output.pre_geno_mean;
    pre_geno_std = output.pre_geno_std;
    post_geno_mean = output.post_geno_mean;
    post_geno_std = output.post_geno_std;

    geno = pre_geno_mean.geno;

    x_label = {'pre','post'};
    
    for i = 1:length(output.parameters)
        activity_name = output.parameters{i};
        %locate the correct column for this parameter
        pre_mean = pre_geno_mean{:,2+i};
        pre_std = pre_geno_std{:,2+i};
        post_mean = post_geno_mean{:,2+i};
        post_std = post_geno_std{:,2+i};
        figure;
        hold on;
        for j = 1:length(geno)
            the_color = color_lab{j};
            x = [1,4];
            y = [pre_mean(j), post_mean(j)];
            err = [pre_std(j), post_std(j)];
            e = errorbar(x,y,err);
            e.Color = the_color;
            e.LineWidth = 2;
            set(gca, 'XTick',[1,4], 'XTickLabel',x_label, 'TickLabelInterpreter','none')
            set(gca, 'XLim',[0 5])
        end
        legend(geno, 'interpreter','none','location',location_lab{i});
        title(activity_name, 'interpreter','none');
        
        set(gcf, 'PaperPosition', [0 0 10 5]); %Position plot at left hand corner with width 5 and height 5.
        set(gcf, 'PaperSize', [10 5]); 
        saveas(gcf,['plots/prepost_' activity_name],'pdf');
        close;
    end
end