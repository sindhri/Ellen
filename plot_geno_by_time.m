function plot_geno_by_time(output)
    if exist('plots','dir')==0
        mkdir('plots');
    end
    
    %generate a plot for each parameter, skipping the first three columns
    %which is geno, start, GroupCount 
    n_params = length(output.parameters);
    for i = 1:n_params
        activity_name = output.parameters{i};
        param_index = i + 3; %skipping geno, start, GroupCount
        
        mean1 = preprocess(output.pre_geno_time_mean, param_index);
        std1 = preprocess(output.pre_geno_time_std, param_index);
        mean2 = preprocess(output.post_geno_time_mean, param_index);
        std2 = preprocess(output.post_geno_time_std, param_index);
        
        figure;
        hold on;
        color_lab = {'b','r','y','m','k',[0.4,0.15,0.3], [0.4, 0.9, 0.76],[0.3,0.55,0.87]};
        % in general put the legend in the top left corner
        location_lab = repmat({'northwest'},[n_params,1]);
        % in a couple special cases put the legend in the bottom left corner
        location_lab{2} = 'southwest';
        location_lab{8} = 'southwest';
        h = [];
        for j = 1:length(output.geno)
            h(j) = plot_m_std(mean1{j,2:end}, std1{j,2:end},...
                mean2{j,2:end}, std2{j,2:end},color_lab{j});
        end
        legend(h, mean1.geno, 'interpreter','none','location',location_lab{i});
        title(activity_name);
        
        set(gcf, 'PaperPosition', [0 0 20 10]); %Position plot at left hand corner with width 5 and height 5.
        set(gcf, 'PaperSize', [20 10]); 
        saveas(gcf,['plots/geno_time_' activity_name],'pdf');
        close;
    end
    %plotting 8 types are too much. need to split or aggregate
    %here split to 
    
end

function table_processed = preprocess(one_table, column_to_plot)
    %column1 geno, column2 start
    table_selected = one_table(:,[1,2,column_to_plot]);
    table_processed = unstack(table_selected, ...
        one_table.Properties.VariableNames(column_to_plot),'start');
end

% input, m_pre, std_pre, m_post, std_post, color
% output a handler for later generating a matching legend
function h = plot_m_std(m_pre, std_pre, m_post, std_post, the_color)

    x1_pre = 1:length(m_pre);
    curve1_pre = m_pre + std_pre;
    curve2_pre = m_pre - std_pre;
    x2_pre = [x1_pre, fliplr(x1_pre)];
    inBetween = [curve1_pre, fliplr(curve2_pre)];
    h = fill(x2_pre, inBetween, the_color);
    set(h,'facealpha',.5);
    plot(x1_pre, m_pre, 'Color', the_color, 'LineWidth', 2);
    
    x1_post = length(m_pre)+1:length(m_pre)*2;
    curve1_post = m_post + std_post;
    curve2_post = m_post - std_post;
    x2_post = [x1_post, fliplr(x1_post)];
    inBetween = [curve1_post, fliplr(curve2_post)];
    h = fill(x2_post, inBetween, the_color);
    set(h,'facealpha',.5);
    plot(x1_post, m_post, 'Color', the_color, 'LineWidth', 2);
end