%20210308, changed colors, split into 4 lines per plot
%20200929, added saving figures in fig
%20200922, used the actual column names from the aggregated table
function plot_geno_by_time(output)
    if exist('plots','dir')==0
        mkdir('plots');
    end
    
    %generate a plot for each parameter, skipping the first three columns
    %number of columns to skip, geno, start, GroupCount 
    ncolumn_to_skip = 3;
    activity_names = output.pre_geno_time_mean.Properties.VariableNames(ncolumn_to_skip+1:end);
    n_activities = length(activity_names);
    % in general put the legend in the top left corner
    location_lab = repmat({'northwest'},[n_activities,1]);
    % in a couple special cases put the legend in the bottom left corner
    location_lab{2} = 'southwest';
    location_lab{8} = 'southwest';
            
    if exist([output.pathname 'line/'],'dir')~=7
       mkdir([output.pathname 'line/']);
    end

    for i = 1:n_activities
        activity_name = activity_names{i};
        activity_index = i + ncolumn_to_skip; %skipping geno, start, GroupCount
        
        [DMSO_mean1, nonDMSO_mean1] = preprocess(output.pre_geno_time_mean, activity_index);
        [DMSO_std1, nonDMSO_std1] = preprocess(output.pre_geno_time_std, activity_index);
        [DMSO_mean2, nonDMSO_mean2] = preprocess(output.post_geno_time_mean, activity_index);
        [DMSO_std2, nonDMSO_std2] = preprocess(output.post_geno_time_std, activity_index);
        plot_DMSO(DMSO_mean1, DMSO_std1, DMSO_mean2, DMSO_std2,...
            'DMSO', location_lab{i}, activity_name, output.pathname);
        if size(nonDMSO_mean1,1) > 1
            plot_DMSO(nonDMSO_mean1, nonDMSO_std1, nonDMSO_mean2, nonDMSO_std2,...
            'nonDMSO', location_lab{i}, activity_name, output.pathname);
        end
    end
    
end

function plot_DMSO(DMSO_mean1, DMSO_std1, DMSO_mean2, DMSO_std2,...
    title_prefix, legend_location, activity_name, pathname)
        figure;
        hold on;
        %color_lab = {[1,0,0],[0,0,0.5],[0.88,0.68,0],[1,0,0],...
        %[0.4,0.4,1],[0.4,1,0.4],[0,1,1],[0,0,1]};
        color_lab = {[0.88,0.68,0],[1,0,0],[0,1,1],[0,0,1]};
        h = [];
        legend_text = cell(1);
        sequence = [4,2,3,1]; 
        %to follow a set sequence: 
        %WT + DMSO + PTZ
        %HOM + DMSO + PTZ
        %WT + DMSO + water
        %HOM + DMSO + water
        for j = 1:length(sequence)
            current_row = sequence(j);
            current_geno = DMSO_mean1.geno{current_row};
            h(j) = plot_m_std(DMSO_mean1{current_row,2:end}, DMSO_std1{current_row,2:end},...
                DMSO_mean2{current_row,2:end}, DMSO_std2{current_row,2:end},color_lab{current_row});
            legend_text{j} = current_geno;
        end
        legend(h, legend_text, 'interpreter','none','location',legend_location);
        title([title_prefix '_' activity_name], 'interpreter','none');


        set(gcf, 'PaperPosition', [0 0 6 4]); %Position plot at left hand corner with width 5 and height 5.
        set(gcf, 'PaperSize', [6 4]); 
        saveas(gcf,[pathname 'line/geno_time_' title_prefix '_' activity_name],'pdf');
        saveas(gcf,[pathname 'line/geno_time_' title_prefix '_' activity_name],'fig');
        close;
end

function [DMSO, non_DMSO] = preprocess(one_table, column_to_plot)
    %column1 geno, column2 start
    table_selected = one_table(:,[1,2,column_to_plot]);
    table_processed = unstack(table_selected, ...
        one_table.Properties.VariableNames(column_to_plot),'start');
    DMSO = table_processed(contains(table_processed.geno, 'DMSO'),:);
    non_DMSO = table_processed(~contains(table_processed.geno, 'DMSO'),:);
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
    set(h,'facealpha',0.8);
    plot(x1_pre, m_pre, 'Color', the_color, 'LineWidth', 2);
    
    x1_post = length(m_pre)+1:length(m_pre)*2;
    curve1_post = m_post + std_post;
    curve2_post = m_post - std_post;
    x2_post = [x1_post, fliplr(x1_post)];
    inBetween = [curve1_post, fliplr(curve2_post)];
    h = fill(x2_post, inBetween, the_color);
    set(h,'facealpha',.8);
    plot(x1_post, m_post, 'Color', the_color, 'LineWidth', 2);
end