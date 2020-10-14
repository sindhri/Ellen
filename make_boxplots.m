%20201013, boxplots for all the effects for a 3-way ANVOVA
% for diff to compare geno x drug1 x drug2
% and for HOM, and WT to compare time x drug1 x drug2

function make_boxplots(output)

    % combine pre and post, add the time column (pre_ or post)
    pre_table = output.pre_fish_mean;
    post_table = output.post_fish_mean;
    n_col_to_skip = 3; %fish, geno, GroupCount
    geno_table_index = 2;
    activity_col_index = n_col_to_skip+1:size(pre_table,2); %skip n_col_to_skip

    combined_table = helper_combine(pre_table, post_table,...
    geno_table_index, activity_col_index);

    combined_table = permutate_all_factors(combined_table);
    
    pre_table = combined_table(strcmp(combined_table.time, 'pre')==1,:);
    post_table = combined_table(strcmp(combined_table.time, 'post')==1,:);
    diff_table = combined_table(strcmp(combined_table.time, 'diff')==1,:);

    pre_post_table = combined_table(strcmp(combined_table.time, 'diff')~=1,:);
    HOM_table = pre_post_table(contains(pre_post_table.factor1, 'HOM')==1,:);
    WT_table = pre_post_table(contains(pre_post_table.factor1, 'WT')==1,:);
    %boxplot(pre_table{:,4},pre_table.factor1);
    destination_folder = [output.pathname 'boxplot/'];
        %save the figures
    if exist(destination_folder,'dir')~=7
       mkdir(destination_folder);
    end
    
    % plot the pre_table
    data_table = pre_table;
    data_table_name = 'pre';
    for i = 1:length(activity_col_index)
        plot_data_index = activity_col_index(i);
        make_one_figure(data_table, plot_data_index,destination_folder, data_table_name);
    end
    
    % plot the post_table
    data_table = post_table;
    data_table_name = 'post';
    for i = 1:length(activity_col_index)
        plot_data_index = activity_col_index(i);
        make_one_figure(data_table, plot_data_index,destination_folder, data_table_name);
    end
    
    % plot the diff_table
    data_table = diff_table;
    data_table_name = 'diff';
    for i = 1:length(activity_col_index)
        plot_data_index = activity_col_index(i);
        make_one_figure(data_table, plot_data_index,destination_folder, data_table_name);
    end
    
    % plot the HOM_table
    data_table = HOM_table;
    data_table_name = 'HOM';
    for i = 1:length(activity_col_index)
        plot_data_index = activity_col_index(i);
        make_one_figure(data_table, plot_data_index,destination_folder, data_table_name);
    end
    
    data_table = WT_table;
    data_table_name = 'WT';
    for i = 1:length(activity_col_index)
        plot_data_index = activity_col_index(i);
        make_one_figure(data_table, plot_data_index,destination_folder, data_table_name);
    end
end

function combined_col = combine_two_columns(col1, col2)
    combined_col = cell(1);
    for i = 1:length(col1)
        combined_col{i,1} = [col1{i}, '_', col2{i}];
    end
end

function combined_col = combine_three_columns(col1, col2, col3)
    combined_col = cell(1);
    for i = 1:length(col1)
        combined_col{i,1} = [col1{i}, '_', col2{i}, '_', col3{i}];
    end
end

function data_table = permutate_all_factors(data_table)
    data_table.factor12 = combine_two_columns(data_table.factor1,...
        data_table.factor2);
   data_table.factor13 = combine_two_columns(data_table.factor1,...
        data_table.factor3);
    data_table.factor23 = combine_two_columns(data_table.factor2,...
        data_table.factor3);
    data_table.factor123 = combine_three_columns(data_table.factor1,...
        data_table.factor2, data_table.factor3);

    data_table.factort2 = combine_two_columns(data_table.time,...
        data_table.factor2);
    data_table.factort3 = combine_two_columns(data_table.time,...
        data_table.factor3);
    data_table.factort23 = combine_three_columns(data_table.time,...
        data_table.factor2, data_table.factor3);
end

function make_one_figure(data_table, plot_data_index,...
    destination_folder, data_table_name)
    figure;
    set(gcf, 'PaperPosition', [0 0 15 8]); %Position plot at left hand corner with width 5 and height 5.
    set(gcf, 'PaperSize', [15 8]);
    titlename = [data_table_name '_' data_table.Properties.VariableNames{plot_data_index}];
    sgtitle(titlename, 'interpreter','none');
    
    % set up the grid
%    pos1 = get(subplot(3,3,1),'Position');
    pos1 = [0.1300    0.7093    0.2134    0.2157];
%    pos2 = get(subplot(3,3,2),'Position');
    pos2 = [0.4108    0.7093    0.2134    0.2157];
%    pos3 = get(subplot(3,3,3),'Position');
    pos3 = [0.6916    0.7093    0.2134    0.2157];
    %pos4 = get(subplot(3,3,4),'Position');
    pos4 = [0.1300    0.4096    0.2134    0.2157];
    %pos5 = get(subplot(3,3,5),'Position');
    pos5 = [0.4108    0.4096    0.2134    0.2157];
    %pos6 = get(subplot(3,3,6),'Position');
    pos6 = [0.6916    0.4096    0.2134    0.2157];
    %pos7 = get(subplot(3,3,[7,8]),'Position');
    pos7 = [0.1300    0.1100    0.4942    0.2157];
    pos1 = adjust_col1(pos1);
    pos2 = adjust_col2(pos2);
    pos3 = adjust_col3(pos3);
    pos4 = adjust_col1(pos4);
    pos5 = adjust_col2(pos5);
    pos6 = adjust_col3(pos6);
    pos7 = adjust_col1(pos7);
    
    
   %plot level 3
    hs = axes;
    set(hs,'OuterPosition',pos7);
    if strcmp(data_table_name,'HOM')==1 || strcmp(data_table_name,'WT')==1
        yl = make_one_boxplot(data_table, plot_data_index, data_table.factort23);
    else
        yl = make_one_boxplot(data_table, plot_data_index, data_table.factor123);
    end
    
    %plot level 1
    hs = axes;
    set(hs,'OuterPosition',pos1);
    if strcmp(data_table_name,'HOM')==1 || strcmp(data_table_name,'WT')==1
       make_one_boxplot(data_table, plot_data_index, data_table.time);
    else
       make_one_boxplot(data_table, plot_data_index, data_table.factor1);
    end
    ylim(yl);
    
    hs = axes;
    set(hs,'OuterPosition',pos2);
    make_one_boxplot(data_table, plot_data_index, data_table.factor2);

    ylim(yl);
    
    hs = axes;
    set(hs,'OuterPosition',pos3);
    pause(1);
    make_one_boxplot(data_table, plot_data_index, data_table.factor3);
    ylim(yl);

    %plot level 2
    hs = axes;
    set(hs,'OuterPosition',pos4);
    if strcmp(data_table_name,'HOM')==1 || strcmp(data_table_name,'WT')==1
        make_one_boxplot(data_table, plot_data_index, data_table.factort2);
    else
        make_one_boxplot(data_table, plot_data_index, data_table.factor12);
    end
    ylim(yl);

    hs = axes;
    set(hs,'OuterPosition',pos5);
    if strcmp(data_table_name,'HOM')==1 || strcmp(data_table_name,'WT')==1
        make_one_boxplot(data_table, plot_data_index, data_table.factort3);
    else
        make_one_boxplot(data_table, plot_data_index, data_table.factor13);
    end
    ylim(yl);

    hs = axes;
    set(hs,'OuterPosition',pos6);
    make_one_boxplot(data_table, plot_data_index, data_table.factor23);
    ylim(yl);

    saveas(gcf,[destination_folder titlename],'pdf');
    saveas(gcf,[destination_folder titlename],'fig');
    close;

end
function pos = adjust_col1(pos)
    pos(1) = pos(1) - 0.15;
    pos(3) = pos(3) * 1.75;
end

function pos = adjust_col2(pos)
    pos(1) = pos(1) - 0.1;
    pos(3) = pos(3) * 1.75;
end

function pos = adjust_col3(pos)
    pos(1) = pos(1) - 0.05;
    pos(3) = pos(3) * 1.75;
end

function final_label = create_friendly_label(xticklabel)
    n_label = length(xticklabel);
    n_factors = length(find(xticklabel{1}=='_'))+1;
        result_split = cell(n_label,n_factors);
        for i = 1:n_label
            row = xticklabel{i};
            result_split(i,:) = strsplit(row, '_');
        end
    final_label = [];
    for i = 1:size(result_split,1)
        for j = 1:size(result_split,2)
            if mod(j,n_factors)~=0
               final_label = [final_label sprintf([result_split{i,j} '\\newline'])];
            else
               final_label = [final_label sprintf([result_split{i,j} '\n'])];
            end
        end
    end
end

function yl = make_one_boxplot(data_table, plot_data_index, label_col)
    boxplot(data_table{:,plot_data_index}, label_col);
    yl = ylim;
    xticklabel = get(gca,'xticklabel');
    final_label = create_friendly_label(xticklabel);
    set(gca,'XTickLabel',final_label,'TickLabelInterpreter','tex');
end