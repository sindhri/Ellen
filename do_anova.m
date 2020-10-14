%todo
%2. make the anova more friendly

%20201013, analysis on HOM alone and ST alone, compare pre and post
%20200929, saved tables
%try use the difference wave of pre and post, then 3-way anova
%only use the mean across the fish in that geno across time

function combined_table = do_anova(output)

    pre_table = output.pre_geno_mean;
    post_table = output.post_geno_mean;
    geno_table_index = 1;
    activity_col_index = 3:size(pre_table,2); %skip the first two columns, geno, GroupCount
    
    % combine the pre and post table, add time(either pre_ or post), 
    % and extract the geno column into three factors
    combined_table = helper_combine(pre_table, post_table,...
    geno_table_index, activity_col_index);

    % extract HOM and WT from the combined table
    pre_post_table = combined_table(strcmp(combined_table.time, 'diff')~=1,:);
    HOM_table = pre_post_table(contains(pre_post_table.factor1, 'HOM'),:);
    WT_table = pre_post_table(contains(pre_post_table.factor1, 'WT'),:);

    diff_table = combined_table(strcmp(combined_table.time, 'diff')==1,:);

    % do anova and make anova plots for each table
    destination_folder = [output.pathname 'anova/'];
    if exist(destination_folder,'dir')~=7
       mkdir(destination_folder);
    end     
    anova_plot(diff_table, 'diff', activity_col_index, destination_folder);
    anova_plot(HOM_table, 'HOM', activity_col_index, destination_folder);
    anova_plot(WT_table, 'WT', activity_col_index, destination_folder);
end

% run three types, diff, HOM, WT
function anova_plot(data_table, type, activity_col_index, destination_folder)
    if strcmp(type, 'diff') == 0 && strcmp(type, 'HOM') == 0 && strcmp(type, 'WT') == 0
        fprintf('type has to be either diff, HOM, or WT, abort\n')
        return
    end
    
    for i = 1:length(activity_col_index)
            activity_name = data_table.Properties.VariableNames{activity_col_index(i)};
            activity_data = data_table{:,activity_col_index(i)};
        % 3-way anova with interaction
        % for diff_table
            %factor1: geno-- HOM vs. WT
            %factor2: drug1 MC vs. DMSO
            %factor3: drug2 H20 vs. PTZ
        % for HOM/WT table
            %factor1: time-- pre_ vs. post
            %factor2: drug1 MC vs. DMSO
            %factor3: drug2 H20 vs. PTZ
        if strcmp(type, 'diff')==1
            [p,tbl,stats] = anovan(activity_data,{data_table.factor1,...
                data_table.factor2, data_table.factor3},...
                'model','interaction','varnames',{'Geno','drug1','drug2'});
        else
            [p,tbl,stats] = anovan(activity_data,{data_table.time,...
                data_table.factor2,data_table.factor3},...
                'model','interaction','varnames',{'Time','drug1','drug2'});
        end
        figure1 = figure(1);
        set(figure1, 'PaperPosition', [0 0 8 4],'PaperPositionMode','auto');
        set(figure1,'PaperSize',[8,4]);
        saveas(figure1,[destination_folder '/anova_' type '_' activity_name],'pdf');
        close(figure1);
        
    end
end
