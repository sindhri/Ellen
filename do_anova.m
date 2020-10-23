%20201022, modify the model to use individual fish
%only run the 3-way interaction
%add a post hoc t test
%friendly output of all the results

%20201013, analysis on HOM alone and ST alone, compare pre and post
%20200929, saved tables
%try use the difference wave of pre and post, then 3-way anova
%only use the mean across the fish in that geno across time

function results = do_anova(output)

    pre_table = output.pre_fish_mean;
    post_table = output.post_fish_mean;
    geno_table_index = 2;
    col_to_skip = 3;
    activity_col_index = col_to_skip+1:size(pre_table,2); %skip the first two columns, geno, GroupCount
    
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
    results = struct;
    results = anova_plot(diff_table, 'diff', activity_col_index, ...
        destination_folder,results);
    results = anova_plot(HOM_table, 'HOM', activity_col_index, ...
        destination_folder,results);
    results = anova_plot(WT_table, 'WT', activity_col_index, ...
        destination_folder,results);
end

% run three types, diff, HOM, WT
function results = anova_plot(data_table, type, activity_col_index, destination_folder,results)
    if strcmp(type, 'diff') == 0 && strcmp(type, 'HOM') == 0 && strcmp(type, 'WT') == 0
        fprintf('type has to be either diff, HOM, or WT, abort\n')
        return
    end
    for i = 1:length(activity_col_index)
            activity_name = data_table.Properties.VariableNames{activity_col_index(i)};
            activity_data = data_table{:,activity_col_index(i)};
            anova_name = [type '_' activity_name];
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
                'model',[1,1,1],'varnames',{'Geno','drug1','drug2'});
            results.data{i,1} = [tbl{2,3}, tbl{3,3}, tbl{2,6},tbl{2,7}];
        else
            [p,tbl,stats] = anovan(activity_data,{data_table.time,...
                data_table.factor2,data_table.factor3},...
                'model',[1,1,1],'varnames',{'Time','drug1','drug2'});
            if strcmp(type,'HOM')==1
               results.data{i,2} = [tbl{2,3}, tbl{3,3}, tbl{2,6},tbl{2,7}];
            else
               results.data{i,3} = [tbl{2,3}, tbl{3,3}, tbl{2,6},tbl{2,7}];
            end

        end
        
        results.names{i} = anova_name;
        figure1 = figure(1);
        set(figure1, 'PaperPosition', [0 0 8 4],'PaperPositionMode','auto');
        set(figure1,'PaperSize',[8,4]);
        saveas(figure1,[destination_folder '/anova_' anova_name],'pdf');
        close(figure1);
    end
end
