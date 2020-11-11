%20201022, modify the model to use individual fish
%only run the 3-way interaction, but run 3 times
%1. difference score: HOM_WT * MC_DMS * H2O_PTZ
%2. HOM: pre/post * MC_DMSO * H2O_PTZ
%3. WT: pre/post * MC_DMSO * H2O_PTZ
%add a post hoc t test
%friendly output of all the results, to a csv file
 
%20201013, analysis on HOM alone and ST alone, compare pre and post
%20200929, saved tables
%try use the difference wave of pre and post, then 3-way anova
%only use the mean across the fish in that geno across time
 
function [output_table, results] = do_anova(output)
 
    pre_table = output.pre_fish_mean;
    post_table = output.post_fish_mean;
    geno_table_index = 2;
    col_to_skip = 3;
    activity_col_index = col_to_skip+1:size(pre_table,2); %skip the first few columns, geno, GroupCount
    
    % combine the pre and post table, add time(either pre_ or post), 
    % and extract the geno column into three factors
    combined_table = helper_combine(pre_table, post_table,...
    geno_table_index, activity_col_index);
 
    % extract HOM and WT from the combined table
    pre_post_table = combined_table(strcmp(combined_table.time, 'diff')~=1,:);
    HOM_table = pre_post_table(contains(pre_post_table.factor1, 'HOM'),:);
    WT_table = pre_post_table(contains(pre_post_table.factor1, 'WT'),:);
 
    diff_table = combined_table(strcmp(combined_table.time, 'diff')==1,:);
    
    results = struct;
    results = anova_plot(diff_table, 'diff', activity_col_index, results);
    results = anova_plot(HOM_table, 'HOM', activity_col_index, results);
    results = anova_plot(WT_table, 'WT', activity_col_index, results);
    
    output_table = write_results_doc(results, output.pathname);
end
 
% run three types, diff, HOM, WT
function results = anova_plot(data_table, type, activity_col_index,results)
    if strcmp(type, 'diff') == 0 && strcmp(type, 'HOM') == 0 && strcmp(type, 'WT') == 0
        fprintf('type has to be either diff, HOM, or WT, abort\n')
        return
    end
    for i = 1:length(activity_col_index)
            
            activity_name = data_table.Properties.VariableNames{activity_col_index(i)};
            activity_data = data_table{:,activity_col_index(i)};
            results.row_names{i,1} = activity_name;
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
            [~,tbl,~] = anovan(activity_data,{data_table.factor1,...
                data_table.factor2, data_table.factor3},...
                'model',[1,1,1],'varnames',{'HOM_WT','MC_DMSO','H2O_PTZ'},...
                'display','off');
            results.data{i,1} = [tbl{2,3}, tbl{3,3}, tbl{2,6},tbl{2,7}];
            
%            figure1 = figure(1);
%            set(figure1, 'PaperPosition', [0 0 8 4],'PaperPositionMode','auto');
%            set(figure1,'PaperSize',[8,4]);
%            saveas(figure1,[destination_folder '/anova_' anova_name],'pdf');
%            close(figure1);
            
            
            %only for validation, results were similar to repeated measures
 
            %HOM_table = data_table(contains(data_table.factor1,'HOM'),:);
            %WT_table = data_table(contains(data_table.factor1,'WT'),:);
            
            %HOM_data = HOM_table{:,activity_col_index(i)};
            %WT_data = WT_table{:,activity_col_index(i)};
            
            %[p,tbl,stats] = anovan(HOM_data,{HOM_table.factor2,...
            %    HOM_table.factor3},...
            %    'model',[1,1],'varnames',{'drug1','drug2'});
            %results.data{i,4} = [tbl{2,3}, tbl{3,3}, tbl{2,6},tbl{2,7}];
            %close(figure(1));
            %[p,tbl,stats] = anovan(WT_data,{WT_table.factor2,...
            %    WT_table.factor3},...
            %    'model',[1,1],'varnames',{'drug1','drug2'});
            %results.data{i,5} = [tbl{2,3}, tbl{3,3}, tbl{2,6},tbl{2,7}];
            %close(figure(1));
        else
            % build the repeated-measures data table
            activity_table = table(activity_data, data_table.time, data_table.fish);
            activity_table_unstacked = unstack(activity_table,'activity_data', 'Var2');
            activity_table_unstacked.Properties.VariableNames{'Var3'} = 'fish';
            data_table_pre = data_table(strcmp(data_table.time,'pre'),:);
            data_table_pre = table(data_table_pre.fish, ...
                data_table_pre.factor2, data_table_pre.factor3);
            data_table_pre.Properties.VariableNames = {'fish','factor2',...
                'factor3'};
            anova_table = join(activity_table_unstacked, ...
                data_table_pre, 'keys','fish');
            Meas = table([1 2]','VariableNames',{'Measurements'});
 
            % run repeated-measures ANOVA
            rm = fitrm(anova_table,'pre-post~factor2*factor3','WithinDesign',Meas); 
            tbl = ranova(rm);
            if strcmp(type,'HOM')==1
               results.data{i,2} = [tbl{4,2}, tbl{5,2}, tbl{4,4},tbl{4,5}];
            else
               results.data{i,3} = [tbl{4,2}, tbl{5,2}, tbl{4,4},tbl{4,5}];
            end
            
            %do t tests for HOM
            if strcmp(type, 'HOM')==1
                HOM_DMSO_PTZ = anova_table(strcmp(anova_table.factor2,...
                    'DMSO') & strcmp(anova_table.factor3,...
                    'PTZ'),2:3);
                HOM_10MC_PTZ = anova_table(strcmp(anova_table.factor2,...
                    '10MC') & strcmp(anova_table.factor3,...
                    'PTZ'),2:3);
                % pre: comparing DMSO vs. 10MC
                [~,p,~,stats] = ttest2(HOM_DMSO_PTZ{:,2},HOM_10MC_PTZ{:,2});
                results.data{i,4} = [stats.df, stats.tstat, p];
                % post: comparing DMSO vs. 10MC
                [~,p,~,stats] = ttest2(HOM_DMSO_PTZ{:,1},HOM_10MC_PTZ{:,1});
                results.data{i,5} = [stats.df, stats.tstat, p];
            else
 
                WT_DMSO_PTZ = anova_table(strcmp(anova_table.factor2,...
                    'DMSO') & strcmp(anova_table.factor3,...
                    'PTZ'),2:3);
                WT_10MC_PTZ = anova_table(strcmp(anova_table.factor2,...
                    '10MC') & strcmp(anova_table.factor3,...
                    'PTZ'),2:3);
                % pre: comparing DMSO vs. 10MC
                [~,p,~,stats] = ttest2(WT_DMSO_PTZ{:,2},WT_10MC_PTZ{:,2});
                results.data{i,6} = [stats.df, stats.tstat, p];
                % post: comparing DMSO vs. 10MC
                [~,p,~,stats] = ttest2(WT_DMSO_PTZ{:,1},WT_10MC_PTZ{:,1});
                results.data{i,7} = [stats.df, stats.tstat, p];   
            end
        end
        
    end
end
 
function output_table = write_results_doc(results, pathname)
%    fid = fopen('anova_results.txt','w');
%    fprintf(fid,'Seizure analysis ANOVA results\n');
    %analysis_titles = {'Difference wave ANOVA 3-way interaction: WT/HOM * DMSO/MC * H2O/PTZ',...
    %    'HOM ANOVA 3-way interaction: pre/post * DMSO/MC * H2O/PTZ',...
    %    'WT ANOVA 3-way interaction: pre/post * DMSO/MC * H2O/PTZ',...
    %    'Independent samples t tests PRE: HOM/DMSO/PTZ vs HOM/MC/PTZ',...
    %    'Independent samples t tests POST: HOM/DMSO/PTZ vs HOM/MC/PTZ'};
    analysis_titles_short = {'Diff_ANOVA', 'HOM_ANOVA', 'WT_ANOVA',...
        'ttest_PRE-HOM_DMSO-PTZ_vs_MC-PTZ','ttest_POST-HOM_DMSO-PTZ_vs_MC-PTZ'};
    output_table = table(results.row_names);
    output_table.Properties.VariableNames = {'Activity'};
    for i = 1:length(analysis_titles_short)
        title = analysis_titles_short{i};
        column = cell(1);
        for j = 1:length(results.row_names)
            cell_content = results.data{j,i};
            if i <=3 %ANOVA
                column{j,1} = anova_to_txt(cell_content);
            else
                column{j,1} = ttest_to_txt(cell_content);
            end
        end
        output_table = [output_table, table(column)];
        output_table.Properties.VariableNames{i+1} = analysis_titles_short{i};
    end
    filename = 'anova_results.csv';
    writetable(output_table,[pathname filename]);
    msgbox(sprintf(['ANOVA results saved in\n' pathname filename]));
end
 
%anova_result follows: df1, df2, F, p
function output = anova_to_txt(anova_result)
    df1 = anova_result(1);
    df2 = anova_result(2);
    F = anova_result(3);
    p = anova_result(4);
    if p < 0.001
        p_string = 'p < .001';
    else
       p_string = sprintf('p = %.2f', p);
    end
    output = sprintf('F(%d, %d) = %.2f, %s', df1, df2, F, p_string);
end
 
%anova_result follows: df1, df2, F, p
function output = ttest_to_txt(ttest_result)
    df = ttest_result(1);
    t = ttest_result(2);
    p = ttest_result(3);
    if p < 0.001
        p_string = 'p < .001';
    else
       p_string = sprintf('p = %.2f', p);
    end
    output = sprintf('t(%d) = %.2f, %s', df, t, p_string);
end
