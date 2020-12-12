function plot_distance_WT

[filename,pathname] = uigetfile('*.csv','select the distance file');
main_table = readtable([pathname filename]);

fig_dir = [pathname '/fig/'];
if exist(fig_dir, 'dir') ~=7
    mkdir(fig_dir);
end

control_genos_original = {'WT + DMSO', 'related WT + DMSO'};
% clean up to match the format of matlab columns
control_genos_cleaned = regexprep(control_genos_original, '_', '');
control_genos_cleaned = regexprep(control_genos_cleaned, ' ', '');
control_genos_cleaned = regexprep(control_genos_cleaned, '+', '_');

% total distance to WT, HOM only
for i = 1:length(control_genos_cleaned)
    control_geno = control_genos_cleaned{i};
   
    distance = main_table.(control_geno);
    label = main_table.euclidean_distance;
       
    table_temp = table(label, distance);
    table_sorted = sortrows(table_temp,'distance');
    table_sorted(1,:) = [];
    
    % plot 1, target
    % take only the target rows
    target_geno = 'HOM';
    table1 = get_target_rows(table_sorted, target_geno);
    writetable(table1,[pathname target_geno '_' control_geno '.csv']);
    plot_distance_bar(table1.distance, table1.label,...
        fig_dir, control_geno, target_geno);
    
    
    % plot 2, by geno initials
    temp = control_genos_original{i};
    temp2 = strsplit(temp, '+');
    target_geno = temp2{1}; % WT or related WT
    table2 = get_target_rows(table_sorted, target_geno);
    table2 = sortrows(table2, 'distance', 'descend');
    writetable(table2,[pathname target_geno '_' control_geno '.csv']);

    plot_distance_bar(table2.distance, table2.label,...
        fig_dir, control_geno, target_geno);
end
end

function plot_distance_bar(distance, labels, ...
    fig_dir, control_geno, target_geno)

%limit = 100;
title_name = [target_geno '_to_' control_geno];
%if length(distance) > limit
%    title_name = ['top ' int2str(limit) ' ' title_name];
%else
    limit = length(distance);
%end

barh(distance(1:limit));
yticklabels(labels(1:limit));
yticks(1:limit);
xticks(1:max(distance(1:limit)));

title(title_name, 'interpreter', 'none', 'fontsize', 15);

set(gcf, 'PaperPosition', [0 0 6 10]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [6 10]); 
set(gca,'Ydir','reverse', 'fontsize', 12, ...
    'xaxisLocation','top', 'TickLabelInterpreter', 'none');
saveas(gcf,[fig_dir target_geno '_to_' control_geno],'pdf');
close(gcf);

end

% take only the rows of column :label that start with target
% the table content is a cell, so had to take the nested inforation out
function T_output = get_target_rows(T, target)
m = 1;
selected_label = cell(1);
selected_distance = [];
for i = 1:size(T,1)
    the_label = T.label{i};
    if length(the_label) < length(target)
        continue;
    end
    
    sub_label = extractBetween(the_label, 1, length(target));
    if strcmp(sub_label, target)==1
        selected_label{m,1} = the_label;
        selected_distance(m,1) = T.(T.Properties.VariableNames{2})(i);
        m = m + 1;
    end    
end
T_output = table(selected_label, selected_distance);
T_output.Properties.VariableNames = T.Properties.VariableNames;
end