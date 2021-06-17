% 20210617
% Take the differential distance file, and average all the doses for each
% given drug
function filename_output = avg_over_doses(pathname, differential_distance_file)

if nargin==0
    [differential_distance_file, pathname] = uigetfile('*.csv','select the differential distance file');    
end

filename_noext = strsplit(differential_distance_file,'.');
filename_noext = filename_noext{1};
filename_output = [filename_noext '_avg_over_doses.csv'];

main_table = readtable([pathname differential_distance_file]);
%remove DMSO
main_table = main_table(strcmp(main_table.drug, 'DMSO')==0,:);
main_table.dose_cell = regexp(main_table.drug, 'M', 'split','once');
for i = 1:size(main_table,1)
    temp = main_table.dose_cell{i}{1};
    main_table.dose{i} = temp(1:length(temp)-1);
    temp = main_table.dose_cell{i}{2};
    temp2 = temp(2:length(temp));
    if strcmp(temp2(1:2), '_-')==1
        main_table.drug_name{i} = temp2(3:length(temp2));
    else
        main_table.drug_name{i} = temp2;
    end
end

[group, id] = findgroups(main_table.drug_name);
func = @(d)[mean(d)];
result = splitapply(func, main_table.differential_distance, group);
drug_names = array2table(id,'VariableNames', {'drug'});
differential_distance_avg = array2table(result,'VariableNames',{'differential_distance_avg'});
output = [drug_names, differential_distance_avg];
writetable(output, [pathname filename_output]);

end