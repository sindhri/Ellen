%20210617, implement variations for switching out target and whether to
%exclude sleep latency and sleep length. put results in subfolders
%find largest differential effect
%1. calculate z-score differentially against WT+DMSO and HOM+DMSO
%2. averege the z-scores across fish
%3. calcualte euclidean distance
%4. pair HOM/HET and WT for the same drug
%5. calculate the averaged differential effect for each drug

function largest_differential_drug(target, to_exclude_LL)
[pathname, filename, output1] = calculate_zscore_subgroup('WT');
[~, ~, output2] = calculate_zscore_subgroup(target, pathname, filename);

if to_exclude_LL == 0
    pathname_subfolder = [pathname target '_all_params' filesep];
else
    pathname_subfolder = [pathname target '_sleep_LL_excluded' filesep];
end
if exist(pathname_subfolder,'dir')~=7
    mkdir(pathname_subfolder)
end

filename_noext = strsplit(filename,'.');
filename_noext = filename_noext{1};
output_z_differential = [filename_noext '_zscore_differential.csv'];

csv1 = readtable([pathname output1]);
csv2 = readtable([pathname output2]);
allCsv = [csv1;csv2]; % Concatenate vertically

writetable(allCsv, [pathname_subfolder output_z_differential]);

mean_by_geno = average_after_zscore(pathname_subfolder, output_z_differential);
distance_file = get_euclidean_distance(pathname_subfolder, mean_by_geno, to_exclude_LL);
differential_distance_file = get_differential_distance(target, pathname_subfolder, distance_file);
plot_differential_distance_bar(pathname_subfolder, differential_distance_file);

differential_distance_avg_file = avg_over_doses(pathname_subfolder, differential_distance_file);
plot_differential_distance_bar(pathname_subfolder, differential_distance_avg_file);

end