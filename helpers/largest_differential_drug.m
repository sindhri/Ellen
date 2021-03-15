function largest_differential_drug
[pathname, filename, output1] = calculate_zscore_subgroup('WT');
[~, ~, output2] = calculate_zscore_subgroup('HOM', pathname, filename);

filename_noext = strsplit(filename,'.');
filename_noext = filename_noext{1};
output_z_differential = [filename_noext '_zscore_differential.csv'];

csv1 = readtable([pathname output1]);
csv2 = readtable([pathname output2]);
allCsv = [csv1;csv2]; % Concatenate vertically
writetable(allCsv, [pathname output_z_differential]);

[mean_by_geno, output_rms] = average_after_zscore(pathname, output_z_differential);

distance_file = get_euclidean_distance(pathname, mean_by_geno);
get_differential_distance(pathname, distance_file);