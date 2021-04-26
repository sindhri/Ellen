%find largest differential effect
%1. calculate z-score differentially against WT+DMSO and HOM+DMSO
%2. averege the z-scores across fish
%3. calcualte euclidean distance
%4. pair HOM and WT for the same drug
addpath('helpers/');
[pathname, filename, output1] = calculate_zscore_subgroup('WT');
[~, ~, output2] = calculate_zscore_subgroup('HOM', pathname, filename);

filename_noext = strsplit(filename,'.');
filename_noext = filename_noext{1};
output_z_differential = [filename_noext '_zscore_differential.csv'];

csv1 = readtable([pathname output1]);
csv2 = readtable([pathname output2]);
allCsv = [csv1;csv2]; % Concatenate vertically
writetable(allCsv, [pathname output_z_differential]);

mean_by_geno = average_after_zscore(pathname, output_z_differential);

distance_file = get_euclidean_distance(pathname, mean_by_geno);
differential_output = get_differential_distance(pathname, distance_file);
plot_differential_distance_bar(pathname, differential_output);