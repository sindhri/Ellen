%20210325, added inputs and outputs to change the functions together
% no longer need to select every input file, only the original data needs
% to be selected and the rest will be generated automatically
addpath('helpers/');
[pathname, zscore_file] = calculate_zscore;
mean_by_geno_file = average_after_zscore(pathname, zscore_file);
avgz_to_clustergram(pathname, mean_by_geno_file);
run_PCA_euclidean(pathname, mean_by_geno_file);
distance_file = get_euclidean_distance(pathname, mean_by_geno_file);
plot_distance_WT(pathname, distance_file);
