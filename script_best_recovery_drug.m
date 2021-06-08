%20200524, added more inputs to run_CPA_euclidean
%20210325, added inputs and outputs to change the functions together
% no longer need to select every input file, only the original data needs
% to be selected and the rest will be generated automatically
addpath('helpers/');
[pathname, zscore_file] = calculate_zscore;
mean_by_geno_file = average_after_zscore(pathname, zscore_file);
avgz_to_clustergram(pathname, mean_by_geno_file);
d1 = run_PCA_euclidean('HOM', pathname, mean_by_geno_file, 0);
d2 = run_PCA_euclidean('HOM', pathname, mean_by_geno_file, 1);
d3 = run_PCA_euclidean('HET', pathname, mean_by_geno_file, 0);
d4 = run_PCA_euclidean('HET', pathname, mean_by_geno_file, 1);

% no need to run this one, this will generate distance between all geno
% pairs, no pre-selection
%distance_file = get_euclidean_distance(pathname, mean_by_geno_file);

% plot distance file
plot_distance_WT(pathname, distance_file);
plot_distance('HOM', pathname, filename)
