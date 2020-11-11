%select the mean_by_geno file
filename = 'individualData_07302020_zscore_mean_by_geno.csv';
main_table = readtable(filename);
data = main_table{:,2:end};

%the current data is geno_by_activity,m x n, 52 x 24
%output
%COEFF n x n
%rows: coefficients for the n activities
%columns: n pricipal components

%SCORE m x n
%for each drug in m, 
%each row is the loading for that drug on component n
%for example, for component 1, row 43 has the largest loadings

%latent: variabilities for each components
%tsquared: Hotelling?s T-squared statistic values
%explained: percentage explained by components, the first componnet took
%94% of the variance!

[coeff,score,latent,tsquared,explained,mu] = pca(data);
plot(SCORE(1:5,:)');