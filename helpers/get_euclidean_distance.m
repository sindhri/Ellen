function filename_output = get_euclidean_distance(pathname, filename)

if nargin==0
    [filename,pathname] = uigetfile('*.csv','select the mean_by_geno file');
end
main_table = readtable([pathname filename]);
filename_noext = strsplit(filename,'.');
filename_noext = filename_noext{1};
filename_output = [filename_noext '_distance.csv'];

%skip the first two columns, genotype and GroupCount
data = main_table{:,3:end}; 

%data, 52 x 24, geno_by_activity
%original program 1st dimension is 120, so it is likely to be geno instead
%of activity
%which way makes more sense? to normalize across geno, or activity?
%If normalize across geno, it would compare the relative strength of geno
%for each individual activity, 
%if normalize across activity, it would compare the relative strength of
%activity, for each individual geno
%I think normalize across geno makes more sense.
%However in papers the plotting has the dimeison of geno x activity
%On the other hand, both papers seems to trying to
%do PCA on the dimension of activity, which would have activity as column
%and study EACH geno, how much it is relfected on a bundle of activites
%I guess the assumption is that activities could be highly correlated
%so trying to find the bundle of activities that best representing the
%whole picture is beneficial
%Then study EACH geno, (instead of bundle geno into groups)


%normalize across activity (2nd dimension), for each geno (1st dimension)
[n_geno, ~] = size(data);
for i = 1:n_geno
  logdata(i,:) = data(i,:)-mean(data(i,:)); % make sure there is no DC offset)
  normdata(i,:) = logdata(i,:)./std(logdata(i,:)); 
end

fig_dir = [pathname '/fig/'];
if exist(fig_dir, 'dir') ~=7
    mkdir(fig_dir);
end

% decompose the data matrix (i.e. do principal components analysis
[u, fulls, fullv] = svd(normdata);

% and get the projections-- but I will try the first 10
fullproj = normdata*fullv(:,1:10);

euclidean_distances = pdist(fullproj);
formated_distances = squareform(euclidean_distances);
distance_table = array2table(formated_distances);
distance_table.Properties.VariableNames = main_table.genotype;

distance_table = [main_table.genotype, distance_table];
distance_table.Properties.VariableNames{1} = 'euclidean_distance';
writetable(distance_table, [pathname filename_output]);
fprintf('Euclidean distance file generated: %s%s\n', pathname, filename_output);
end