function filename_output = get_differential_distance(pathname, distance_file)
if nargin==0
    [distance_file, pathname] = uigetfile('*.csv','select the distance file');    
end

filename_noext = strsplit(distance_file,'.');
filename_noext = filename_noext{1};
filename_output = [filename_noext '_output.csv'];


main_table = readtable([pathname distance_file]);

for i = 1:size(main_table,1)
    genotype_splits = strsplit(main_table.drug{i},' + ');
    main_table.geno{i} = genotype_splits{1};
    main_table.drug{i} = genotype_splits{2};
end

total_drug = length(unique(main_table.drug));
drug = cell(1);
differential_distance = zeros(total_drug, 1);
for i = 1:total_drug
    geno_name = main_table.drug{i};
    row_name = ['HOM + ' geno_name];
    col_name = ['WT + ' geno_name];
    col_name = regexprep(col_name, ' ', '');
    col_name = regexprep(col_name, '-', '_');
    col_name = regexprep(col_name, '\.', '_');
    col_name = regexprep(col_name, '+', '_');

    one_distance = main_table{strcmp(main_table.euclidean_distance,row_name)==1, col_name};
    drug{i,1} = geno_name;
    differential_distance(i,1) = one_distance;
end

differential_table = table(drug, differential_distance);
differential_table = sortrows(differential_table, 'differential_distance',...
    'descend');
writetable(differential_table,[pathname filename_output]);
