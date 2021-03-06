%20210314, modified calculate_zscore. It will select a subgroup of the
%whole file, with the genotype starging with the input geno_name, then use
%geno_name + DMSO as the control

%202103, change Biochanin A to Biochanin-A, related WT to WT
%20201103, include all the controls in the zscore outputs
%20200619
%load the data
%remove NaN
%take the genotypes (HET + xxx) and (HOM + xxx)
%normalized by the mean and std of either (WT + DMSO) or (related WT + DMSO) by
%date

% 202103
% zscore has less rows than original data
% because original data has HOm_del44 and WT/HET types 
% with experiment date 200912, and they were not recognized in the zscore
% step
function [pathname, filename, filename_output] = calculate_zscore_subgroup(geno_name, pathname, filename)
if nargin==1
    [filename,pathname] = uigetfile('*.csv','select the raw data file');
end
filename_noext = strsplit(filename,'.');
filename_noext = filename_noext{1};
filename_output = [filename_noext '_zscore_' geno_name '.csv'];

fprintf('Reading file: %s%s\n', pathname, filename);

main_table = readtable([pathname filename]);

main_table.genotype = regexprep(main_table.genotype, 'Biochanin A', 'Biochanin-A');
main_table.genotype = regexprep(main_table.genotype, 'related WT', 'WT');
%decide which type of file it is, 
%file type = 1: "individualData.csv", using "plate"
%file type = 2: "scn1lab_rw_split_means.csv" using "experiment"
exist_column = strcmp('plate',main_table.Properties.VariableNames);
if(find(exist_column==1))
    file_type = 1;
else
    file_type = 2;
    main_table(:,1)=[];
    main_table.Properties.VariableNames(1) = {'plate'};
end
parameter_start_column = 4;

% Only take the rows that starts with geno_name
fprintf('picking the rows whose genotype starts with %s\n', geno_name);
main_table = main_table(contains(main_table.genotype, geno_name),:);

%replace Inf with Nan
fprintf('replacing Inf with Nan, will take a minute .....\n');
for i = 1:size(main_table, 1)
    for j = 3:size(main_table, 2)
        if isinf(main_table{i,j})
            main_table{i,j} = NaN;
        end
    end
end

%remove whole columns with missing values
fprintf('remove whole empty columns.....\n');
missing_matrix = ismissing(main_table);
to_remove = [];
for i = 1:size(missing_matrix,2)
    if isempty(find(missing_matrix(:,i)==0))
        to_remove = [to_remove,i];
    end
end
column_index_to_keep = setdiff(1:size(main_table,2),to_remove);
main_table = main_table(:,column_index_to_keep);


%create 3 new columns:
%type
%date
%genotype_plate_combined
%for information aggregation
for i = 1:size(main_table,1)
    genotype_splits = strsplit(main_table.genotype{i},' ');
    main_table.type{i} = genotype_splits{1};
    plate_splits = strsplit(main_table.plate{i},'_');
    main_table.date{i} = plate_splits{1};
end
main_table.genotype_plate_combined = strcat(main_table.genotype,...
    '__',main_table.date);

%get all the data columns
parameters = main_table.Properties.VariableNames(parameter_start_column:end-3);

%find the WT and 
geno_control = [geno_name ' + DMSO'];

fprintf('Calculate the mean and std for the control types.\n');

table_WT = main_table(strcmp(main_table.genotype,geno_control),:);

mean_function = @(x) mean(x,'omitnan');
std_function = @(x) std(x,'omitnan');

WTmeanByGenotypePlate = varfun(mean_function,table_WT,...
    'GroupingVariables',{'genotype','date'},...
    'InputVariables',{parameters{:}});

WTstdByGenotypePlate = varfun(std_function,table_WT,...
    'GroupingVariables',{'genotype','date'},...
    'InputVariables',{parameters{:}});

%add the first part of genotype as its own column,'WT','HOM','HET','related' 

fprintf('Calculate the z-scores based on the control type mean and std by date.\n');

table_analysis = main_table;

categories = unique(table_analysis.genotype_plate_combined);
table_transferred = [];
for i = 1:length(categories)

    categories_split = strsplit(categories{i},'__');
    genotype = categories_split{1};
    date = categories_split{2};
    %skip 2020 data
    if strcmp(date(1:2),'20')
        continue;
    end
    data = table_analysis(strcmp(table_analysis.genotype,...
    genotype) & strcmp(table_analysis.date,date),:);
    WT_mean = WTmeanByGenotypePlate(strcmp(WTmeanByGenotypePlate.date,...
        date),:);
    WT_std = WTstdByGenotypePlate(strcmp(WTstdByGenotypePlate.date,...
        date),:);
    data_zscore = rdivide(data{:,...
        parameter_start_column:end-3}-WT_mean{1,...
        parameter_start_column:end},WT_std{1,parameter_start_column:end});
    data_transferred = data;
    data_transferred{:,parameter_start_column:end-3} = data_zscore;
    table_transferred = [table_transferred;data_transferred];
end
writetable(table_transferred,[pathname filename_output]);
fprintf('zscore file generated: %s%s\n', pathname, filename_output);
end