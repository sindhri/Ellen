%20200619
%load the data
%remove NaN
%take the genotypes (HET + xxx) and (HOM + xxx)
%normalized by the mean and std of either (WT + DMSO) or (related WT + DMSO) by
%date

%20200814, unified the output to be the similar structure as in the
%nonburst outpt and individualData.csv

[filename,pathname] = uigetfile('*.csv','select the zscore file');
filename_noext = strsplit(filename,'.');
filename_noext = filename_noext{1};
filename_output = [filename_noext '_zscore_burst.csv'];
main_table = readtable([pathname filename]);

%removed first column and renamed 'experiment' to be 'plate'
%so it is consistent with 'individualData.csv'
main_table(:,1)=[];
main_table.Properties.VariableNames(1) = {'plate'};
    
%remove whole columns with missing values
missing_matrix = ismissing(main_table);
to_remove = [];
for i = 1:size(missing_matrix,2)
    if isempty(find(missing_matrix(:,i)==0))
        to_remove = [to_remove,i];
    end
end
column_index_to_keep = setdiff(1:size(main_table,2),to_remove);
main_table = main_table(:,column_index_to_keep);

for i = 1:size(main_table,1)
    temp = strsplit(main_table.genotype{i},' ');
    main_table.type{i} = temp{1};
    temp = strsplit(main_table.plate{i},'_');
    main_table.date{i} = temp{1};
end
main_table.genotype_plate_combined = strcat(main_table.genotype,...
    '__',main_table.date);

%genotype = unique(main_table.genotype);
index_start = 4;
index_end = 17;
parameters = main_table.Properties.VariableNames(index_start:index_end);

%find the WT and 
geno_control1 = 'WT + DMSO';
geno_control2 = 'related WT + DMSO';
geno_exp1 = 'HET';
geno_exp2 = 'HOM';

table_WT = main_table(strcmp(main_table.genotype,...
    geno_control1) | strcmp(main_table.genotype,geno_control2),:);


mean_function = @(x) mean(x,'omitnan');
std_function = @(x) std(x,'omitnan');

%calculate mean and std, 

WTmeanByGenotypePlate = varfun(mean_function,table_WT,...
    'InputVariables',{parameters{:}});

WTstdByGenotypePlate = varfun(std_function,table_WT,...
    'InputVariables',{parameters{:}});
WTstdByGenotypePlate.Fun_burct_night2_mean = WTstdByGenotypePlate.Fun_burct_night3_mean;
WTstdByGenotypePlate.Fun_burdur_night2_mean = WTstdByGenotypePlate.Fun_burdur_night3_mean;

%add the first part of genotype as its own column,'WT','HOM','HET','related' 

table_analysis = main_table(strcmp(main_table.type,...
    geno_exp1) | strcmp(main_table.type,geno_exp2),:);

categories = unique(table_analysis.genotype);
table_transferred = [];
for i = 1:length(categories)
    genotype = categories{i};
    data = table_analysis(strcmp(table_analysis.genotype,genotype),1:index_end);
    WT_mean = WTmeanByGenotypePlate;
    WT_std = WTstdByGenotypePlate;
    data_temp = rdivide(data{:,index_start:index_end}-WT_mean{:,:},WT_std{:,:});
    data_transferred = data;
    data_transferred{:,index_start:index_end} = data_temp;
    table_transferred = [table_transferred;data_transferred];
end
writetable(table_transferred,[pathname filename_output]);