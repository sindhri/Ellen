%input zscore file, average across fish and date for each geno
%output two files, 1, geno x activyty, 2, geno x activity summary

%select the file 
%and set up name for the output file by adding '_averaged' to
%the filename
function average_after_zscore
[filename,pathname] = uigetfile('*.csv','select the zscore file');
filename_noext = strsplit(filename,'.');
%remove the '.csv' extention from the original filename
%so that we can add '_averaged' in it.
filename_noext = filename_noext{1};

filename_output_noaggregration = [filename_noext '_mean_by_geno.csv'];

main_table = readtable([pathname filename]);
fprintf('Reading file: %s%s\n', pathname, filename);


parameters = main_table.Properties.VariableNames(4:end-3);
mean_function = @(x) mean(x,'omitnan');

MeanByGenotype = varfun(mean_function,main_table,...
    'GroupingVariables',{'genotype'},...
    'InputVariables',{parameters{:}});
parameters2 = MeanByGenotype.Properties.VariableNames(3:end);


%remove columns that end with 1 or 4 in the parameters
fprintf('Remove columns that end with 1 or 4 in the parameters\n');
to_remove = cell(1);
m = 1;
for i = 1:length(parameters2)
    parameter = parameters2{i};
    ending_letter = parameter(end);
    temp = split(parameter, '_');
    ending_letter2 = temp{3}(end);
    if(strcmp(ending_letter,'1') || strcmp(ending_letter,'4')) || (strcmp(ending_letter2,'1') || strcmp(ending_letter2,'4'))
        to_remove{m} = parameter;
        m = m + 1;
    end
end

MeanByGenotype_trimmed = removevars(MeanByGenotype,to_remove);

fprintf('mean by geno zscore file generated: %s%s\n', pathname, filename_output_noaggregration);
writetable(MeanByGenotype_trimmed,[pathname filename_output_noaggregration]);

%run RMS separately for HOM and HET
HOM = get_target_rows_start(MeanByGenotype_trimmed, 'genotype', 'HOM');
HET = get_target_rows_start(MeanByGenotype_trimmed, 'genotype', 'HET');
filename_output = [pathname filename_noext '_rms_HOM.csv'];
make_rms(HOM, filename_output);
filename_output = [pathname filename_noext '_rms_HET.csv'];
make_rms(HET, filename_output);
end

function final_T = make_rms(T, filename_output)
fprintf('\nCalculating the averages for %s.....\n', filename_output);
%only export the aggregated file for nonburst zscores
if ~contains(filename_output,'burst')
%only select certain columns to caculate rms
%after replacing Inf in the original 'Individualdata.csv'
%there are no more Inf columns in MeanByGenotype_HOM
%inf_columns = [15, 16, 17, 18, 19, 20, 21, 22];
    selection(1).number = 11:14;
    selection(1).name = 'bout';
    selection(2).number = 3:10;
    selection(2).name = 'activity';
    selection(3).number = 23:26;
    selection(3).name = 'sleep';
    selection(4).number = [3:14, 23:size(T,2)];
    selection(4).name = 'all_remove_sleep_latency_duration';
    selection(5).number = 3:size(T,2);
    selection(5).name = 'all';

%allocate a matrics for rms and mean
    n_selection = length(selection);
    matrics_rms_mean = zeros(size(T,1),n_selection*2);

    for i = 1:n_selection
        subtable = T(:,selection(i).number);
        fprintf('selection %d %s columns: \n',i, selection(i).name);
        for j = 1:length(selection(i).number)
            fprintf('%s\n', T.Properties.VariableNames{selection(i).number(j)});
        end
        fprintf('\n');
        matrics = subtable{:,:};
        matrics_rms_mean(:,i) = rms(matrics,2);
        matrics_rms_mean(:,i + n_selection) = mean(matrics,2);
    end

%add the variable names to the output
%and add the genotype to the table
    T_output = array2table(matrics_rms_mean,'VariableNames',...
        {['rms_' selection(1).name],['rms_' selection(2).name],...
        ['rms_' selection(3).name],['rms_' selection(4).name],...
        ['rms_' selection(5).name],...
        ['mean_' selection(1).name],['mean_' selection(2).name],...
        ['mean_' selection(3).name],['mean_' selection(4).name],...
        ['mean_' selection(5).name]});
    final_T = [T(:,1),T_output];
    writetable(final_T, filename_output);

    fprintf('RMS and Mean zscore file generated: %s\n', filename_output);
end
end