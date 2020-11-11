% 20201013, fixed message box multiple line bug
% 20200929, added saving an intermediate file of individual fish aggregated
% across time
% 20200922, removed 'Fun_' from the column names in the aggregated files
% pre and post drug experiment data
% 96 fish, measured 1 time per minute (60 seconds) for a total of 30 minutes
% corresponding geno type for each fish is in the text file '200725_0Bgenotype.txt'
function output = pre_post()
%post = readtable('files/200725_0B_pre.xlsx');
%post = readtable('files/200725_0B_post.xlsx');

%Mac OS is not able to show the title when opening the files
%Just have to remember the sequence of the files
%pre, post, geno
[filename_pre,pathname] = uigetfile('*.xlsx','select the pre file in excel format');
pre = readtable(([pathname filename_pre]));

output.pathname = [pathname 'output/'];
if exist(output.pathname,'dir')~=7
   mkdir(output.pathname);
end
pathname_intermediate = output.pathname;
%set up filename and path for the intermediate export 
%replace pre with intermediate, or add '_intermediate' at the end of the
%file name
if contains(filename_pre,'pre')
    filename_intermediate = regexprep(filename_pre, 'pre', 'intermediate');
else
    filename_intermediate = regexprep(filename_pre,'.xlsx','_intermediate.xlsx');
end

[filename_post,pathname] = uigetfile('*.xlsx','select the post file in excel format');
post = readtable(([pathname filename_post]));

% remove the extra seconds from post after 1800 seconds, so it's the same size as pre
post(post.start>=1800,:)=[];

% add a fish column based on the last two digits of location
pre = add_fish_column(pre);
post = add_fish_column(post);

% import geno file and convert it to a table with columns 'geno' and 'fish'
% mapping
%geno_file = importdata('files/200725_0Bgenotype.txt', '\t', 2);
[filename_geno,pathname] = uigetfile('*.txt','select the geno text file');
geno_file = importdata([pathname, filename_geno], '\t', 2);

fprintf('Preprocess files......\n');
geno_table = import_geno_to_table(geno_file);

% add the geno column to both pre and post dataset
pre_with_geno = innerjoin(geno_table, pre, 'keys','fish');
post_with_geno = innerjoin(geno_table, post, 'keys','fish');

% gather the parameters to measure in this new dataset
parameters = pre_with_geno.Properties.VariableNames(13:21);

output.geno = geno_file.colheaders;
output.parameters = parameters;

% aggregation across fish and time
% The first aggregation is almost unnecerry, because the original file is 
% already the format of one row per fishxstartxgeno 
% but this step removed irrelevant variables and
% sorted the variables by the order of fish,start,geno, 
% so it is consistent with later aggregated files

fprintf('restructure fish, time geno file......\n');
[output.pre_fish_time_mean, output.pre_fish_time_std] = aggr_func(pre_with_geno, ...
    {'fish','start','geno'}, parameters);
[output.post_fish_time_mean, output.post_fish_time_std] = aggr_func(post_with_geno, ...
    {'fish','start','geno'}, parameters);

fprintf('creating intermediate file......\n');
% aggregation across time
[output.pre_fish_mean, output.pre_fish_std] = aggr_func(pre_with_geno, ...
    {'fish','geno'}, parameters);
[output.post_fish_mean, output.post_fish_std] = aggr_func(post_with_geno, ...
    {'fish','geno'}, parameters);
export_fish_intermediate(output.pre_fish_mean,output.post_fish_mean,...
    pathname_intermediate, filename_intermediate);

fprintf('aggregate across fish......\n');
% aggregation across geno and time based on individual fish results
[output.pre_geno_time_mean, output.pre_geno_time_std] = aggr_func(output.pre_fish_time_mean, ...
    {'geno','start'}, parameters);
[output.post_geno_time_mean, output.post_geno_time_std] = aggr_func(output.post_fish_time_mean, ...
    {'geno','start'}, parameters);

% aggregation across geno only based on individual fish results
fprintf('aggregate across fish and time......\n');
[output.pre_geno_mean, output.pre_geno_std] = aggr_func(output.pre_fish_mean, ...
    {'geno'}, parameters);
[output.post_geno_mean, output.post_geno_std] = aggr_func(output.post_fish_mean, ...
    {'geno'}, parameters);

fprintf('Output generated for plotting and stats\n');
end

function data_table = add_fish_column(data_table)
% add a fish column to both the pre and post data
for i = 1:size(data_table,1)
    temp = data_table.location{i};
    %location W001 transfers to 01
    %location W010 transfers to 10
    fish_temp = temp(end-1:end);
    %convert to integer
    data_table.fish(i) = str2double(fish_temp);
end
end

function geno_table = import_geno_to_table(geno_file)
% convert geno text file into a table
geno_table = table();
for i = 1:size(geno_file.data,2)
    fish = geno_file.data(:,i);
    %remove NaN in the fish column
    fish(find(isnan(fish)==1))=[];
    geno = repmat(geno_file.colheaders(i),size(fish));
    geno_table = [geno_table; table(geno, fish)];
end
end

function [tablem, tables] = aggr_func(data_table, ...
    aggr_var, parameters)

% create metrics
mean_function = @(x) mean(x,'omitnan');
std_function = @(x) std(x,'omitnan');

tablem = varfun(mean_function,data_table,...
    'GroupingVariables',aggr_var,...
    'InputVariables',{parameters{:}});

tables = varfun(std_function,data_table,...
    'GroupingVariables',aggr_var,...
    'InputVariables',{parameters{:}});

%remove 'Fun_' in the columns of the aggregated tables
tablem.Properties.VariableNames=regexprep(tablem.Properties.VariableNames, 'Fun_', '');
tables.Properties.VariableNames=regexprep(tables.Properties.VariableNames, 'Fun_', '');

end

function export_fish_intermediate(pre,post,pathname,filename)
%matlab requires the column to have the same number of
%digits in order to concatinate the tables    
time = repmat('pre_',[size(pre,1),1]);
    pre = addvars(pre,time,'Before','geno');
    time = repmat('post',[size(post,1),1]);
    post = addvars(post,time,'Before','geno');
    final = [pre;post];
    writetable(final,[pathname filename]);
    msgbox(sprintf(['Intermediate file saved at\n' pathname filename]));
end
