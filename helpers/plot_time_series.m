% 20210124, added parameter_name, which could be either 'avewaking', or
% 'tenminute'. Also calcualte mean instead of using the _mean, because or
% avewaking and tenminute the _mean variables don't exist.
% 20210110
% select a combined .mat file
% plot time series for wildtype (+/+) in blue, HET(geno/+) in green, HOM(geno/geno) in red
% For transhet, it has 2 HET (geno1/+, geno2/+) plot in green and
% yellow
function plot_time_series(parameter_name)
% select the combined mat file
[filename,pathname] = uigetfile('*.mat','select the combined mat file');
load([pathname filename]);
all_variables = who;
times = [];
data_type = ''; %to distinguish old data naming convention 
%in scn1lab_del44 and the rest of the combined files 
%from the 2021/11 batch
for i = 1:length(all_variables)
    varname = all_variables{i};
    if length(varname) >=3
        if strcmp(varname(1:3), 'exp')==1
            times = eval(varname).time;
            data_type = 'new';
        end
    end
end
if isempty(times) % for scn1lab_del44
    data_type = 'old';
    if exist('C','var')
        times = C.time;
    else
        times = A.time;
    end
end

for i = 1:length(all_variables)
    variable_name = all_variables{i};
    if contains(variable_name, ['_' parameter_name])
        temp = split(variable_name,'_');
        geno = temp{1};
        break;
    else
        if contains(variable_name, [parameter_name '_']) 
           temp = split(variable_name,'_');
           geno = temp{3};
           break;
        end
    end
end

if strcmp(data_type,'old')==1
    if strcmp(geno, 'transhet')~=1
        type_names = {'wt', 'het', 'hom'};
    else
        type_names = {'wt', 'het_del44', 'het_del5', 'hom'};
    end
else
    type_names = {'WT', 'HET', 'HOM'};
end

%initiate the mean matrix
varname = [geno '_' parameter_name '_' type_names{1}];
if exist(varname,'var')
    data = eval([geno '_' parameter_name '_' type_names{1}]);
else %to accommodate new naming convention
    if strcmp(geno, 'chd8') ~=1
        data = eval([parameter_name '_' type_names{1} '_' geno]);
    else
        data = eval([parameter_name '_' type_names{1} '_' geno '_merged']);
    end
end
data_mean = zeros(size(data,1), length(type_names));
data_sem = zeros(size(data,1), length(type_names));

for i = 1:length(type_names)
    varname = [geno '_' parameter_name '_' type_names{i}];
    if exist(varname,'var')
        data = eval([geno '_' parameter_name '_' type_names{i}]);
    else %to accommodate new naming convention
    if strcmp(geno, 'chd8') ~=1
        data = eval([parameter_name '_' type_names{i} '_' geno]);
    else
        data = eval([parameter_name '_' type_names{i} '_' geno '_merged']);
    end
    end
    data_mean(:,i) = mean(data,2,'omitnan');    
    data_sem(:,i) = std(data,0, 2,'omitnan')/sqrt(size(data,2));            
end
%make_one_plot(geno, data_mean, times, data_sem,...
%        pathname, type_names, parameter_name);
if strcmp(geno, 'del44') == 1
        zoomin = 2000:4000; % roughly hour 35-65
else
    zoomin = 800:length(times);
end
make_one_plot(geno, data_mean, times, data_sem,...
        pathname, type_names, parameter_name, zoomin);
end

function make_one_plot(geno, data, times, data_std, ...
    pathname, legend_names, parameter_name, zoomin)
if nargin==8
    data = data(zoomin,:);
    times = times(zoomin);
    data_std = data_std(zoomin,:);
    title_name = [geno '_' parameter_name ' zoomed in'];
else
    zoomin = '';
    title_name = [geno '_' parameter_name];
end

if size(data,2)==3
    if ~isempty(zoomin)
%        color_lib = {[0, 0, 1],[0, 1, 0],[1, 0, 0]};
        color_lib = {[24, 110, 181]/255,...
            [18, 99, 26]/255,...
            [196, 79, 6]/255};
    else
%        color_lib = {[0, 0, 1],[0, 1, 0],[1, 0, 0]};
        color_lib = {[24, 110, 181]/255,...
            [18, 99, 26]/255,...
            [196, 79, 6]/255};
    end
else
    if ~isempty(zoomin)
        color_lib = {[0, 0, 1],[0, 1, 0],[1, 165/255,0],[1, 0, 0]};
    else
        color_lib = {[0, 0, 1],[0, 1, 0],[1, 165/255,0],[1, 0, 0]};
%        color_lib = {[211/255,223/255,242/255], ...
%            [242/255,204/255,196/255], [241/255,242/255,196/255],...
%            [242/255,235/255,196/255]};
    end
end
figure;
set(gcf, 'PaperPosition', [0 0 10 5]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [10 5]);

p = [];
for i = 1:size(data,2)
    y = data(:,i)'; %important to transpose!
%    if ~isempty(zoomin)
        y = movmean(y, 50);
%    end
    x = times;
    std_dev = data_std(:,i)';
    curve1 = y + std_dev;
    curve2 = y - std_dev;
    x2 = [x, fliplr(x)];
    inBetween = [curve1, fliplr(curve2)];
    fill(x2, inBetween, color_lib{i}, 'facealpha',0.1, 'edgecolor', [1,1,1]);
    hold on;
    p(i) = plot(x, y, 'color',color_lib{i}, 'LineWidth', 1);
end
if i==3
    l = legend([p(1), p(2), p(3)],legend_names);
else
    l = legend([p(1), p(2), p(3), p(4)],legend_names);
end
set(l, 'interpreter','none');
title(title_name, 'interpreter','none');
xlabel('hours');
xlim([times(1),times(end)]);
ylabel('data');

figure_dir = [pathname 'time_series_plots' filesep];
%save the figures
if exist(figure_dir,'dir')~=7
   mkdir(figure_dir);
end
saveas(gcf,[figure_dir title_name],'pdf');
saveas(gcf,[figure_dir title_name],'fig');
close;

end