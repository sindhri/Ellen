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
if exist('C','var')
    times = C.time;
else
    times = A.time;
end

for i = 1:length(all_variables)
    variable_name = all_variables{i};
    if contains(variable_name, ['_' parameter_name])
        temp = split(variable_name,'_');
        geno = temp{1};
        break;
    end
end

if strcmp(geno, 'transhet')~=1
    type_names = {'wt', 'het', 'hom'};
else
    type_names = {'wt', 'het_del44', 'het_del5', 'hom'};
end

data = eval([geno '_' parameter_name '_wt']);
data_mean = zeros(size(data,1), length(type_names));
data_sem = zeros(size(data,1), length(type_names));

for i = 1:length(type_names)
    data = eval([geno '_' parameter_name '_' type_names{i}]);
    data_mean(:,i) = mean(data,2,'omitnan');    
    data_sem(:,i) = std(data,0, 2,'omitnan')/sqrt(size(data,2));            
end
    make_one_plot(geno, data_mean, times, data_sem,...
        pathname, type_names, parameter_name);
    zoomin = length(times)-600:length(times);
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
        color_lib = {[0, 0, 1],[1, 0, 0],[1,165/255,0]};
    else
        color_lib = {[0, 0, 1],[1, 0, 0],[1,165/255,0]};
%        color_lib = {[201,229,248]/255, ...
%            [255,222,220]/255, [235,201,129]/255};
    end
else
    if ~isempty(zoomin)
        color_lib = {[0, 0, 1],[1, 0, 0],[1, 165/255,0],[0, 1, 0]};
    else
        color_lib = {[0, 0, 1],[1, 0, 0],[1, 165/255,0],[0, 1, 0]};
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
    if ~isempty(zoomin)
        y = movmean(y, 50);
    end
    x = times;
    std_dev = data_std(:,i)';
    curve1 = y + std_dev;
    curve2 = y - std_dev;
    x2 = [x, fliplr(x)];
    inBetween = [curve1, fliplr(curve2)];
    fill(x2, inBetween, color_lib{i}, 'facealpha',0.1, 'edgecolor', color_lib{i});
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