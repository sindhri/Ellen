% 20100110
% select a combined .mat file
% plot time series for wildtype (+/+) in blue, HET(geno/+) in green, HOM(geno/geno) in red
% For transhet, it has 2 HET (geno1/+, geno2/+) plot in green and
% yellow
function plot_time_series
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
    if contains(variable_name, '_data')
        temp = split(variable_name,'_');
        geno = temp{1};
        break;
    end
end

if strcmp(geno, 'transhet')~=1
    data = eval([geno '_data_wt_mean']);
    data(:,2) = eval([geno '_data_het_mean']);
    data(:,3) = eval([geno '_data_hom_mean']);

    data_std = zeros(size(data));
    data_std(:,1) = eval(['std(' geno '_data_wt,0,2)']);
    data_std(:,2) = eval(['std(' geno '_data_het,0,2)']);
    data_std(:,3) = eval(['std(' geno '_data_hom,0,2)']);

    color_lib = {'b','g','r'};
    legend_names = {'wt', 'het', 'hom'};
else
    data = eval([geno '_data_wt_mean']);
    data(:,2) = eval([geno '_data_het_del44_mean']);
    data(:,3) = eval([geno '_data_het_del5_mean']);
    data(:,4) = eval([geno '_data_hom_mean']);

    data_std = zeros(size(data));
    data_std(:,1) = eval(['std(' geno '_data_wt,0,2)']);
    data_std(:,2) = eval(['std(' geno '_data_het_del44,0,2)']);
    data_std(:,3) = eval(['std(' geno '_data_het_del5,0,2)']);
    data_std(:,4) = eval(['std(' geno '_data_hom,0,2)']);

    color_lib = {'b','g','y','r'};  
    legend_names = {'wt', 'het_del44','het_del5', 'hom'};
end
    make_one_plot(geno, data, times, data_std,...
        color_lib, pathname, legend_names);
    zoomin = length(times)-600:length(times);
    make_one_plot(geno, data, times, data_std,...
        color_lib, pathname, legend_names, zoomin);
end

function make_one_plot(geno, data, times, data_std, ...
    color_lib, pathname, legend_names, zoomin)
if nargin==8
    data = data(zoomin,:);
    times = times(zoomin);
    data_std = data_std(zoomin,:);
    title_name = [geno ' zoomed in'];
else
    title_name = geno;
end
figure;
set(gcf, 'PaperPosition', [0 0 5 5]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [5 5]);

p = [];
for i = 1:size(data,2)
    y = data(:,i);
    x = times;
    std_dev = data_std(:,i);
    curve1 = y + std_dev;
    curve2 = y - std_dev;
    x2 = [x, fliplr(x)];
    inBetween = [curve1; fliplr(curve2)];
    fill(x2, inBetween, color_lib{i});
    hold on;
    p(i) = plot(x, y, color_lib{i}, 'LineWidth', 2);
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