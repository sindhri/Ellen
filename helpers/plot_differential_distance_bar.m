function plot_differential_distance_bar(pathname, filename)

distance_table = readtable([pathname filename]);
distance = distance_table.differential_distance;
labels = distance_table.drug;

limit = 100;
title_name = 'largest_differential_distance';
if length(distance) > limit
    title_name = ['top ' int2str(limit) ' ' title_name];
else
    limit = length(distance);
end

figure;
barh(distance(1:limit));
yticklabels(labels(1:limit));
yticks(1:limit);
xticks(0:5:max(distance(1:limit)));

title(title_name, 'interpreter', 'none', 'fontsize', 15);

fig_dir = 'fig/';
if exist(fig_dir,'dir')~=7
    mkdir(fig_dir);
end

set(gcf, 'PaperPosition', [0 0 6 10]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [6 10]); 
set(gca,'Ydir','reverse', 'fontsize', 12, ...
    'xaxisLocation','top', 'TickLabelInterpreter', 'none');
saveas(gcf,[pathname fig_dir 'bar_' title_name],'pdf');
close(gcf);

end