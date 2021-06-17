%20210617, can only be used for plotting differential distance or
%differential distance avg
% added sortrows to sort in descending order
function plot_differential_distance_bar(pathname, filename)

distance_table = readtable([pathname filename]);
Exist_Column = strcmp('differential_distance',distance_table.Properties.VariableNames);
val = Exist_Column(Exist_Column==1) ;
if ~isempty(val)
    col_name = 'differential_distance';
    figure_name = 'Differential Drug Effects by Doses';
else
    col_name = 'differential_distance_avg';
    figure_name = 'Differential Drug Effects';
end
distance_table = sortrows(distance_table, col_name,'descend');
distance = distance_table.(col_name);
labels = distance_table.drug;

limit = 100;
%title_name = 'largest_differential_distance';
if length(distance) > limit
    figure_name = ['top ' int2str(limit) ' ' figure_name];
else
    limit = length(distance);
end

figure;
barh(distance(1:limit));
yticklabels(labels(1:limit));
yticks(1:limit);
xticks(0:5:max(distance(1:limit)));

title(figure_name, 'interpreter', 'none', 'fontsize', 15);

fig_dir = 'fig/';
if exist(fig_dir,'dir')~=7
    mkdir(fig_dir);
end

set(gcf, 'PaperPosition', [0 0 6 10]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [6 10]); 
set(gca,'Ydir','reverse', 'fontsize', 12, ...
    'xaxisLocation','top', 'TickLabelInterpreter', 'none');
saveas(gcf,[pathname fig_dir 'bar_' figure_name],'pdf');
close(gcf);

end