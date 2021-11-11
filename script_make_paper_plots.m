% 20210110, made time-series plots
% make a folder of the combined files
% the program will ask you to select which combined file to plot
% then the final time-series plots will be in the same folder as the combined files

addpath('helpers/');
plot_time_series('data');

plot_time_series('avewaking');
% only del5 has tenminute. need to fix bugs
% plot_time_series('tenminute');