% 20210524, added target_geno, so the code can be used to generate either
% HOM or HET distance, compared to WT+DMSO
% added target_geno and whether to exclude sleep latency and length into
% filename
% organize results in subfolders
% refactorized
% 20200430, changed PCA to only include HOM + drugs and WT+DMSO
% changed 24 to nparams, able to limit sleep latency and length
% 20210402, put the normalization back
% 20210326, pick onlly HOMM and WT
% 20200325, removed the normalization step, added optional inputs
% 20210103, closed all the figures at the end
% 20201203, places the figures in a folder 'fig' in the source file folder
% input, geno x activity file
% output, graphs, euclidean distance
function filename_output = run_PCA_euclidean(target_geno, pathname, filename, remove_sleep_latency_length)

if nargin==1
    [filename,pathname] = uigetfile('*.csv','select the mean_by_geno file');
    remove_sleep_latency_length = 0;
end
main_table = readtable([pathname filename]);

%only use HOM
%main_table = main_table(contains(main_table.genotype, 'HOM') | contains(main_table.genotype, 'WT'),:);
main_table = main_table(contains(main_table.genotype, target_geno) | contains(main_table.genotype, 'WT + DMSO'),:);
parameters = main_table.Properties.VariableNames;
nparams = length(parameters)-2;%exclude genotype and GroupCount
% remove sleep latency and sleep length in the parameters

if remove_sleep_latency_length == 1
    subfolder_name = [target_geno '_sleep_LL_excluded'];
    for i = 3:length(parameters) %skip genotype and GroupCount
        param = parameters(i);
        if contains(param, 'sleepLatency') || contains(param, 'sleepLength')
            main_table = removevars(main_table, param);
            nparams = nparams - 1;
        end
    end
else
    subfolder_name = [target_geno '_all_params'];
end
filename_noext = strsplit(filename,'.');
filename_noext = filename_noext{1};
filename_output = [filename_noext '_distance_' subfolder_name '.csv'];

%skip the first two columns, genotype and GroupCount
data = main_table{:,3:end}; 

%data, 52 x nparams, geno_by_activity
%original program 1st dimension is 120, so it is likely to be geno instead
%of activity
%which way makes more sense? to normalize across geno, or activity?
%If normalize across geno, it would compare the relative strength of geno
%for each individual activity, 
%if normalize across activity, it would compare the relative strength of
%activity, for each individual geno
%I think normalize across geno makes more sense.
%However in papers the plotting has the dimeison of geno x activity
%On the other hand, both papers seems to trying to
%do PCA on the dimension of activity, which would have activity as column
%and study EACH geno, how much it is relfected on a bundle of activites
%I guess the assumption is that activities could be highly correlated
%so trying to find the bundle of activities that best representing the
%whole picture is beneficial
%Then study EACH geno, (instead of bundle geno into groups)


%normalize across activity (2nd dimension), for each geno (1st dimension)
[n_geno, n_activity] = size(data);
for i = 1:n_geno
  logdata(i,:) = data(i,:)-mean(data(i,:)); % make sure there is no DC offset)
  normdata(i,:) = logdata(i,:)./std(logdata(i,:)); 
end

fig_dir = [pathname '/' subfolder_name '/fig/'];
if exist(fig_dir, 'dir') ~=7
    mkdir(fig_dir);
end

fullproj = PCA_plots(normdata, fig_dir, nparams, n_geno);
filename_output = create_distance(fullproj, main_table, pathname, subfolder_name, filename_output);
plot_distance(target_geno, [pathname subfolder_name], filename_output);

end

function filename_output = create_distance(fullproj, main_table, pathname, subfolder_name, filename_output)
euclidean_distances = pdist(fullproj);
formated_distances = squareform(euclidean_distances);
distance_table = array2table(formated_distances);
distance_table.Properties.VariableNames = main_table.genotype;

distance_table = [main_table.genotype, distance_table];
distance_table.Properties.VariableNames{1} = 'euclidean_distance';
distance_file = [pathname subfolder_name '/' filename_output];
writetable(distance_table, distance_file);
fprintf('Euclidean distance file generated: %s%s\n', [pathname subfolder_name], filename_output);
end

function fullproj = PCA_plots(normdata, fig_dir, nparams, n_geno)
% decompose the data matrix (i.e. do principal components analysis
[u fulls fullv] = svd(normdata);

% see how much of the variance we explain (i.e. a scree plot)
tmp = 100*(diag(fulls)./sum(diag(fulls)));

% and get the projections-- but I will try the first 10
fullproj = normdata*fullv(:,1:10);

% and the vector lengths
fulllen = sqrt(fullproj(:,1).^2 + fullproj(:,2).^2 + fullproj(:,3).^2+ fullproj(:,4).^2 + fullproj(:,5).^2 + fullproj(:,6).^2+ fullproj(:,7).^2 + fullproj(:,8).^2 + fullproj(:,9).^2+ fullproj(:,10).^2);

% a figure describing the raw data
figure
load my_colormap;
imagesc(normdata,[-3 3])
set(gca,'position',[0 0 1 1],'visible','off')
set(gcf,'position',[100 100 800 350],'paperposition',[1 1 8 3.5])
colormap(my_colormap)
print(gcf,'-depsc2',[fig_dir 'woods1.eps'])
% better to annotavte/reshape this in illustrator


figure
hold on
plot(1,tmp(1),'r.','markersize',18)
plot(2,tmp(2),'g.','markersize',18)
plot(3,tmp(3),'b.','markersize',18)
plot(4:length(tmp),tmp(4:end),'k.','markersize',18)
plot([0 10],[10 10],'k:')
axis square,box off
set(gca,'xlim',[0 10],'ylim',[0 60],'ytick',[0 20 40],'xtick',[],'fontsize',12)
xlabel('Eigenvector')
ylabel('% Variance')
%stevify
print(gcf,'-depsc2',[fig_dir 'woods2.eps'])

% Does NOT explain that much of the data? 
% (83.16% for the first 10 factors, how much is enough

% now, how did the first three do?
sprintf('The first three components account for %g of the variance',sum(tmp(1:3)))
%'The first three components account for 44.3158 of the variance'
% Jia: that's not a lot of variance

% ok, now let's plot the first three eigenvectors
figure
hold on
%what is order? why only starting at 10th?
%plot(fullv(order,10),'r','linewidth',2);
%plot(fullv(order,11),'g','linewidth',2);
%plot(fullv(order,12),'b','linewidth',2);
plot(fullv(1:nparams,10),'r','linewidth',2);
plot(fullv(1:nparams,11),'g','linewidth',2);
plot(fullv(1:nparams,12),'b','linewidth',2);
plot([0 26],[0 0],'k:')
plot([8 8],[-.8 .8],'k:')
axis square, box off
set(gca,'xlim',[0 26],'ylim',[-.8 .8],'visible','off')
%stevify
print(gcf,'-depsc2',[fig_dir 'woods3.eps'])
% probably good to annotate the interesting peaks with text detailing what they are
% but that is easier to do in illustrator



% Jia: According to the note, putting geno as the first dimension is the
% correct call
% This is to test the model is stable across all genos, not driven by a few
% genos, it is a sense of bootstrapping

% now, let's get to testing the generality of this space
% if we leave one peptide out, and redo the analysis
% then the spaces should be similar.
% the way we'll test that is to see 
% a) first, are the spaces the same (i.e. are the vectors similar?)
% b) is each peptide similarly far from the origin?
% c) is each peptide similarly far from its neighbor
% we'll use the Euclidian distance (i.e. sqrt of the sum of squares) as a
% distance metric.

for i = 1:n_geno
 dex = circshift(1:n_geno,[0 i]); % choose a unique set for each run
 [u s v] = svd(normdata(dex(1:n_geno-1),:)); % take the first six
 
 % how much variance is accounted for?
 tmp = 100*(diag(s)./sum(diag(s)));
 varacct(i) = sum(tmp(1:10)); 
 
 % how similar are the vectors that define the space?
 eigvecproj(i,:) = diag(v(:,1:10)'*fullv(:,1:10));
 
 % save the indiviudal projections
 skiponev(i,:,:) = v;
 
 % get the new projections
 proj(:,:,i) = normdata*v(:,1:10);
 
 % now, how far is each point from the origin in this new space?
 veclen(i,:) = sqrt(proj(:,1,i).^2 + proj(:,2,i).^2 + proj(:,3,i).^2+proj(:,4,i).^2 + proj(:,5,i).^2 + proj(:,6,i).^2+proj(:,7,i).^2 + proj(:,8,i).^2 + proj(:,9,i).^2+proj(:,10,i).^2);
 
 % and how far is each point from each others?
 spacing(i,:) = pdist(proj(:,:,i));
end


% Jia: not sure what is it trying to do. 
% we have nparams activities so this will do nothing to our results

% on the fourth run, eigenvectors 2 and 3 end up switched. so we'll
% fix that manually
if size(normdata,2) == 16
eigvecproj(4,2) = skiponev(4,:,2)*fullv(:,3);
eigvecproj(4,3) = skiponev(4,:,3)*fullv(:,2);
end

if size(normdata,2) == 17
eigvecproj(5,2) = skiponev(5,:,2)*fullv(:,3);
eigvecproj(5,3) = skiponev(5,:,3)*fullv(:,2);
end

% now, how did the first three do?
sprintf('The first three components in the leave one out analysis\naccount for %g of the variance, SD = %g',mean(varacct),std(varacct))
% The first three components in the leave one out analysis
% account for 83.2439 of the variance, SD = 0.182833
% Jia: how come bootstraping increased so much variance?


numruns = 10000;
randveclen = zeros(numruns,n_geno);
% 1326 = 52* 51 /2, number of pairwise distances between pairs of
% observations

randspacing = zeros(numruns,n_geno * (n_geno-1)/2);
randeigvecproj = zeros(numruns,10);

for i = 1:numruns
  randnormdata = reshape(normdata(randperm(numel(normdata))),[size(normdata,1) size(normdata,2)]);
  [u s v] = svd(randnormdata);
  randeigvecproj(i,:) = diag(v(:,1:10)'*fullv(:,1:10));
  randproj = normdata*v(:,1:10);
  randveclen(i,:) = sqrt(randproj(:,1).^2 + randproj(:,2).^2 + randproj(:,3).^2+randproj(:,4).^2 + randproj(:,5).^2 + randproj(:,6).^2+randproj(:,7).^2 + randproj(:,8).^2 + randproj(:,9).^2+randproj(:,10).^2);
  % and how far is each point from each others?
  randspacing(i,:) = pdist(randproj);
end

% For random data, we'll just shuffle the original data so
% that we use the same values, but without any structure.

% first, the projections
bins = 0:.01:1;
P = hist(abs(eigvecproj(:)),bins);
randP = hist(abs(randeigvecproj(:)),bins);
P = P./sum(P);
randP = randP./sum(randP);

%black 'k' is the ACTUAL egen vector
%gray '.6,.6,.6' is the random vector
figure
hold on
plot(bins,cumsum(P),'k','linewidth',2);
%plot(bins,cumsum(randP),'color',[.6 .6 .6],'linewidth',2)
plot(bins,cumsum(randP),'color',[1, .6, 0],'linewidth',2)
axis square
set(gca,'xlim',[0 1],'ylim',[0 1], 'fontsize',12,'xtick',[0 .5 1],'ytick',[0 .5 1])
xlabel('Similarity of Eigenvectors')
ylabel('Cumulative Probability')
% stevify 
print(gcf,'-depsc2',[fig_dir 'woods4.eps'])



% that seems pretty clear; the vectors are largely the same even when we've
% left out one peptide. now, is each peptide about the same distance from
% the origin as it would have been in the original space?
bins = 0:.1:5;
randy = hist(randveclen(:),bins);
y = hist(veclen(:),bins);

%orange is random vector
%black is left-out model
figure
hold on
%plot(bins,cumsum(randy./sum(randy)),'color',[.6 .6 .6],'linewidth',2)
plot(bins,cumsum(randy./sum(randy)),'color',[1 .6 0],'linewidth',2)
plot(bins,cumsum(y./sum(y)),'color','k','linewidth',2)
axis square
box off
set(gca,'xlim',[0 5],'xtick',[0  5],'ylim',[0 1],'ytick',[0 .5 1],'fontsize',12)
xlabel('Length')
ylabel('Cumulative Probability')
%stevify
print(gcf,'-depsc2',[fig_dir 'woods5.eps'])




% looks pretty good too; now, what about the inter-point distance?
% here, we want to see how each value compares to the actual distance in
% the full space

% first, the random values
randx = repmat(pdist(fullproj),[numruns 1]);
x = repmat(pdist(fullproj),[n_geno 1]);

% there's really no reason to plot all 10000 points
dex = randperm(numel(randx));
dex = dex(1:1000);

figure
hold on
%plot(randx(dex),randspacing(dex),'.','color',[.6 .6 .6])
plot(randx(dex),randspacing(dex),'.','color',[1 .6 0])
plot(x(:),spacing(:),'k.')
plot([0 10],[0 10],'k:')
axis square
box off
set(gca,'xlim',[0 10],'ylim',[0 8],'xtick',[0 10],'ytick',[0 8],'fontsize',12)
xlabel('True Distance')
ylabel('Distance')
% stevify
print(gcf,'-depsc2',[fig_dir 'woods6.eps'])

close all;
end

%20210524, plot distance based on WT+DMSO, input distance file with single
%category of target_geno
%202103, removed separate control for 'related WT' because it was changed
%to 'WT' while creating z scores
%20210103, fixed a bug of overlapping figures when run right after the last
%a script that has a figure open
function plot_distance(target_geno, pathname, filename)

if nargin==1
    [filename,pathname] = uigetfile('*.csv','select the distance file');
end

filename_noext = strsplit(filename,'.');
filename_noext = filename_noext{1};

main_table = readtable([pathname '/' filename]);

fig_dir = [pathname '/fig/'];
if exist(fig_dir, 'dir') ~=7
    mkdir(fig_dir);
end

control_genos_original = {'WT + DMSO'};
% clean up to match the format of matlab columns
control_genos_cleaned = regexprep(control_genos_original, '_', '');
control_genos_cleaned = regexprep(control_genos_cleaned, ' ', '');
control_genos_cleaned = regexprep(control_genos_cleaned, '+', '_');

% total distance to WT+DMSO
for i = 1:length(control_genos_cleaned)
    control_geno = control_genos_cleaned{i};
   
    distance = main_table.(control_geno);
    label = main_table.euclidean_distance;
       
    table_temp = table(label, distance);
    table_sorted = sortrows(table_temp,'distance');
    table_sorted(1,:) = [];
    
    % take only the target rows
    table1 = get_target_rows(table_sorted, target_geno);
    
    writetable(table1,[pathname '/'  filename_noext '_compare_' target_geno '_' control_geno '.csv']);
    plot_distance_bar(table1.distance, table1.label,...
        fig_dir, control_geno, target_geno);
end
end

function plot_distance_bar(distance, labels, ...
    fig_dir, control_geno, target_geno)

limit = 100;
title_name = [target_geno '_to_' control_geno];
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

set(gcf, 'PaperPosition', [0 0 6 10]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [6 10]); 
set(gca,'Ydir','reverse', 'fontsize', 12, ...
    'xaxisLocation','top', 'TickLabelInterpreter', 'none');
saveas(gcf,[fig_dir 'bar_' target_geno '_to_' control_geno],'pdf');
close(gcf);

end

% take only the rows of column :label that start with target
% the table content is a cell, so had to take the nested inforation out
function T_output = get_target_rows(T, target)
m = 1;
selected_label = cell(1);
selected_distance = [];
for i = 1:size(T,1)
    the_label = T.label{i};
    if length(the_label) < length(target)
        continue;
    end
    
    sub_label = extractBetween(the_label, 1, length(target));
    if strcmp(sub_label, target)==1
        selected_label{m,1} = the_label;
        selected_distance(m,1) = T.(T.Properties.VariableNames{2})(i);
        m = m + 1;
    end    
end
T_output = table(selected_label, selected_distance);
T_output.Properties.VariableNames = T.Properties.VariableNames;
end

