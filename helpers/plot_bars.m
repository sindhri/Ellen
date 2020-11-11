%make bar plot for all the effects for a 3-way ANVOVA

function plot_bars(output)
    pre = output.pre_fish_mean;
    post = output.post_fish_mean;
    
    ncolumn_to_skip = 3; %skip fish, geno, GroupCount
    activity_names = pre.Properties.VariableNames(ncolumn_to_skip+1:end);
    
    %make a difference table, diff = post - pre
    diff_data = post{:,ncolumn_to_skip+1:end} - pre{:,ncolumn_to_skip+1:end};
    diff = pre;
    diff{:,ncolumn_to_skip+1:end} = diff_data;
    
    %organize the data to reflect the 3-way factors
    factor_names = {'muta','drug1','drug2'};
    
    geno_table = pre(:,2);
    for i = 1:size(geno_table,1)
        %remove space and underscore
        temp = split(geno_table.geno{i}, '+');
        for j = 1:length(temp)
            factor{i,j} = regexprep(temp{j},' ','');
            factor{i,j} = regexprep(factor{i,j},'_','');
        end
    end
    for i = 1:size(factor,2)
        pre.(factor_names{i})=factor(:,i);
        post.(factor_names{i})=factor(:,i);
        diff.(factor_names{i})=factor(:,i);
    end
    
    plotdata_pre = prepare_plotdata(pre, factor_names, activity_names,'pre');
    plotdata_post = prepare_plotdata(post, factor_names, activity_names,'post');
    plotdata_diff = prepare_plotdata(diff, factor_names, activity_names,'diff');
    color_lab = {[0,1,0],[0,0,0.5],[1,0.4,0]};
    for i = 1:length(activity_names)
        make_bars_one_figure(plotdata_pre,output.pathname,color_lab{1},i);
        make_bars_one_figure(plotdata_post,output.pathname,color_lab{2},i);
        make_bars_one_figure(plotdata_diff,output.pathname,color_lab{3},i);
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

function plotdata = prepare_plotdata(data_table, factor_names,...
    activity_names,data_name)
    plotdata.data_name = data_name;
    
    %prepare mean and std for level 1: main effect
    [plotdata.level1m1, plotdata.level1s1] = aggr_func(data_table,factor_names{1},activity_names);
    [plotdata.level1m2, plotdata.level1s2] = aggr_func(data_table,factor_names{2},activity_names);
    [plotdata.level1m3, plotdata.level1s3] = aggr_func(data_table,factor_names{3},activity_names);

    %prepare mean and std for level 2: interaction between 2 factors
    [plotdata.level2m1, plotdata.level2s1] = aggr_func(data_table,{factor_names{1},factor_names{2}},activity_names);
    [plotdata.level2m2, plotdata.level2s2] = aggr_func(data_table,{factor_names{1},factor_names{3}},activity_names);
    [plotdata.level2m3, plotdata.level2s3] = aggr_func(data_table,{factor_names{2},factor_names{3}},activity_names);

    %prepare mean and std for level 3: interaction between 3 factors
    [plotdata.level3m1, plotdata.level3s1] = aggr_func(data_table,factor_names,activity_names);

end

function make_bars_one_figure(plotdata,pathname,color,param_index)
    i=param_index;
    figure;
    set(gcf, 'PaperPosition', [0 0 15 8]); %Position plot at left hand corner with width 5 and height 5.
    set(gcf, 'PaperSize', [15 8]);
    n_column_to_skip = [2,3,4]; %skip different number of variables for different levels
    titlename = [plotdata.data_name '_' plotdata.level1m1.Properties.VariableNames{n_column_to_skip(1)+i}];
    sgtitle(titlename, 'interpreter','none');
    
    % set up the grid
%    pos1 = get(subplot(3,3,1),'Position');
    pos1 = [0.1300    0.7093    0.2134    0.2157];
%    pos2 = get(subplot(3,3,2),'Position');
    pos2 = [0.4108    0.7093    0.2134    0.2157];
%    pos3 = get(subplot(3,3,3),'Position');
    pos3 = [0.6916    0.7093    0.2134    0.2157];
    %pos4 = get(subplot(3,3,4),'Position');
    pos4 = [0.1300    0.4096    0.2134    0.2157];
    %pos5 = get(subplot(3,3,5),'Position');
    pos5 = [0.4108    0.4096    0.2134    0.2157];
    %pos6 = get(subplot(3,3,6),'Position');
    pos6 = [0.6916    0.4096    0.2134    0.2157];
    %pos7 = get(subplot(3,3,[7,8]),'Position');
    pos7 = [0.1300    0.1100    0.4942    0.2157];
    pos1 = adjust_col1(pos1);
    pos2 = adjust_col2(pos2);
    pos3 = adjust_col3(pos3);
    pos4 = adjust_col1(pos4);
    pos5 = adjust_col2(pos5);
    pos6 = adjust_col3(pos6);
    pos7 = adjust_col1(pos7);
    
    
   %plot level 3
    hs = axes;
    set(hs,'OuterPosition',pos7);
    make_one_bar(plotdata.level3m1,plotdata.level3s1,...
        n_column_to_skip(3)+i,3,color);
    yl = ylim;
    
    %plot level 1
    hs = axes;
    set(hs,'OuterPosition',pos1);
    make_one_bar(plotdata.level1m1,plotdata.level1s1,...
        n_column_to_skip(1)+i,1,color);
    ylim(yl);
    
    hs = axes;
    set(hs,'OuterPosition',pos2);
    make_one_bar(plotdata.level1m2,plotdata.level1s2,...
        n_column_to_skip(1)+i,1,color);
    ylim(yl);
    
    hs = axes;
    set(hs,'OuterPosition',pos3);
    pause(1);
    make_one_bar(plotdata.level1m3,plotdata.level1s3,...
        n_column_to_skip(1)+i,1,color);
    ylim(yl);

    %plot level 2
    hs = axes;
    set(hs,'OuterPosition',pos4);
    make_one_bar(plotdata.level2m1,plotdata.level2s1,...
        n_column_to_skip(2)+i,2,color);
    ylim(yl);

    hs = axes;
    set(hs,'OuterPosition',pos5);
    make_one_bar(plotdata.level2m2,plotdata.level2s2,...
        n_column_to_skip(2)+i,2,color);
    ylim(yl);

    hs = axes;
    set(hs,'OuterPosition',pos6);
    make_one_bar(plotdata.level2m3,plotdata.level2s3,...
        n_column_to_skip(2)+i,2,color);
    ylim(yl);

    %save the figures
    if exist([pathname 'bar/'],'dir')~=7
       mkdir([pathname 'bar/']);
    end
    saveas(gcf,[pathname 'bar/' titlename],'pdf');
    saveas(gcf,[pathname 'bar/' titlename],'fig');
    close;

end

function make_one_bar(datam,datas,col_index,level,color)
    n_bar = size(datam,1);
    x = 1:n_bar;
    m = datam{:,col_index};
    s = datas{:,col_index};    
    bar(x,m,'FaceColor',color);
    hold on;
    er = errorbar(x,m,s,s);    
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';

    if level==1
        labelArray = datam{:,1:level};
        bar_labels = labelArray;
    end
    labelArray = datam{:,1:level}';
    if level==2
        bar_labels = strtrim(sprintf('%s\\newline%s\n', labelArray{:}));
    end
    if level==3
        labelArray = datam{:,1:level}';
        bar_labels = strtrim(sprintf('%s\\newline%s\\newline%s\n', labelArray{:}));
    end
    set(gca,'XTickLabel',bar_labels,'TickLabelInterpreter','tex');
    
    hold off
end

function pos = adjust_col1(pos)
    pos(1) = pos(1) - 0.1;
    pos(3) = pos(3) * 1.75;
end

function pos = adjust_col2(pos)
    pos(1) = pos(1) - 0.05;
    pos(3) = pos(3) * 1.75;
end

function pos = adjust_col3(pos)
    pos(1) = pos(1);
    pos(3) = pos(3) * 1.75;
end