%todo
%1. analysis on HOM alone and WT alone, check the interaction:
% MC_DMSO * H20_PTZ
%2. make the anova more friendly

%20200929, saved tables
%try use the difference wave of pre and post, then 3-way anova
%only use the mean across the fish in that geno across time

function [p,tbl,stats] = do_anova(output)
    ncolumn_to_skip = 2; %geno, GroupCount
    pre = output.pre_geno_mean{:,ncolumn_to_skip+1:end};
    post = output.post_geno_mean{:,ncolumn_to_skip+1:end};
    diff = post - pre;
    
    activity_names = output.pre_geno_mean.Properties.VariableNames(ncolumn_to_skip+1:end);
    geno_table = output.pre_geno_mean(:,1);
    for i = 1:size(geno_table,1)
        %remove space and underscore
        temp = split(geno_table.geno{i}, '+');
        factor1{i,1} = regexprep(temp{1},' ','');
        factor2{i,1} = regexprep(temp{2},' ','');
        factor3{i,1} = regexprep(temp{3},' ','');
        factor1{i,1} = regexprep(factor1{i,1},'_','');
        factor2{i,1} = regexprep(factor2{i,1},'_','');
        factor3{i,1} = regexprep(factor3{i,1},'_','');
    end
    if exist([output.pathname 'anova/'],'dir')~=7
       mkdir([output.pathname 'anova/']);
    end    
    for i = 1:length(activity_names)
        activity_name = activity_names{i};
        activity_diff = diff(:,i);
        % 3-way anova with interaction
        %geno-- HOM vs. WT
        %drug1 MC vs. DMSO
        %drug2 H20 vs. PTZ
        [p,tbl,stats] = anovan(activity_diff,{factor1,factor2,factor3},...
            'model','interaction','varnames',{'Geno','drug1','drug2'});
        figure1 = figure(1);
        set(figure1, 'PaperPosition', [0 0 8 4],'PaperPositionMode','auto');
        set(figure1,'PaperSize',[8,4]);
        saveas(figure1,[output.pathname 'anova/anova_' activity_name],'pdf');
        close(figure1);
        
        % multicompare interactive plot
        %multcompare(stats,'Dimension',[1 2 3]);
        %title(activity_name);
        
        % 3-way anova without interaction
        %[p,tbl,stats,terms] = anovan(data,{factor1,factor2,factor3});
    end
        
end