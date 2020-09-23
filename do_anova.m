%try use the difference wave of pre and post, then 3-way anova
%only use the mean across the fish in that geno across time

function do_anova(output)
    ncolumn_to_skip = 2; %geno, GroupCount
    pre = output.pre_geno_mean{:,ncolumn_to_skip+1:end};
    post = output.post_geno_mean{:,ncolumn_to_skip+1:end};
    diff = post - pre;
    
    activity_names = output.pre_geno_mean.Properties.VariableNames(ncolumn_to_skip+1:end);
    geno_table = output.pre_geno_mean(:,1);
    for i = 1:size(geno_table,1)
        temp = split(geno_table.geno{i}, '+');
        geno_table.factor1{i} = regexprep(temp{1},' ','');
        geno_table.factor2{i} = regexprep(temp{2},' ','');
        geno_table.factor3{i} = regexprep(temp{3},' ','');
    end
    
    factor1 = geno_table.factor1;
    factor2 = geno_table.factor2;
    factor3 = geno_table.factor3;
    for i = 1:length(activity_names)
        activity_diff = diff(:,i);
        % 3-way anova with interaction
        [p,tbl,stats,terms] = anovan(activity_diff,{factor1,factor2,factor3},...
            'model','interaction','varnames',{'HOM_WT','MC_DMSO','H2O_PTZ'});
        % 3-way anova withou interaction
       %[p,tbl,stats,terms] = anovan(data,{factor1,factor2,factor3});
    end
        
end