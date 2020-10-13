% 20201013
% combine the pre and post table, add time(either pre_ or post), 
% and extract the geno column into three factors

function combined_table = combine_add_time_factors(pre_table, post_table,...
    geno_table_index)

    
    % fpr both pre and post tables, split the geno column by + 
    % and extract the original three levels of
    % factors: geno x drug1 x drug2
    geno_table = pre_table(:,geno_table_index);
    for i = 1:size(geno_table,1)
    %remove space and underscore
        temp = split(geno_table.geno{i}, '+');
        factor1{i,1} = regexprep(temp{1},' |_','');
        factor2{i,1} = regexprep(temp{2},' |_','');
        factor3{i,1} = regexprep(temp{3},' |_','');
    end
    pre_table.factor1 = factor1;
    pre_table.factor2 = factor2; 
    pre_table.factor3 = factor3; 
    % add the time column for pre
    time = cell(1);
    for i = 1:size(geno_table,1)
        time{i,1} = 'pre';
    end
    pre_table.time = time;
    
    post_table.factor1 = factor1;
    post_table.factor2 = factor2; 
    post_table.factor3 = factor3;   
    %add the time column for post
    time = cell(1);
    for i = 1:size(geno_table,1)
        time{i,1} = 'post';
    end
    post_table.time = time;
    
    % combine pre and post table
    combined_table = [pre_table; post_table];

end