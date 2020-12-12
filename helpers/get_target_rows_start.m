% take only the rows on colname that start with target
function T_output = get_target_rows_start(T, colname, target)
selected_rows = [];
for i = 1:size(T,1)
    the_label = T.(colname){i};
    if length(the_label) < length(target)
        continue;
    end
    
    sub_label = extractBetween(the_label, 1, length(target));
    if strcmp(sub_label, target)==1
        selected_rows = [selected_rows, i];
    end    
end
T_output = T(selected_rows,:);
end