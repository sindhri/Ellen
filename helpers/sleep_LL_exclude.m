%20210617, remove sleeplatency and sleeplength parameters
function main_table = sleep_LL_exclude(main_table)
    parameters = main_table.Properties.VariableNames;
    for i = 1:length(parameters) %skip genotype and GroupCount
        param = parameters(i);
        if contains(param, 'sleepLatency') || contains(param, 'sleepLength')
            main_table = removevars(main_table, param);
        end
    end
    
end
