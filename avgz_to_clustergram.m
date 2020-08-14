[filename,pathname] = uigetfile('*.csv','select the avg zscore file');
avgz = readtable([pathname filename]);

data = avgz{:,2:end};
genotypes = avgz{:,1};
parameters = avgz.Properties.VariableNames(2:end);

%replaced underscore with space in the row and column names
genotypes = strrep(genotypes,'_', ' ');
parameters = strrep(parameters,'_', ' ');

load my_colormap;
cgz = clustergram(data, 'RowLabels', genotypes,...
                             'ColumnLabels', parameters,...
                             'RowPdist', 'correlation',...
                             'ColumnPdist', 'correlation',...
                             'ImputeFun', @knnimpute)
cgz.Colormap = my_colormap;           
cgz.DisplayRange = 3;

