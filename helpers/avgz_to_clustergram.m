%20210325, added optional inputs
% input, averaged z score
% output, clustergram

%todo, colorbar, verify values, burst data
% 20201211, 
% separate HOM and HET into its own clustergrams

function avgz_to_clustergram(pathname, filename)
if nargin==0
    [filename,pathname] = uigetfile('*.csv','select the mean_by_geno file');
end
avgz = readtable([pathname filename]);

% choose only the rows starting with 'HOM' or 'HET'
HOM_avgz = get_target_rows_start(avgz, 'genotype', 'HOM');
HET_avgz = get_target_rows_start(avgz, 'genotype', 'HET');

% REMOVE HOM_DMSO and HET_DMSO
HOM_avgz = HOM_avgz(strcmp(HOM_avgz.genotype, 'HOM + DMSO') ~=1,:);
HET_avgz = HET_avgz(strcmp(HET_avgz.genotype, 'HET + DMSO') ~=1,:);

fig_dir = [pathname '/fig/'];
if exist(fig_dir, 'dir') ~=7
    mkdir(fig_dir);
end

make_clustergram(HOM_avgz, fig_dir, 'HOM');
make_clustergram(HET_avgz, fig_dir, 'HET');

end

function make_clustergram(avgz, pathname, filename)

%3, removed 'group_count' from the parameters
column_start = 3;
data = avgz{:,column_start:end};
genotypes = avgz{:,1};
parameters = avgz.Properties.VariableNames(column_start:end); 
%replaced underscore with space in the row and column names
genotypes = strrep(genotypes,'_', ' ');
parameters = strrep(parameters,'_', ' ');

load my_colormap;
cgz = clustergram(data, 'RowLabels', genotypes,...
                             'ColumnLabels', parameters,...
                             'RowPdist', 'correlation',...
                             'ColumnPdist', 'correlation',...
                             'ImputeFun', @knnimpute);
cgz.Colormap = my_colormap;           
cgz.DisplayRange = 3;

% saveas did not work
% saveas(cgz, [pathname 'clustergram_' filename], 'fig');
end