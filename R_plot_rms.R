library(readr)
library(ggplot2)

# Define the pathname for the data and path for the figures
pathname <- 'drug_discovery_1210/'
pathname_fig <- paste0(pathname, 'fig/')

# Make a directory if the subfolder 'fig' does not exist
if (!dir.exists(pathname_fig)) {
  mkdir(pathname_fig)
}

# read the data
data <- read_csv(paste(pathname, 'individualData_121020_zscore_rms_HOM.csv', sep=""))

# generate a rms plot for each column that starts with 'rms'
col_names <- colnames(data)
for(col_name in col_names) {
  if (substr(col_name, 1, 3) == 'rms') {
    p<-ggplot(data=data, aes(x=reorder(genotype, -!!sym(col_name)), y=!!sym(col_name))) +
      geom_bar(stat="identity")
    
    p + coord_flip() + scale_y_continuous(position = 'right')
    
    ggsave(paste0(pathname_fig, col_name, ".png"), width = 10, height = 8, device='png')
  }
}