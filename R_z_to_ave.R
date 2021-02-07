library(readr)
library(stringr)
library(dplyr)


# select the zscore files
fname <- file.choose()
filename <- basename(fname)

# find the parent directory as the pathname
pathname <- gsub(filename, "", fname)
# defind the name of the output file
if (substr(filename, nchar(filename) - 9, nchar(filename) - 4) != 'zscore'){
  print('Please choose an zscore file.Abort')
}
filename_new <- paste0(substr(filename,1,nchar(filename)-4), '_ave.csv')
output_name <- paste0(pathname, filename_new)

df <- read_csv(fname)
genotypes <- sort(unique(df$genotype))
z_ave_table_all <- data.frame()
for(genotype1 in genotypes){
  geno_selected <- df %>%
    filter(genotype == genotype1) 
    ncols <- ncol(geno_selected)
    z_ave <- apply(geno_selected[,5:ncols], 2, mean, na.rm=TRUE)
    table_seed <- data.frame(genotype=genotype1)
    z_ave_table <- cbind(table_seed, t(z_ave))
    z_ave_table_all <- rbind(z_ave_table_all, z_ave_table)
}

write_csv(z_ave_table_all, output_name)