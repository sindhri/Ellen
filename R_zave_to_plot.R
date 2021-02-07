# heatmap1, HOM (del5/del5)
# heatmap2, compare HOM and HET in the same figure

library(readr)
library(ggplot2)
library(stringr)
library(dplyr)
library(tidyr)

prepare_plot_data <- function(data){
  ncols <- ncol(data)
  data <- data %>%
    select(1:(ncols-1))
  data_t <- t(data[,2:ncol(data)])
  rownames(data_t) <- NULL
  colnames(data_t) <- gsub("/", "_",data$genotype)
  
  activity_col <- colnames(data)[2:ncol(data)]
  data_t <- cbind(data.frame(activity=activity_col), data_t)
  data_t$activity <- factor(data_t$activity)
  cols_to_plot <- colnames(data_t)[2:ncol(data_t)]
  
  data_t_longer <- data_t %>%
    pivot_longer(cols = colnames(data_t)[2:ncol(data_t)],names_to = "geno", values_to = "zscore")
  return(data_t_longer)
}

# select the zscore files
fname <- file.choose()
filename <- basename(fname)

# find the parent directory as the pathname
pathname <- gsub(filename, "", fname)

filelist <- list.files(pathname)
zscore_ave_files <- list()

for (filename in filelist){
  if(nchar(filename) > 14){
    if(substr(filename, nchar(filename)-13, nchar(filename)-4)=="zscore_ave"){
      zscore_ave_files <- append(zscore_ave_files, filename)
    }
  }
}

data <- data.frame()
for (filename in zscore_ave_files){
  zscore_ave_table <- read_csv(paste0(pathname, filename)) %>%
    filter(genotype != '+/+')
  if(str_detect(filename, 'transhet')){
    zscore_ave_table$genotype <- paste0(zscore_ave_table$genotype, 't')
  }
  data <- rbind(data, zscore_ave_table)
}

genotype_part2=str_split(data$genotype,'/',simplify=TRUE)

HOM_data <- data %>%
  mutate(geno2 = genotype_part2[,2]) %>%
  filter(geno2 != '+' & geno2 != '+t') 

HOM_plot <- prepare_plot_data(HOM_data)
  
ggplot(data = HOM_plot, aes(x=activity, y = zscore, group=geno,color = geno)) +
  geom_line()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ggsave('HOM.png', height = 7 , width = 20)

HOM_HET_plot <- prepare_plot_data(data)

ggplot(data = HOM_HET_plot, aes(x=activity, y = zscore, group=geno,color = geno)) +
  geom_line()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ggsave('HOM_HET.png', height = 7 , width = 20)