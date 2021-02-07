# step 1, select a raw file and output a file with the z-score based on wild type of the same date

# input data, 4 csvs, each for a different geno (del5, del44, didy, transshet)
# For each geno, it has wild type (+/+), HOM (del44/del44), HET(del44/+), 
# heamap1, HOM
# plot z-scores (HOM contrasted with wild type) for all of them on the same figure, 
# four lines: del44/del44, del5/del5, del5/del44, didy/didy
# on the following parameters: averageActivity_day/night2/3, averageWaking day/night2/3, sleepBout_day/night/2/3
# SleepLatency_day/night/2/3, sleepLength_day/night/2/3, sleep_day/night/2/3
# exclude 171106_)A, 171106_0B
# heatmap 2, compare HOM and HET in the same figure
library(readr)
library(stringr)
library(dplyr)


# read the data
fname <- file.choose()
filename <- basename(fname)
# find the parent directory as the pathname
pathname <- gsub(filename, "", fname)
# defind the name of the output file
filename_new <- paste0(substr(filename,1,nchar(filename)-4), '_zscore.csv')
output_name <- paste0(pathname, filename_new)

df <- read_csv(fname)

# remove 'scn1lab ' before geno to create simplified geno types
extra_string = str_split(df$genotype[1], " ", simplify = TRUE)[,1]

df$genotype <- gsub(paste0(extra_string, " "), "", df$genotype)
# trim the leading space
df$genotype <- str_trim(df$genotype, side="left")

# get the unique genotypes
genotypes <- sort(unique(df$genotype))
geno_WT <- genotypes[1]

# get the dates of each plate and make a new columns for date
date <- str_split(df$plate, "_", simplify = TRUE)[,1]
df$date <- date

df_filtered <- df %>%
  filter(plate!= '171106_0A' & plate!='171106_0B')

activities_keep <- c('averageActivity_day2', 'averageActivity_day3', 
                  'averageActivity_night2', 'averageActivity_night3', 
                  'averageWaking_day2', 'averageWaking_day3', 
                  'averageWaking_night2', 'averageWaking_night3', 
                  'sleepBout_day2', 'sleepBout_day3', 
                  'sleepBout_night2', 'sleepBout_night3', 
                  'sleepLatency_day2', 'sleepLatency_day3', 
                  'sleepLatency_night2', 'sleepLatency_night3', 
                  'sleepLength_day2', 'sleepLength_day3', 
                  'sleepLength_night2', 'sleepLength_night3', 
                  'sleep_day2', 'sleep_day3', 
                  'sleep_night2', 'sleep_night3')

columns_keep <- c('genotype', 'plate','fish','date',activities_keep)

df_prepared <- df_filtered %>%
  select(columns_keep)
df_prepared[mapply(is.infinite, df_prepared)] <- NA

# create the table_means and table_stds tables
# the first column of each table is the date
# the rest of the columns are the activity columns

date_unique <- unique(df_prepared$date)
mean_table <- data.frame()
std_table <- data.frame()
for(i in 1:length(date_unique)) {
  table_seed <- data.frame("date" = date_unique[i])
  WT_date_selected <- df_prepared %>%
    filter(date == date_unique[i] & genotype == geno_WT)
  
  WT_date_selected_activities <- WT_date_selected %>%
    select(activities_keep)
  
  WT_date_selected_mean <- apply(WT_date_selected_activities, 2, mean, na.rm=TRUE)
  WT_date_selected_mean_table <- cbind(table_seed, t(data.frame(WT_date_selected_mean)))
  
  WT_date_selected_std <- apply(WT_date_selected_activities, 2, sd, na.rm=TRUE)
  WT_date_selected_std_table <- cbind(table_seed, t(data.frame(WT_date_selected_std)))

  mean_table <- rbind(mean_table, WT_date_selected_mean_table)
  std_table <- rbind(std_table, WT_date_selected_std_table)
}

# for each unique date + genotype combination, create the equivalent z score based on the same date wild type

z_table_all <- data.frame()
for(i in 1:length(date_unique)){
  for(j in 1:length(genotypes)){
    geno_date_activities_selected <- df_prepared %>%
      filter(date == date_unique[i] & genotype == genotypes[j]) 
    
    nrows <- nrow(geno_date_activities_selected)
    
    mean_selected <- mean_table %>%
      filter(date == date_unique[i]) %>%
      select(activities_keep)%>% 
      slice(rep(1:n(), each = nrows))
    
    std_selected <- std_table %>%
      filter(date == date_unique[i]) %>%
      select(activities_keep) %>%
      slice(rep(1:n(), each = nrows))
    
    z_data <- (select(geno_date_activities_selected, activities_keep) - mean_selected)/std_selected
    z_table <- geno_date_activities_selected
    ncols <- ncol(z_table)
    z_table[,5:ncols] <- z_data
    z_table_all <- rbind(z_table_all, z_table)
  }
  
  write_csv(z_table_all, output_name)
}