library(readr)
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stats)
library(kableExtra)
library(plyr)

prepare_raw <- function(df){
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
  
  columns_keep <- c('genotype', 'plate','fish',activities_keep)
  df_filtered <- df %>%
    filter(plate!= '171106_0A' & plate!='171106_0B')
  df_prepared <- df_filtered %>%
    select(columns_keep)
  df_prepared[mapply(is.infinite, df_prepared)] <- NA
  # if the name of the activites change
  #activity_colnames <- colnames(df_prepared)[4:length(colnames(df_prepared))]
  #activities <- unique(str_split(activity_colnames,'_',simplify = TRUE)[,1])
  
#  activities <- c('sleep','sleepBout', 'sleepLength', 'sleepLatency', 'averageActivity', 'averageWaking')
#  for(activity in activities){
#    for(day_or_night in c('night', 'day')){
#      vname <- paste0(activity, '_', day_or_night, '23')  
#      vname1 <- paste0(activity, '_', day_or_night, '2')  
#      vname2 <- paste0(activity, '_', day_or_night, '3')

      df_prepared <- df_prepared %>%
        mutate(sleep_night23=sleep_night2+sleep_night3, sleepBout_night23=sleepBout_night2+sleepBout_night3, 
               sleepLength_night23=sleepLength_night2+sleepLength_night3, sleepLatency_night23=sleepLatency_night2+sleepLatency_night3, 
               averageActivity_night23=averageActivity_night2+averageActivity_night3, averageWaking_night23=averageWaking_night2+averageWaking_night3, 
              sleep_day23=sleep_day2+sleep_day3, sleepBout_day23=sleepBout_day2+sleepBout_day3, 
              sleepLength_day23=sleepLength_day2+sleepLength_day3, sleepLatency_day23=sleepLatency_day2+sleepLatency_day3, 
              averageActivity_day23=averageActivity_day2+averageActivity_day3, averageWaking_day23=averageWaking_day2+averageWaking_day3)
#    }
#  }
  return(df_prepared)
}

#expr <- paste0('df_prepared <- df_prepared %>% mutate(', vname, '= (', vname1, '+', vname2, ')/2)')
#eval(expr)

calculate_values <- function(df){
  col_start <- 4 # start statistics from column 4 for raw data, 5 for zscore
  col_end <- ncol(df)
  p_all <- vector()
  for(i in c(col_start:col_end)){
    # Compute the analysis of variance
    res.aov <- aov(eval(as.symbol(colnames(df)[i])) ~ eval(as.symbol(colnames(df)[1])), data = df)
    data_p <- summary(res.aov)[[1]][["Pr(>F)"]][1]
    p_all <- append(p_all, data_p)
  }
  fileinfos <- list()
  
  p_table <- data.frame(t(p_all))
  data_col_names <- colnames(df[col_start:col_end])
  colnames(p_table) <- data_col_names
  
  fileinfos$p <- p_table
  fileinfos$mean <- data.frame(df_prepared %>%
    group_by(genotype) %>%
    summarise_at(vars(-plate,-fish), list(~ mean(., na.rm=TRUE))))
  fileinfos$std <- data.frame(df_prepared %>%
    group_by(genotype) %>%
    summarise_at(vars(-plate,-fish), list(~ sd(., na.rm=TRUE))))

  return(fileinfos)
}

get_filenames <- function(){
  # select the files
  fname <- file.choose()
  filename <- basename(fname)
  
  # find the parent directory as the pathname
  pathname <- gsub(filename, "", fname)
  
  filelist <- list.files(pathname)
  filenames <- list()
  
  for (filename in filelist){
    if(nchar(filename) > 14){
      # select z score files
      #if(substr(filename, nchar(filename)-9, nchar(filename)-4)=="zscore"){
      filename_part1 <- substr(filename, nchar(filename)-7, nchar(filename)-4)
      filename_part2 <- substr(filename, nchar(filename)-3, nchar(filename))
      if((filename_part1 !="core") & (filename_part1 != "_ave") & (filename_part2 == '.csv')){
        filenames <- append(filenames, filename)
      }
    }
  }
  fileinfo <- list()
  fileinfo$filenames <- filenames
  fileinfo$pathname <- pathname
  return(fileinfo)
}

process_one_file <- function(pathname, filename){
  df <- read_csv(paste0(pathname, filename))
  
  # remove 'scn1lab ' before geno to create simplified geno types
  extra_string = str_split(df$genotype[1], " ", simplify = TRUE)[,1]
  df$genotype <- gsub(paste0(extra_string, " "), "", df$genotype)
  # trim the leading space
  df$genotype <- str_trim(df$genotype, side="left")
  # get the unique genotypes
  genotypes <- sort(unique(df$genotype))
  geno_WT <- genotypes[1]
  if(length(genotypes)==4){
    geno_HET1 <- genotypes[2]
    geno_HET2 <- genotypes[3]
    geno_HOM <- genotypes[4]
    geno_name <- geno_HOM
  }else{
    geno_HET <- genotypes[2]
    geno_HOM <- genotypes[3]
    geno_name <- str_split(geno_HET, '/', simplify = TRUE)[,1]
  }
  
  df_prepared <- prepare_raw(df)
  results <- calculate_values(df_prepared)
  
  write_csv(results$p, paste0(pathname, 'results_', geno_name, '_p.csv'))
  write_csv(results$mean, paste0(pathname, 'results_', geno_name, '_mean.csv'))
  write_csv(results$std, paste0(pathname, 'results_', geno_name, '_std.csv'))
}

fileinfo <- get_filenames()
for (filename in fileinfo$filenames) {
  process_one_file(fileinfo$pathname, filename)
}


#genotypes <- sort(unique(df$genotype))
#geno_WT <- genotypes[1]
#geno_HET <- genotypes[2]
#if(str_detect(filename, 'transhet') == TRUE){
#  geno_HOM1 <- genotypes[3]
#  geno_HOM2 <- genotypes[4]
#  } else {
#  geno_HOM <- genotypes[3]
#  }


# for later, a formmated table
#table1_colnames <- c('sleep', 'sleepBout', 'sleepLength')
#table1_colnames_friendly <- c('Sleep Length (min)', 'Number of Sleep Bouts', 'Sleep Bout Length (min)')
#table_colnames <- c('night2', 'night2', 'night23', 'day2', 'day3', 'day23')
#table_colnames_friendly <- c('Night 2', 'Night 3', 'Night 23', 'Day 2', 'Day 3', 'Day 23')

#table1_data <- data.frame(row.names=table1_colnames_friendly)
#for(i in c(1:length(table_colnames))){
#  vname <- paste0(table1_colnames[1],'_', table_colnames[i])
#  for(j in c(4:length(colnames(df_prepared)))){
#    if(colnames[j]==vname){
#      index <- j
#      break
#    }
#    index <- index- 3
#  }
#  cell_mean <- fileinfos$mean[index]
#  cell_std <- fileinfos$std[index]
#}

# in case you want to plot some variables
#ggplot(df, aes_string(x=colnames(df)[1], y=colnames(df)[5], color=colnames(df)[1])) + 
#  geom_boxplot()

