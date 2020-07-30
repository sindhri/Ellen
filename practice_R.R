library(readr)
library(dplyr)
data <- read_csv("individualData.csv")
data_summary <- data %>%
  group_by(genotype,plate) %>%
  summarize(mean = mean(averageActivity_day1), sd = sd(averageActivity_day1))