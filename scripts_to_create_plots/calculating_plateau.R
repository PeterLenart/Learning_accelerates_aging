# libraries - if missing they can be installed by install.packages("name of the library")

library(ggplot2)
library(dplyr)

folder2= "" #folder where the processed results are saved 

name = ""  #name of the folder with results
plot_data<- as.data.frame(read.csv(paste(folder2,name,".csv",sep = "")))

# Calculate the difference between the average of the first 1250 values and the average of the last 1250 values in each interval
sum_data <- plot_data %>%
  group_by(interval = floor(x / 2500)) %>%
  summarize(mean_start = mean(mean[seq(1, 1250)]), mean_end = mean(mean[seq(1250, 2500)]))

# Identify timepoints at which the mean value at the start of the interval is different from the mean value at the end of the interval by 0.5% of starting value
stable_periods <- sum_data %>%
  filter(abs(mean_start - mean_end) < 0.000350) %>%
  select(interval)

# Print the timepoints at which the mean was stable
print(stable_periods)


