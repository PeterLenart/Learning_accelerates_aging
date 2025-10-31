library(ggplot2)
library(dplyr)
library(tidyverse)
library("lattice")

rm(list = ls())

source("path_tocommon_theme.R")

dir <- "path_to_results_with_different_amplitudes_of_one_d\\"


files <- list.files(dir)

# Opens and binds all the files
for (file in files  ) {
  
  file_path <- file.path(dir, file)
  # Open the file
  file_data <- read.csv( file_path, header=TRUE)
  file_data$start_d <- file_data$mean_lmax
  file_data$start_b <- file_data$mean_b[1] #average b in first timestep is the starting step
  file_data <- subset(file_data, time == 999) #we are comparing slope after 1000 timesteps
  # file_data <- file_data[,-c(2:5)]
  file_data <- na.omit(file_data)
  
  if(exists("plot_data"))
  {
    plot_data <- rbind(plot_data ,  file_data)
  }else
  {
    plot_data<- file_data
  }
  rm(file_data)
}

plot_data$start_d <- round(plot_data$start_d ,3)
plot_data$start_b <- round(plot_data$start_b,3)

plot_datab<- plot_data %>%
  group_by(improvement_strength,start_b) %>%
  summarise(z = mean(mean_b))


# Heatmap relative to no change

reference <- subset(plot_datab, improvement_strength == 0)
reference <- reference[,-1]
reference <- reference%>%
  rename(control = z)
plot_data2 <- merge(plot_datab,reference, by = "start_b")
plot_data2$value <- plot_data2$z - plot_data2$control

ggplot(plot_data2, aes(improvement_strength, start_b, fill= value)) + 
  geom_tile()+
  scale_fill_gradient2(low = "red", high = "blue", mid = "white", 
                       midpoint = 0) +
  common_theme() 


# control - no improvements
dir2 <- "path_to_control_results_with_no_improvement\\"


files <- list.files(dir2)

# Opens and binds all the files
for (file in files  ) {
  
  file_path <- file.path(dir2, file)
  # Open the file
  file_data <- read.csv( file_path, header=TRUE)
  file_data$start_b <- round(file_data$mean_b[1],3)
  file_data <- subset(file_data, time == 999)
  file_data$improvement_strength <- 0
  file_data <- na.omit(file_data)
  
  if(exists("plot_data3"))
  {
    plot_data3 <- rbind(plot_data3, file_data)
  }else
  {
    plot_data3<- file_data
  }
  rm(file_data)
}

plot_data4<- plot_data3 %>%
  group_by(start_b) %>%
  summarise(control = mean(mean_b))

# Heatmap relative to control with no survival improvements

final_plot <- merge(plot_datab,plot_data4, by = "start_b")
final_plot$value <- final_plot$z - final_plot$control

ggplot(final_plot, aes(improvement_strength, start_b, fill= value)) + 
  geom_tile()+
  scale_fill_gradient2(low = "red", high = "blue", mid = "white", 
                       midpoint = 0) +
  common_theme() 
