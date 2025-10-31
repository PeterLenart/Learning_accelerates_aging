library(ggplot2)
library(dplyr)
library(tidyverse)

rm(list = ls())

source("path_to_the_common_theme")

dir <- "path_to_the_result_files"

files <- list.files(dir)


for (file in files ) {
  
  file_path <- file.path(dir, file)
  # Open the file
  file_data <- read.csv( file_path, header=TRUE)
  file_data$scenario <- as.factor(file_data$mean_gmax)
  
  if(exists("plot_data"))
  {
    plot_data <- rbind(plot_data ,  file_data)
  }else
  {
    plot_data<- file_data
  }
  rm(file_data)
}

# Convert factor levels to numeric, handling NA values properly
numeric_scenarios <- as.numeric(levels(plot_data$scenario))
numeric_scenarios[is.na(numeric_scenarios)] <- Inf # Replace NA with Inf for correct sorting

# Get the sorted order of numeric values
sorted_order <- order(numeric_scenarios)

# Apply this sorted order to the original levels
sorted_levels <- levels(plot_data$scenario)[sorted_order]

# Create an ordered factor with the correctly sorted levels
plot_data$scenario <- factor(plot_data$scenario, levels = sorted_levels, ordered = TRUE)


result <- plot_data %>%
  group_by(time, scenario) %>%
  summarise(
    b = mean(mean_b, na.rm = TRUE),
    sd_b = sd(mean_b, na.rm = TRUE),  # Calculate standard deviation, remove NA values
    n = n()  # Count the number of observations
  ) %>%
  mutate(
    se_b = ifelse(n > 1, sd_b / sqrt(n), NA),  # Calculate standard error if n > 1
    ci_lower = ifelse(n > 1, b - qt(0.975, df=n-1) * se_b, NA),  # Lower bound of 95% CI if n > 1
    ci_upper = ifelse(n > 1, b + qt(0.975, df=n-1) * se_b, NA)   # Upper bound of 95% CI if n > 1
  ) %>%
  ungroup()


levels_scenario <- levels(result$scenario)

# Create the ggplot
ggplot(result, aes(x = time, y = b, color = scenario)) +
  geom_line() +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = scenario), alpha = 0.2) +
  scale_fill_discrete() +
  scale_color_discrete() +
  common_theme() +  # Replace with your theme function
  labs(
    title = "Mean and Confidence Intervals Over Time by Scenario",
    x = "Time",
    y = "Mean Value",
    color = "Scenario",
    fill = "Scenario"
  )

# choose example scenarios to plot
low <- subset(result,scenario == 0.055) 
high <- subset(result,scenario == 0.175)

two <- rbind(low,high)

levels_scenario <- levels(two$scenario)
color_palette <- c("red", "blue")
# 
# Create the ggplot
ggplot(two, aes(x = time, y = b, color = scenario)) +
  geom_line() +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = scenario), alpha = 0.2) +
  scale_fill_manual(values = color_palette) +
  scale_color_manual(values = color_palette) +
  common_theme() +  # Replace with your theme function
  labs(
    title = "Mean and Confidence Intervals Over Time by Scenario",
    x = "Time",
    y = "Mean Value",
    color = "Scenario",
    fill = "Scenario"
  )
