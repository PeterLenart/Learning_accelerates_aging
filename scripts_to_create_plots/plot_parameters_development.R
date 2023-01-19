# libraries - if missing they can be installed by install.packages("name of the library")

library(ggplot2)
library(reshape2)
library(dplyr)

# Set the directory path

dir <- "C:\\Users\\example\\desktop\\work\\parameters\\" 


# Get a list of all the files in the directory that match the pattern "data"

files <- list.files(dir)

files <- grep("parameters", files, value = TRUE, invert = FALSE)
files <- grep("gm", files, value = TRUE, invert = TRUE)  # In case you want to plot development of parameters from GM model - set invert to FALSE

# Iterate over the list of files and open each one
for (file in files) {
  
  file_path <- file.path(dir, file)
  # Open the file
  file <- sub("_parameters.csv", "", file) # scripts expect parameters to be saved as country_parameters.csv
  file_data <- read.csv( file_path)
  file_data$country <- file
  
  if(exists("plot_data"))
    {
    plot_data <- rbind(plot_data ,  file_data)
      }else
      {
        plot_data <- file_data
      }
      rm(file_data)
}

# graphical settings of plots

common_theme <- function() {  
  ptcolor <- 'black' # plot text color
  theme(
    plot.title=element_text(size=20, lineheight=0.8, color="black", hjust=0.5),
    axis.title.x=element_text(color=ptcolor, size =18),
    axis.title.y=element_text(color=ptcolor,size=18),
    axis.text = element_text(size = 18),
    legend.text= element_text(size=16),
    legend.title= element_text(size=18),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(size = 0.5, linetype = "solid",
                             colour = "black"),
    panel.background = element_blank())
}

# plot parameters 

ggplot(plot_data, aes(x= year, y=log10(Gmax), color=country)) + 
  geom_line(size=1)+
  labs(y="Gmax (log10)")+
  ggtitle("Gmax development")+
  common_theme()

