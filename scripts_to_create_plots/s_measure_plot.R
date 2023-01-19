# libraries - if missing they can be installed by install.packages("name of the library")

library(lattice)
library(ggplot2)
library(dplyr)

# read files with calculated AIC values

cohort_gm <- read.csv("C:\\Users\\example\\desktop\\work\\parameters\\analysis\\smeasure_summary_gm_cohort.csv")
population_gm <- read.csv("C:\\Users\\example\\desktop\\work\\parameters\\analysis\\smeasure_summary_gm_population.csv")

cohort <- read.csv("C:\\Users\\example\\desktop\\work\\parameters\\analysis\\smeasure_summary_cohort.csv")
population <- read.csv("C:\\Users\\example\\desktop\\work\\parameters\\analysis\\smeasure_summary_population.csv")

# read names of countries
countries <- c(colnames(cohort[2:9]))

# create list of data frame names
list_data <-c("cohort_gm","population_gm","cohort","population")

# iterate through the list of data frames
for (i in 1:4)
  {
  # read in the data frame
  full <- get(list_data[i])
  # remove the original data frame
  remove(list = paste(list_data[i],sep = ""), envir = .GlobalEnv)
  #iterate through the columns of the data frame
  for (f in 2:ncol(full))
    {
    #reformat the data frame as a data frame with 3 columns: "AIC", "year", and "country"
    temp <- data.frame(smeasure=full[,f],year= full[,1], country = countries[f-1])
    # check if the data frame already exists
    if (!exists(paste(list_data[i],sep = "")))
      {
      # if not, assign the temp data frame to the original name
      assign((list_data[i]), temp, envir = .GlobalEnv)
       }
    else{
      # if yes, then bind the temp data frame to the original data frame
      assign((list_data[i]), rbind(get(list_data[i]),temp), envir = .GlobalEnv)
    }
  }
}


# distinguish types of fitted data 
cohort_gm$type <- "cohort"
population_gm$type <- "population"
cohort$type <- "cohort"
population$type <- "population"

# create two new data frames, one for GMM and one for GLA
data_gm<- rbind(cohort_gm,population_gm)
data_gm$model <- "GMM"
data <- rbind(cohort,population)
data$model <- "GLA"

colnames(data)<- colnames(data_gm) 
# combine the two data frames into a single data frame
plot_data  <- rbind(data,data_gm) 

# define a function to set various visual elements for the final plot
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


# create a jitter plot of GLA and GMM s measures
ggplot(plot_data, aes(x= model, y=smeasure, color=country)) + 
  geom_jitter(position=position_jitter(0.4), size =2,alpha=0.7)+
  labs(y="s measure, in %")+
  xlab(NULL)+
  theme(legend.position = "none")+
  ggtitle("Comparison of GLA and GMM models")+
  common_theme()

# create a jitter comparing s measures of GLA model fitted on population and cohort data
ggplot(data, aes(x= type, y=smeasure, color=country)) + 
  geom_jitter(position=position_jitter(0.4), size =2,alpha=0.7)+
  labs(y="s measure, in %")+
  xlab(NULL)+
  ggtitle("Comparison of cohort and population data")+
  common_theme()
