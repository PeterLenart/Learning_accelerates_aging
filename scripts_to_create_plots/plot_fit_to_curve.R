# libraries - if missing they can be installed by install.packages("name of the library")

library(ggplot2)

# Select year and country you want to plot

selected_year = 1951
country = "France"

# select directory with GLA fitted parameters in format country_parameters.csv

dir <- "C:\\Users\\example\\desktop\\work\\parameters\\"

# select directory with GLA fitted parameters in format country_parameters_gm.csv

dir2 <- "C:\\Users\\example\\desktop\\work\\parameters\\"

# select directory with mortality data (from human mortality database https://www.mortality.org/)  the file should contain per year mortality and have name equal to country name e.g., franceMx_1x1.csv you need to import the data into excell and save as csv for it to work well, (the file should have 5 collums, first is year and last is totoal mortality)

dir3 <- "C:\\Users\\example\\desktop\\work\\mortality_data\\"


# read GLA parameters, GM parameters and  mortality data

parameters <- read.csv(paste(dir,country,"_parameters.csv",sep=""))
gompertz <- read.csv(paste(dir2,country,"_parameters_gm.csv",sep=""))
real <- read.csv(paste(dir3,tolower(country),"Mx_1x1.csv",sep=""))

colnames(real) <- c("year","age","f_mx","m_mx","mx")

# subset data frames to year you are interested in
parameters <- subset(parameters, year == selected_year)
gompertz <- subset(gompertz, year ==selected_year)
real <- subset(real, year ==selected_year)

# Assign selected GLA parameters 

a = parameters$a
b = parameters$b
c = parameters$c
Lmax = parameters$Lmax
Gmax = parameters$Gmax
n = parameters$n
k = parameters$k_learning
h = parameters$growth_speed

# GLA function (alternative equivalent form)

GLA = function(age){ifelse(age < 30,((c+a*exp((b*age)))+(((Lmax/(1 + exp(n*(age-k)))) - Lmax))+ (Gmax/(1 + age**h) - Gmax)),((c+a*exp((b*age)))+(((Lmax/(1 + exp(n*(age-k)))) - Lmax))+(Gmax/(1 + 30**h) - Gmax)))}

# Create GLA prediction

plot_data <-data.frame(age=c(0:65),mx=NA)
plot_data$mx <- GLA(plot_data$age)
plot_data$model <- "GLA"

# Assign selected  parameters of GMM model
a <-gompertz$a
b <-gompertz$b
c <-gompertz$c

# Create GMM prediction

GMM= function(age){(c+a*exp((b*age)))}
gomp_data <-data.frame(age=c(0:65),mx=NA)
gomp_data$mx <- GMM(gomp_data$age)
gomp_data$model <- "GMM"

# bind GLA and GMM prediction into one data frame

plot_data <- rbind(plot_data,gomp_data)


# Select relevant data  - select ages of interest (those on which models were fitted) and total mortality
real <- real[0:65,c(1:2,5)]
real$model <- "data"
real$mx <- as.numeric(real$mx)
real$age <- as.numeric(real$age)

# Graphical setting

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

# plot

ggplot(plot_data, aes(x= age, y=log10(mx), color = model)) + 
  geom_line(size =1.5, linetype=1)+
  labs(y="mortality rate by year (log10)")+
  geom_point(data = real,aes(x= age,y=log10(mx),group = model),size=1.5, color= "black")+
  ggtitle(paste(country," ", selected_year, sep=""))+
  common_theme()

