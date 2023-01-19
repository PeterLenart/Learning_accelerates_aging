############
# before you start
############

# load or install libraries
library(ggplot2) # if missing it can be installed by install.packages("ggplot2")
library(dplyr) # if missing it can be installed by install.packages("dplyr")

# set directory where you want to save results 

folder <- "C:\\Users\\example\\desktop\\work\\" 

############
#Simulation settings
############

no_individuals = 20000 #number of individuals per scenario 
no_simulations = 20 #number of replicates
maxtime = 5000 #number of steps simulation is run

############
#Functions settings
############

a = 1.15106610e-02  #a parameter from the Gompertz-Makeham model
c = 4.26805906e-02 #constant c parameter from the Gompertz-Makeham model 
k = 3.94526937e+01 #age at which the effect of learning is half
Lmax = 0.125 #  MAXIMAL benefit of learning 
Gmax = 1.00399978e-01  # MAXIMAL benefit of growth
n = 1 #steepness of age-dependence
h = 9.23916941e-02 #growth exponent
b_starting_mean = 0.15 #starting mean of b = rate of aging
mutation_rate = 0.02  #mutation rate
mutation_SD = 0.0075 #mutation SD

#Functions

aging = function(age,b){c+a*exp((b*age))} # Gompertz-Makeham model
learning = function(age){(Lmax/(1 + exp(n*(age-k)))) - Lmax} # Function corresponding to learning
growth = function(age){ ifelse (age < 30,(Gmax/(1 + age**h) - Gmax),(Gmax/(1 + 30**h) - Gmax)) } # corresponds to growth
death = function(hr) {ifelse (runif(1) < hr, 0,1)} #function determining death

############
# simulation
############

s_count =0

#while loop to run the simulation for a given number of replicates

while (s_count < no_simulations){

#With Learning - create data frame with simulated population

population <- data.frame(ID = c(1:no_individuals), Alive = 1, b = 1, age = 0, hr = 0)
population$b  <- rnorm(n= length(population$b) ,mean= b_starting_mean ,sd=0.001)
population$age <- rnorm(n= population$age, mean = 25, sd =10)
population$age <- round(population$age,0)
population$age[population$age < 0] <- 0


#Without Learning - create data frame with simulated population

population2 <- data.frame(ID = c(1:no_individuals), Alive = 1, b = 1, age = 0, hr = 0)
population2$b  <- rnorm(n= length(population2$b) ,mean=b_starting_mean,sd=0.001)
population2$age <- rnorm(n= population2$age, mean = 25, sd =10)
population2$age <- round(population2$age,0)
population2$age[population2$age < 0] <- 0

#runs a scenario with learning

time =0

while (time <= maxtime){
  
  #hazard rate and death
  population$hr <-  aging(population$age,population$b) + learning(population$age) + growth(population$age)
  
  #boundary check for hazard rate
  
  population$hr[population$hr > 1] <- 1
  population$hr[population$hr < 0] <- 0
  
  #update the alive status of the individuals
  
  population <- population %>%  rowwise() %>% mutate(Alive =  ifelse (runif(1) < hr, 0,1) )
  
  
  death_step <- length(which(population$Alive==0)) #counts number of deaths per time step
  
  if(exists("death_final")) #saves number of deaths per time step
  {
    death_final <- rbind(death_final,  death_step)
  }else 
  {
    death_final <-  death_step
  }
  
  
  b_step <- c(summary(population$b),(s_count+0.5)) #summary data about aging rate for a given timestep and replicate
  
  if(exists("b_final")) #saves summary data about aging rate for a given timestep and replicate
  {
    b_final <- rbind(b_final,  b_step)
  }else 
  {
    b_final <-  b_step
  }
  
  rows_to_keep <- which(!population$Alive == 0) 
  population = population[rows_to_keep,]
  
  #Reproduction
  
  
if (length(population$Alive) < no_individuals) #Check if the number of individuals that are alive is less than the target population size.
    {
  population_kids <- data.frame(ID = c((length(population$Alive)+1):no_individuals), Alive = 1, b = 1, age = 0, hr = 0) #create a new data frame 'population_kids' to hold the data of the new offspring.
  
  kids = 0  
  while (kids < (no_individuals - (length(population$Alive)))){     #creates offspring in a way that keeps population size constant
    kids = kids +1
    mutation <-  ifelse (runif(1) < mutation_rate ,rnorm(n=1, mean = 0, sd = mutation_SD),0) #Apply mutation to the aging rate of the offspring.
    b_kids_t <- sample(population$b,1) + mutation
    if(exists("b_new")) #create a variable 'b_new' to hold the aging rate data of the offspring.
    {
      b_new <- rbind(b_new,  b_kids_t)
    }else 
    {
      b_new <-b_kids_t
    }
  }
  
  if (exists("b_new"))
  {
    population_kids$b <-  b_new
  }
  population_kids$b[population_kids$b < 0] <- 0 # set aging rate values less than 0 to 0
  rm(b_new)
  population <- rbind(population, population_kids) # combine the parent population and the offspring population
  rm(population_kids)
  } 
  #Preparing for the next round
  population$age <- population$age +1
  cat(time, "learning replicate",s_count+1,"\n")
  time=time+1
}


#runs a scenario without learning

time =0

while (time <= maxtime){
  
  #hazard rate and death
  population2$hr <-  aging(population2$age,population2$b) # +  growth(population2$age)
  
  #boundary check for hazard rate,
  
  population2$hr[population2$hr > 1] <- 1
  population2$hr[population2$hr < 0] <- 0
  
  #update the alive status of the individuals
  
  population2 <- population2 %>%  rowwise() %>% mutate(Alive =  ifelse (runif(1) < hr, 0,1) )
  
  
  death_step2 <- length(which(population2$Alive==0)) #counts number of deaths per time step
  
  if(exists("death_final2")) #saves number of deaths per time step
  {
    death_final2 <- rbind(death_final2,  death_step2)
  }else 
  {
    death_final2 <-  death_step2
  }
  
  
  b_step2 <- c(summary(population2$b),s_count) #summary data about aging rate for a given timestep and replicate
 
  if(exists("b_final2")) #saves summary data about aging rate for a given timestep and replicate
  {
    b_final2 <- rbind(b_final2,  b_step2)
  }else 
  {
    b_final2 <-  b_step2
  }
  
  rows_to_keep2 <- which(!population2$Alive == 0) 
  population2 = population2[rows_to_keep2,]  
  
  #Reproduction

if (length(population2$Alive) < no_individuals)  #Check if the number of individuals that are alive is less than the target population size.
  {
  population_kids2 <- data.frame(ID = c((length(population2$Alive)+1):no_individuals), Alive = 1, b = 1, age = 0, hr = 0) #create a new data frame 'population_kids2' to hold the data of the new offspring.
  
  kids2 = 0
  while (kids2 < (no_individuals - (length(population2$Alive)))){  #creates offspring in a way that keeps population size constant
    kids2 = kids2 +1
    mutation2 <-  ifelse (runif(1) < mutation_rate ,rnorm(n=1, mean = 0, sd = mutation_SD),0) #Apply mutation to the aging rate of the offspring.
    b_kids_t2 <- sample(population2$b,1) + mutation2
    if(exists("b_new2")) #create a variable 'b_new' to hold the aging rate data of the offspring.
    {
      b_new2 <- rbind(b_new2,  b_kids_t2)
    }else 
    {
      b_new2 <-b_kids_t2
    }
  }
  
  if (exists("b_new2"))
    {
    population_kids2$b <-  b_new2
      }
  population_kids2$b[population_kids2$b < 0] <- 0 # set aging rate values less than 0 to 0
  population2 <- rbind(population2, population_kids2) # combine the parent population and the offspring population
  rm(population_kids2)
} 
  
  #Preparing for the next round
  population2$age <- population2$age +1
  cat(time, "non-learning replicate",s_count+1,"\n")
  time=time+1
  rm(b_new2)
}
s_count = s_count +1
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

# creates data.frame for easier ploting

b_final<- as.data.frame(b_final)
b_final$name <- factor("learning")
b_final$x <- 1:(maxtime+1)
b_final2<- as.data.frame(b_final2)
b_final2$name <- factor("without learning")
b_final2$x <- 1:(maxtime+1)
plot_data <- rbind(b_final,b_final2)
names(plot_data)[names(plot_data) == 'Mean'] <- 'y'

# plot each simulation run 


ggplot(plot_data,mapping = aes(x, y,group = V7, col=name, linetype=name,size=name)) + 
  geom_line()+
  scale_linetype_manual(values = c(1,1))+
  scale_color_manual(values=c("blue","red"))+
  scale_size_manual(values=c(1,1))+ 
  ggtitle("b evolution")+
  labs(x="time", y="b", col="scenarios", linetype = "scenarios", size= "scenarios")+
  common_theme()

# plot mean of simulation runs with 95% Confidence intervals, labels need to be manually adjusted


ggplot(plot_data,mapping = aes(x, y=y, col=name, group = name)) +
  stat_summary(geom ="line", fun = "mean", linetype = 3)+
  stat_summary(geom ="ribbon", fun.data = "mean_cl_boot",mapping = aes(x, y=y,fill=name), alpha =0.3)+
  scale_color_manual(values=c("blue","red"),labels = c(bquote("Learning ("~L[max]*" = 0.125)"), bquote("No Learning ("~L[max]*" = 0)")))+
  scale_fill_manual(values=c("blue","red"),labels = c(bquote("Learning ("~L[max]*" = 0.125)"), bquote("No Learning ("~L[max]*" = 0)")))+
  scale_size_manual(values=c(1,1,1))+ 
  labs(x="time", y="aging rate (b)", fill="scenarios", group ="scenarios", col="scenarios")+
  ggtitle("Asexual reproduction")+
  common_theme()

# save data

name = paste("Lmax_",Lmax," ",no_individuals,".csv",sep = "")
write.csv(plot_data,paste(folder,name,".csv",sep = ""), row.names = FALSE)

# delete all data - recommended before trying new simulation settings

rm(list = ls()) # deletes everything from global environment
    