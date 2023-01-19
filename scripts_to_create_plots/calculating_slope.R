# libraries - if missing they can be installed by install.packages("name of the library")

library(ggplot2)
library(data.table)

# Unlike the other scripts this one works for unprocessed results from sexual simulation containing only one scenario e.g., learning Lmax = 0.1 

seting =0.150 # seting of Lmax
name = "625_100timesteps015" #name of the raw txt file with the results
time= 625 #number of replictes

# set folder containing results
folder = "C:\\Users\\example\\desktop\\work\\simulations\\"

# set folder where you want to save information about slopes
folder2= "C:\\Users\\example\\desktop\\work\\simulations\\slopes\\"

# loads results
results <- read.csv(paste(folder,name,".txt",sep = ""), sep ="",header=FALSE)

learning <-results[c(TRUE, TRUE),] # Unlike the previous scripts this one works for unprocessed results from sexual simulation containing only one scenario e.g., learning Lmax = 0.1 
learning <-t(learning)

#calculates slope for each replicate
 for (i in 1:ncol(learning)) 
{
  
  data_temp <-data.frame(x = c(1:nrow(learning)), y = learning[,i])
  
  lm_temp <- lm(y ~ x,data_temp)
  slope_temp <- coef(lm_temp)[2]

  if(exists("slopes"))
  {
    slopes <- rbind(slopes , slope_temp)
  }else 
  {
    slopes <- slope_temp
  }
  rm(slope_temp)
  cat("colum",i,"learning","\n")
}


slopes_summ <-data.frame(x = mean(slopes), sd = sd(slopes), seting = seting, time = time) #mean slope, sd of mean slope

# saves the data about individual slopes

write.csv(slopes,paste(folder2,name,".csv",sep = ""), row.names = FALSE)

# saves the data in a file called slopes_summ_small - if this file with results already exists it adds new results to the previous ones

if (file.exists(paste(folder2,"slopes_summ_small.csv",sep = "")) == FALSE)
{ 
  write.csv(slopes_summ,paste(folder2,"slopes_summ_small.csv",sep = ""), row.names = FALSE)
} else{
  comparison_old <-read.csv(paste(folder2,"slopes_summ_small.csv",sep = ""))
  comparison_new <-rbind(comparison_old, slopes_summ )
  write.csv(comparison_new,paste(folder2,"slopes_summ_small.csv",sep = ""), row.names = FALSE)
}


rm(list = ls()) # empties global environment