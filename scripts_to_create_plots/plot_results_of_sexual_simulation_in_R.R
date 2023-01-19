# libraries - if missing they can be installed by install.packages("name of the library")

library(ggplot2)
library(data.table)

rm(list = ls()) # cleans global environment before analyzing new data

name = "example" # name of the txt document with results of sexual simulation you want to open
folder = "C:\\Users\\example\\desktop\\work\\" # folder with results of the sexual simulations

folder2= "C:\\Users\\example\\desktop\\work\\simulations\\" # folder where you want to save processed results
name2 = "example" # name of the saved processed results


sexual_results <- read.csv(paste(folder,name,".txt",sep = ""), sep ="",header=FALSE) # reads the text file

learning <-sexual_results[c(TRUE, FALSE),]  # saves odd rows as learning scenario

no_learning<-sexual_results[c(FALSE, TRUE),] # saves even rows as non-learning scenario


#function that calculates confidence intervals for the plot

CI_t <- function (x, ci = 0.95)
{
  `%>%` <- magrittr::`%>%`
  Margin_Error <- qt(ci + (1 - ci)/2, df = length(x) - 1) * sd(x)/sqrt(length(x))
  df_out <- data.frame( sample_size=length(x), Mean=mean(x), sd=sd(x),
                        Margin_Error=Margin_Error,
                        'CI lower limit'=(mean(x) - Margin_Error),
                        'CI Upper limit'=(mean(x) + Margin_Error)) %>%
    tidyr::pivot_longer(names_to = "Measurements", values_to ="values", 1:6 )
  return(df_out)
}


for (i in 1:ncol(learning)) #calculates statistic for each timestep in learning scenario
{
  learning_temp <- CI_t(c(learning[,i]))
  learning_temp <- data.frame(t(learning_temp))
  colnames(learning_temp) <- c(learning_temp[1,])
  learning_temp<-learning_temp[2,]
  learning_temp$x <- i

  if(exists("learning_f"))
  {
    learning_f <- rbind(learning_f , learning_temp)
  }else 
  {
    learning_f <- learning_temp
  }
  rm(no_learning_temp)
  cat("time point",i,"learning","\n")
}

for (i in 1:ncol(no_learning)) #calculates statistic for each timestep in non-learning scenario
{
  no_learning_temp <- CI_t(c(no_learning[,i]))
  no_learning_temp <- data.frame(t(no_learning_temp))
  colnames(no_learning_temp) <- c(no_learning_temp[1,])
  no_learning_temp<-no_learning_temp[2,]
  no_learning_temp$x <- i

  if(exists("no_learning_f"))
  {
    no_learning_f <- rbind(no_learning_f ,  no_learning_temp)
  }else
  {
    no_learning_f <- no_learning_temp
  }
  rm(no_learning_temp)
  cat("time point",i,"no_learning","\n")
}

learning
learning_f$scenario <- "learning"
no_learning_f$scenario <- "no learning"


plot_data <- rbind(learning_f,no_learning_f)
plot_data$Mean  <- as.numeric(plot_data$Mean)
plot_data$CI.Upper.limit <- as.numeric(plot_data$CI.Upper.limit)
plot_data$CI.lower.limit <- as.numeric(plot_data$CI.lower.limit)
colnames(plot_data) <- c("sample_size","mean","SD","Margin_error","CI_L","CI_U","x","scenario")

#sets grahpical theme of the plot
common_theme <- function() {  
  ptcolor <- 'black' # plot text color
  theme(
    plot.title=element_text(size=20, lineheight=0.8, color=ptcolor, hjust=0.5),
    axis.title.x=element_text(color=ptcolor, size =18),
    axis.title.y=element_text(color=ptcolor,size=18),
    axis.text = element_text(size = 18),
    legend.text= element_text(size=14),
    legend.title= element_text(size=14),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(size = 0.5, linetype = "solid",
                             colour = "black"),
    panel.background = element_blank())
}

#plots results
ggplot(plot_data,mapping = aes(x=x, y=mean, col=scenario, group = scenario)) +  
  geom_line(linetype = 3)+
  geom_ribbon(aes(ymin=CI_L, ymax=CI_U,fill=scenario),alpha =0.3)+
  scale_color_manual(values=c("blue","red"))+
  scale_fill_manual(values=c("blue","red"))+
  scale_size_manual(values=c(1,1))+ 
  ggtitle("Sexual model")+
  labs(x="time", y="aging rate (b)", group ="scenario", col="scenario")+
  # annotate(geom="text", x=50000, y=0.069, label="Lmax = 0.125",color="black",size = 6)+
  common_theme()

#saves results
write.csv(plot_data,paste(folder2,name2,".csv",sep = ""), row.names = FALSE) #saves proccesed results