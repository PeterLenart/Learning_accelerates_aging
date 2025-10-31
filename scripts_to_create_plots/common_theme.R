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
