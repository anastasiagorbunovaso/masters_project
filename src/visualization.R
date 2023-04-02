viz_dim_red <- function(data,Clusters,alg_name,file_to_work) {
  
  plot <- ggplot(data = data, aes(x = data[,1], y = data[,2],color=Clusters))+geom_point(size = 2)+
    labs(title = alg_name, 
         x = "dimension 1", 
         y = "dimension 2")+
    theme(plot.title=element_text(hjust = 0.5))+
    theme(axis.text.y = element_text(size = 30, face = "italic")) +
    theme(axis.text.x = element_text(size = 30, face = "italic"))+ 
    theme(plot.title=element_text(hjust = 0.5))+
    theme(legend.title = element_text(size = 25)) +
    theme(legend.text = element_text(size = 20)) +
    guides(color = guide_legend(override.aes = list(size = 15)))+
    theme(plot.title = element_text(size = 45),
          axis.title.x = element_text(size = 40),
          axis.title.y = element_text(size = 40))
  
  png(paste0(output_file_dir,file_to_work,'/',alg_name,".png"), width = 1000, height = 600, res = 96)
  print(plot)
  dev.off()
  
}
