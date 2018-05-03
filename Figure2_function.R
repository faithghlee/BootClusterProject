#' Figure 2 (and misspecification plot generation)
#' @param  specify 1 for Non_overlap, 2 for Overlap
#' @param  pred  vector of predicted clusters from BootCluster()$out   
#' @param  n_clusters number of latent clusters 
#' @return returns Figure 2 (or any other number of clusters specified) and plots the classification 

Figure_2<-function(data_set, pred, n_clusters){
  if(missing(data_set)){
    stop("No dataset specification provided!")
  } else if (data_set == 1){
    data_set <-Non_overlap
    plot_title<-c("Trajectories from Non-overlap dataset (Scenario I)")
  } else {
    data_set <-Overlap
    plot_title<-c("Trajectories from Overlap dataset (Scenario II)")
  } 
  n_obs = 10
  data_set$est_cl<-rep(pred, each = n_obs)
  data_set$`Predicted Cluster`<- as.factor(ifelse(data_set$cluster!= data_set$est_cl, "Misspecified", data_set$est_cl))
  
  plot_2_no<-ggplot(data = data_set, aes(x = time, y = response, group = id))+ 
    geom_line() + xlab("Time") + ylab("Response Measure") + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                                               panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + ggtitle("plot_title")
  plot_2_true<-ggplot(data = data_set, aes(x = time, y = response, group = id, color = cluster))+ 
    geom_line() + xlab("Time") + ylab("Response") + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +ggtitle("Trajectories based on true clusters") +scale_color_manual(values = c("1" = "red", "2" = "orange"))
  plot_2_predict<-ggplot(data = data_set, aes(x = time, y = response, group = id, color = `Predicted Cluster`))+ 
    geom_line() + xlab("Time") + ylab("Response") + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +ggtitle("Trajectories based on predicted clusters")+ scale_color_manual(values = c("1" = "red", "2" = "orange", "Misspecified" = "black"))
  
  grid.arrange(plot_2_no, plot_2_true, plot_2_predict, layout_matrix = cbind(c(1,2), c(1,3)))
  
}