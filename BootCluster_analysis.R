# cluster estimator from Boostrap method
B = 1000#Number of bootstrap samples 
N_B = N #Sample Size of each Bootstraped sample
source("R/BootCluster_functions.R")
results <- vector("list", length = 4)
for( l in 2:5){
  
  n_cluster = l # number of clusters
  fitr <- BootCluster(B, N, N_B, n_cluster=2)
  B_clust <- fitr$out
  conv <- fitr$iter/B
  c_mat <- fitr$c_mat
  
  #---------------------------------- CLUSTERING EVALUATION -------------------------------------#
  # true cluster
  true <- unique(data_set[,c(1,6)])[,2]
  # estimated cluster
  out <- B_clust #apply(c_mat, 1, which.max) #apply(apply(final_rc,2, table)[-1,],2,which.max)
  testp <- data.frame(true=true, out=as.factor(out))
  testp_2 <- testp %>% group_by(true, out) %>% summarize(count = n())
  
  crsp_cluster <- matrix(0, nrow=n_cluster, ncol=n_cluster)
  crsp_cluster[cbind(testp_2$true, testp_2$out)] <- testp_2$count
  print(crsp_cluster)

  #---------------------------------- MODEL COMPARISON ------------------------------------#
  
  data_set$B_clust <- rep(B_clust, each=n_obs)
  
  ctrl <- lmeControl(opt='optim')
  lme.ori <- lme(fixed = response~ X1 + X2 + time, random = ~ 1 + time|id, data = data_set
                 , control = ctrl, method = "ML")
  lme.tc <-  lme(fixed = response~ (X1 + X2 + time)*cluster, random = ~ 1 + time|id, 
                 data = data_set, control = ctrl, method = "ML")
  lme.bc <- lme(fixed = response~ (X1 + X2 + time)*B_clust, random = ~ 1 + time|id, 
                data = data_set, control = ctrl, method = "ML")
  grid.arrange(
    plot(lme.ori,type=c("p","smooth")),
    plot(lme.tc,type=c("p","smooth")),
    plot(lme.bc,type=c("p","smooth")))
  
  # LRT: model selection/comparison
  liket <- anova(lme.ori,lme.bc, lme.tc)
  
  results[[l]] <- list("B_clust"=B_clust, "crsp_cluster"=crsp_cluster, "liket"=liket, "convn"=conv, "c_mat"=c_mat)
}
#----------------------------------------------------------------------------------------#