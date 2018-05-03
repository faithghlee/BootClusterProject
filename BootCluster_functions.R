#' Generate a Bootstrap sample
#' @param N_B the sample size of the generated Bootstrap sample
#' @param n_obs number of observations per individual, default passed in as 10 
#' @param data_set the original sample where the Bootstrap sample is generated from
#' @return a data set of the Bootstrap sample
genBootDataSet <- function(N_B, get_id, data_set, n_obs){
  bs_set<-c()
  #Create bootstrap dataset 
  for( j in 1:N_B){
    bs_set<-rbind(bs_set, subset(data_set, data_set$id == get_id[j]))
  }
  #Bootstrap id, so each individual re-sampled twice are seen as "different" individuals in LME fit. 
  bs_set$bs_id <- rep(1:N_B, each = n_obs)
  
  return(bs_set)
}

#' Generate comparing-basis data set to be used in clustering procedure
#' @param lme.bs fitted result from lme model
#' @param N_B sample size, which is the row number of the comparing-basis data set
#' @return a data set describing the comparing basis, i.e., the data set will be partitioned via coclustering
genClustDataSet <- function(lme.bs,bs_set, N_B){
  
  # assign fixed-effect coeff.
  beta_int<-fixed.effects(lme.bs)[1]
  beta_x1<-fixed.effects(lme.bs)[2]
  beta_x2<-fixed.effects(lme.bs)[3]
  beta_t<-fixed.effects(lme.bs)[4]
  
  Sigma_bs <-getVarCov(lme.bs) #Returns the covariance structure of the random effects 
  # generate random-effect coeff. draws (once) for each individual
  bs_theta_j <- mvrnorm(n = N_B, mu = c(mean(random.effects(lme.bs)[,1]), mean(random.effects(lme.bs)[,2])), Sigma = Sigma_bs)
  b_bs_int <- bs_theta_j[,1]
  b_bs_t <- bs_theta_j[,2]
  
  #Get the X_1 and X_2 values for each unique bs_id individual in bootstrap dataset 
  bs_x1 <- bs_set[!duplicated(bs_set$bs_id),'X1'] 
  bs_x2 <- bs_set[!duplicated(bs_set$bs_id),'X2']
  
  #Calculate the Time-dependent effect and #Time independent effect 
  # for each individual in the dataset 
  
  #Time independent effect 
  #B1X1 + B2X2 + alpha_0 + alpha_i 
  time_indep <- beta_int + b_bs_int + beta_x1*bs_x1 + beta_x2*bs_x2 
  #Time dependent effect 
  #(beta_T +b_i) 
  time_dep <- beta_t + b_bs_t
  
  #Make a dataset on the bootstrap individual based on the time-independent and 
  #time- dependent effect, the dataset will be used in coclustering function 
  cluster_bs_data <- data.frame( id = unique(bs_set$bs_id), independent = time_indep, dependent = time_dep )
  
  return(cluster_bs_data)
  
}

#' Reorganize the clustering result obtained from coclustering procedure by cluster mean (with increasing order) and generate a cluster-indicator matrix 
#' @param coc.fit fitted results returned by a cocluster function
#' @param get_id individuals which are selected into the Boostrap sample
#' @param N the sample size of the original sample
#' @param n_cluster number of clusters used to subset the sample
#' @return a matrix of cluster indicators, each row vector contains 0 (not in the specific cluster) and 1, representing 
#' an individual's cluster estimates. 
orgnCluster <- function(coc.fit, get_id, N, cluster_bs_data, n_cluster){
  
  rc <- coc.fit@rowclass
  cc <- coc.fit@colclass
  coc.fit@classmean
  
  # Reorganize the clustering result obtained from coclustering procedure by cluster mean (with increasing order)
  rank_clmean <- rank(coc.fit@classmean[,1])
  correct_rc <- unlist(sapply(rc+1, FUN = function(i){rank_clmean[i]}))
  newcorrect_rc <- numeric(N)
  newcorrect_rc[cluster_bs_data$id] <- correct_rc
  final_rc <- newcorrect_rc
  
  # this very part is for Bootstrap sampling specifically, since individual i can be selected into the B sample more than once.
  # Thus, it is possible that individual i can be classfied into different cluster in the same clustering procedure. We would like
  # to record all cluster estimates.
  test<-data.frame(get_id, final_rc = as.factor(final_rc))
  test_2 <-test %>% group_by(get_id, final_rc) %>% summarize(count = n())
  test<-data.frame(get_id, final_rc = final_rc)
  test_2 <-test %>% group_by(get_id, final_rc) %>% summarize(count = n())
  premat <- matrix(0, ncol=n_cluster, nrow=N)
  premat[as.matrix(unique(test[order(test$get_id),]))] <- test_2$count
  
  return(premat)
}  

#' Generate a final cluster estimator for a sample of individuals using Bootstrap method
#' @param B the total number of Boostrap samples
#' @param data_set specify 1 for Non_overlap, 2 for Overlap 
#' @param N_B the sample size of each Boostrap sample
#' @param n_cluster number of clusters used to subset the sample
#' @return a vector of cluster estimator indicating each individual's cluster
BootCluster <- function(B=1e3, N, N_B=N, data_set = 1, n_cluster=2){
  if(missing(data_set)){
    stop("No dataset specification provided!")
  } else if (data_set == 1){
    data_set <-Non_overlap
  } else {
    data_set <-Overlap
  } 
  n_obs = 10
  # to accumulate the cluster estimates in each Bootstrap step for each individual
  c_mat <- matrix(0, nrow = N, ncol = n_cluster)
  
  b = 0
  iter = 0
  while(b <= B){
    
    iter <- iter + 1
    get_id = sample(1:N, size = N_B, replace = TRUE)
    
    #LME FIT
    argList <- list("N_B"=N_B,
                    "get_id"=get_id,
                    "n_obs" = n_obs,
                    "data_set"=data_set)
    bs_set <- do.call(genBootDataSet, argList)
    ctrl <- lmeControl(opt='optim')
    lme.bs <- try(lme(fixed = response~ X1 + X2 + time, random = ~ 1 + time|bs_id, data = bs_set, control = ctrl), TRUE)
    if("try-error" %in% class(lme.bs)) next
    
    #CLUSTERING 
    argList <- list("N_B"=N_B,
                    "bs_set"=bs_set,
                    "lme.bs"=lme.bs)
    cluster_bs_data <- do.call(genClustDataSet, argList)
    coc.fit <- blockcluster::coclusterContinuous(as.matrix(cluster_bs_data[,c(2,3)]), nbcocluster = c(n_cluster,1), model = 'pik_rhol_sigma2kl')
    
    if(coc.fit@message != "Co-Clustering successfully terminated!") next
    
    premat <- orgnCluster(coc.fit, get_id, N, cluster_bs_data, n_cluster)
    
    c_mat <- c_mat + premat
    
    print(sum(premat))
    
    b <- b + 1
    cat(b,"\n")
  }
  
  # decide the cluster for each individual by mode
  out <- apply(c_mat, 1, which.max)
  result <- list("out"=out, "c_mat"=c_mat, "iter"=iter)
  return(result)
  
}