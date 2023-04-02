metrics <- function(res,Clusters){
  
  res_new <- cbind(res,Clusters)
  
  
  colnames(res_new)<- c('dim_1','dim_2','Cluster')
  
  res_new <- res_new %>%
    group_by(Cluster) %>%
    mutate(Cluster_id = cur_group_id()) %>%
    ungroup()
  
  mean_dist <- 0
  for(cl in unique(res_new$Cluster_id)){
    
    res_temp <- res_new %>%
      filter(Cluster_id == cl) %>%
      dplyr::select(dim_1,dim_2)
    
    res_temp <- as.matrix(rdist(res_temp))
    
    mean_dist <- mean_dist + abs(sum(res_temp[upper.tri(res_temp)])/ncol(res_temp))
    
  }
  mean_dist <- mean_dist/n_distinct(res_new$Cluster_id)
  
  
  
  mean_center <- res_new %>%
    group_by(Cluster_id) %>%
    summarise(CENTER_1 = median(dim_1),
              CENTER_2 = median(dim_2)) %>%
    ungroup() %>%
    dplyr::select(-Cluster_id)
  
  mean_center <- as.matrix(rdist(mean_center))
  
  mean_center <- abs(sum(mean_center[upper.tri(mean_center)])/ncol(mean_center))
  
  Q1 = mean_dist/mean_center
  
  
  
  ######
  mean_dist <- 0
  for(cl in unique(res_new$Cluster_id)){
    
    res_temp <- res_new %>%
      filter(Cluster_id == cl) %>%
      dplyr::select(dim_1,dim_2)
    
    res_temp <- as.matrix(rdist(res_temp)^2)
    
    mean_dist <- mean_dist + abs(sum(res_temp[upper.tri(res_temp)])/ncol(res_temp))
    
  }
  
  
  mean_center <- res_new %>%
    group_by(Cluster_id) %>%
    summarise(CENTER_1 = median(dim_1),
              CENTER_2 = median(dim_2)) %>%
    ungroup() %>%
    dplyr::select(-Cluster_id)
  
  mean_center <- as.matrix(rdist(mean_center)^2)
  
  mean_center <- sum(mean_center[upper.tri(mean_center)])
  
  Q2 = mean_dist/mean_center
  
  return(list(Q1,Q2))
  
}
