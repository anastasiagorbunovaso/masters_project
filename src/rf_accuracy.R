rf_acc <- function(res,Clusters){
  set.seed(1000)
  t0 <- Sys.time()
  
  res_new <- cbind(res,Clusters)
  res_new <- res_new %>%
    rename(Cluster = labels)
  #colnames(res_new)<- c('dim_1','dim_2','Cluster')
  
  res_new <- res_new %>%
    group_by(Cluster) %>%
    mutate(Cluster_id = cur_group_id()) %>%
    ungroup() %>%
    dplyr::select(-Cluster) %>%
    mutate(Cluster_id = as.factor(Cluster_id))
  
  intrain <- createDataPartition(y = res_new$Cluster_id, p = 0.7, list = FALSE)
  training <- res_new[intrain,]
  testing <- res_new[- intrain,]
  
  # control <- trainControl(method='repeatedcv', 
  #                         number=10, 
  #                         repeats=3)
  
  control <- trainControl(method='cv', 
                          number=5)
  
  #Number randomely variable selected is mtry
  tunegrid <- expand.grid(.mtry=c(10,20,30))
  
  rf_gridsearch <- train(Cluster_id~., 
                         data=training, 
                         method='rf', 
                         metric='accuracy', 
                         tuneGrid=tunegrid, 
                         trControl=control)
  
  pred_test_rf <- predict(rf_gridsearch,testing)
  
  cm <- confusionMatrix(as.factor(pred_test_rf) , as.factor(testing$Cluster_id))
  
  ACC <- cm$overall[1]
  
  t_end <- Sys.time()
  exeTime <- difftime(t_end, t0, units='sec')
  
  return(list(ACC,exeTime))
  
}
