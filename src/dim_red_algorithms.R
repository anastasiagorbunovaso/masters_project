#PCA
PCA <- function(data)
{
  t0 <- Sys.time()
  mean_x = apply(data, MARGIN = 2, FUN = mean)
  sd_x= apply(data, MARGIN = 2, FUN = sd)
  norm_x = base::matrix(0, nrow = nrow(data), ncol = ncol(data))
  
  for (i in 1:nrow(data))
  {
    for (j in 1:ncol(data))
    {
      norm_x[i,j] = (data[i,j] - mean_x[j])/sd_x[j];
    }
  }
  norm_x_t = t(norm_x)
  R = (norm_x_t%*%norm_x)/(nrow(data)-1)

  mean_norm_x = apply(norm_x, MARGIN = 2, FUN = mean)
  sd_norm_x= apply(norm_x, MARGIN = 2, FUN = sd)
  Eig_R = eigen(R)
  Delta = Eig_R$values
  Delta = matrix(Delta,nrow=ncol(data),ncol=ncol(data))
  for (i in 1:ncol(data))
  {
    for (j in 1:ncol(data))
    {
      if(i!=j)
        Delta[i,j]=0
    }
  }
  
  A = Eig_R$vectors
  Z = norm_x %*% A
  Z <- as.data.frame(scale(Z))
  Z <- Z[,1:2]
  
  t_end <- Sys.time()
  exeTime <- difftime(t_end, t0, units='sec')
  return(list(Z,exeTime))
}

#PCA with dopar
PCA_dopar <- function(data) {
  t0 <- Sys.time()
  cl <- makeCluster(detectCores())
  registerDoParallel(cl)
  
  mean_x <- foreach(j = 1:ncol(data), .combine = 'c') %dopar% {
    mean(data[, j])
  }
  sd_x <- foreach(j = 1:ncol(data), .combine = 'c') %dopar% {
    sd(data[, j])
  }
  
  norm_x = matrix(0, nrow = nrow(data), ncol = ncol(data))
  
  norm_x <- as.matrix(foreach(i = 1:nrow(data), .combine = rbind) %dopar% {
    (data[i,] - mean_x) / sd_x
  })
  
  norm_x_t = t(norm_x)
  
  R = (norm_x_t%*%norm_x)/(nrow(data)-1)
  
  mean_norm_x <- foreach(j = 1:ncol(norm_x), .combine = 'c') %dopar% {
    mean(norm_x[, j])
  }
  
  sd_norm_x <- foreach(j = 1:ncol(norm_x), .combine = 'c') %dopar% {
    sd(norm_x[, j])
  }
  
  Eig_R = eigen(R)
  Delta = Eig_R$values
  Delta = matrix(Delta,nrow=ncol(data),ncol=ncol(data))
  
  Delta <- foreach(i = 1:nrow(data), .combine = 'cbind') %dopar% {
    Delta_row <- rep(0, ncol(data))
    for (j in 1:ncol(data)) {
      if (i != j) {
        Delta_row[j] <- 0
      } else {
        Delta_row[j] <- Delta[i, j]
      }
    }
    Delta_row
  }
  
  
  A = Eig_R$vectors
  
  Z <- foreach(i = 1:nrow(norm_x), .combine = rbind) %dopar% {
    norm_x[i,] %*% A
  }
  
  Z <- as.data.frame(scale(Z))
  Z <- Z[,1:2]
  stopCluster(cl)
  
  t_end <- Sys.time()
  exeTime <- difftime(t_end, t0, units='sec')
  return(list(Z,exeTime))
}


#PCA with snow
PCA_snow <- function(data, labels) {
  t0 <- Sys.time()
  data <- cbind(data,labels)
  
  # Set up cluster
  #cl <- makeCluster(num_cores)
  cl <- makeCluster(detectCores())
  clusterEvalQ(cl, library(Matrix))
  
  # Split data into chunks
  chunks <- split(data, 1)
  
  # Distribute chunks to workers
  clusterExport(cl, "apply")
  
  Z <- parLapply(cl, chunks, function(chunk) {
    labels <- chunk[,ncol(chunk)]
    chunk <- chunk[,-ncol(chunk)]
    # Perform PCA on chunk
    mean_x = apply(chunk, MARGIN = 2, FUN = mean)
    sd_x= apply(chunk, MARGIN = 2, FUN = sd)
    norm_x = matrix(0, nrow = nrow(chunk), ncol = ncol(chunk))
    
    for (i in 1:nrow(chunk))
    {
      for (j in 1:ncol(chunk))
      {
        norm_x[i,j] = (chunk[i,j] - mean_x[j])/sd_x[j];
      }
    }
    norm_x_t = t(norm_x)
    R = (norm_x_t%*%norm_x)/(nrow(chunk)-1)
    Eig_R = eigen(R)
    A = Eig_R[2]
    A = as.matrix(A[[1]])
    Z = norm_x %*% A
    Z <- as.data.frame(scale(Z))
    Z <- cbind(Z,labels)
    Z <- Z[,c(1:2,ncol(Z))]
  })
  
  # Stop cluster
  stopCluster(cl)
  
  # Combine results
  full_projected_data <- do.call(rbind, Z)
  
  # Return projected data
  full_projected_data
  
  t_end <- Sys.time()
  exeTime <- difftime(t_end, t0, units='sec')
  return(list(full_projected_data,exeTime))
  
}

#ICA
ICA <- function(data){
  t0 <- Sys.time()
  ica_res <- fastICA(data, 2, alg.typ = "parallel", fun = "logcosh", alpha = 1, 
          method = "C", 
          row.norm = FALSE, maxit = 20, 
          tol = 1e-5, verbose = FALSE)
  ica_res <- as.data.frame(scale(ica_res$S))
  ica_res <- ica_res[,1:2]
  t_end <- Sys.time()
  exeTime <- difftime(t_end, t0, units='sec')
  return(list(ica_res,exeTime))
}

#tSNE
tSNE <- function(data){
  t0 <- Sys.time()
  results <- Rtsne(data, dims = 2, perplexity = round(nrow(data)/5), max_iter = 1000, verbose = FALSE)
  results <- as.data.frame(scale(results$Y))
  results <- results[,1:2]
  t_end <- Sys.time()
  exeTime <- difftime(t_end, t0, units='sec')
  return(list(results,exeTime))
}

#tSNE_dopar
tSNE_dopar <- function(data){
  t0 <- Sys.time()
  cl <- makeCluster(detectCores())
  registerDoParallel(cl)
  #clusterEvalQ(cl, library(Rtsne))
  #clusterEvalQ(cl, library(snow))
  tsne <- Rtsne(data, dims = 2, perplexity = round(nrow(data)/10), max_iter = 1000, verbose = FALSE, do_parallel = TRUE)
  tsne <- as.data.frame(scale(tsne$Y))
  tsne <- tsne[,1:2]
  stopCluster(cl)
  t_end <- Sys.time()
  exeTime <- difftime(t_end, t0, units='sec')
  return(list(tsne,exeTime))
}

#tSNE_snow
tSNE_snow <- function(data){
  t0 <- Sys.time()
  cl <- makeCluster(detectCores())
  #clusterEvalQ(cl, library(Rtsne))
  #clusterEvalQ(cl, library(snow))
  tsne <- Rtsne(data, dims = 2, perplexity = round(nrow(data)/10), max_iter = 1000, verbose = FALSE, parallel = TRUE, num_threads = length(cl))
  tsne <- as.data.frame(scale(tsne$Y))
  tsne <- tsne[,1:2]
  stopCluster(cl)
  t_end <- Sys.time()
  exeTime <- difftime(t_end, t0, units='sec')
  return(list(tsne,exeTime))
}

#UMAP
UMAP <- function(data){
  t0 <- Sys.time()
  umap_result <- umap(data, metric = 'cosine')
  umap_result <- as.data.frame(scale(as.matrix(umap_result$layout)))
  umap_result <- umap_result[,1:2]
  t_end <- Sys.time()
  exeTime <- difftime(t_end, t0, units='sec')
  return(list(umap_result,exeTime))
}

#UMAP_n_thread
UMAP_n_thread <- function(data){
  t0 <- Sys.time()
  umap_result <- umap(data, metric = 'cosine', n_threads = 8)
  umap_result <- as.data.frame(scale(as.matrix(umap_result$layout)))
  umap_result <- umap_result[,1:2]
  t_end <- Sys.time()
  exeTime <- difftime(t_end, t0, units='sec')
  return(list(umap_result,exeTime))
}

#UMAP_snow
UMAP_snow <- function(data){
  t0 <- Sys.time()
  cl <- makeCluster(detectCores())
  clusterEvalQ(cl, {
    library(umap)
  })
  data_split <- split(data, 1:detectCores())
  umap_results <- clusterApply(cl, data_split, function(x) umap(as.matrix(x)))
  umap_results <- do.call(rbind, umap_results)
  umap_results <- umap_results[1:detectCores()]
  umap_results <- as.data.frame(do.call(rbind, umap_results))
  umap_results <- umap_results[order(as.integer(rownames(umap_results))), ]
  umap_results <- as.data.frame(scale(as.matrix(umap_results$layout)))
  umap_results <- umap_results[,1:2]
  stopCluster(cl)
  t_end <- Sys.time()
  exeTime <- difftime(t_end, t0, units='sec')
  return(list(umap_results,exeTime))
}

#UMAP_dopar
UMAP_dopar <- function(data){
  t0 <- Sys.time()
  cl <- makeCluster(detectCores())
  registerDoParallel(cl)
  umap_result <- umap(data, metric = 'cosine', n_threads = 8)
  umap_result <- as.data.frame(scale(as.matrix(umap_result$layout)))
  umap_result <- umap_result[,1:2]
  stopCluster(cl)
  t_end <- Sys.time()
  exeTime <- difftime(t_end, t0, units='sec')
  return(list(umap_result,exeTime))
}

#MDS
MDS <- function(data)
{
  t0 <- Sys.time()
  N = nrow(data)
  labels <- labels
  data_mds <- data/N
  k = 2
  D <- rdist(data_mds)^2 
  D <- as.matrix(D)
  I <- diag(nrow(data_mds))
  i <- rep(1,nrow(data_mds))
  J <- I - 1/nrow(data_mds) * i %*% t(i)
  B <- -1/2*J%*%D%*%J 
  Eig_B <-  eigen(B)
  E <- Eig_B$vectors 
  Delta <-  diag(Eig_B$values)
  mds <- data.frame( 'mds_1' = E[,1] * sqrt(Delta[1,1]), 'mds_2'=E[,2] * sqrt(Delta[2,2]))
  mds <- as.data.frame(scale(mds))
  mds <- mds[,1:2]
  t_end <- Sys.time()
  exeTime <- difftime(t_end, t0, units='sec')
  return(list(mds,exeTime))
}


#NMF
NMF <- function(data){
  t0 <- Sys.time()
  results <- nmf(data,2)
  results <- basis(results)
  results <- as.data.frame(scale(results))
  results <- results[,1:2]
  t_end <- Sys.time()
  exeTime <- difftime(t_end, t0, units='sec')
  return(list(results,exeTime))
}


#H2O autoencoder

H2O <- function(data){
  
  
  # Инициализация h2o
  h2o.init(nthreads = detectCores())
  t0 <- Sys.time()
  # Преобразование в h2o.frame
  df_h2o <- as.h2o(data)
  
  # Разделение на обучающий и тестовый наборы
  split <- h2o.splitFrame(df_h2o, ratios = c(0.8), seed = 42)
  train <- split[[1]]
  test <- split[[2]]
  
  # Определение размерности входных данных
  input_dim <- ncol(data)
  
  # Определение сетки гиперпараметров
  param_grid <- expand.grid(
    units1 = c(round(input_dim/2), round(input_dim/1.5)),
    units2 = c(round(input_dim/4), round(input_dim/3)),
    activation = c("relu", "tanh"),
    dropout = c(0.5)
  )
  
  # Функция для создания модели автоэнкодера с заданными гиперпараметрами
  create_model <- function(units1, units2, activation, dropout) {
    model <- keras_model_sequential()
    model %>%
      layer_dense(units = units1, activation = activation, input_shape = c(input_dim)) %>%
      layer_dropout(rate = dropout) %>%
      layer_dense(units = units2, activation = activation) %>%
      layer_dropout(rate = dropout) %>%
      layer_dense(units = units1, activation = activation) %>%
      layer_dropout(rate = dropout) %>%
      layer_dense(units = input_dim, activation = "sigmoid")
    model %>% compile(optimizer = "adam", loss = "mse")
    return(model)
  }
  
  # Функция для обучения модели с заданными гиперпараметрами
  train_model <- function(model, train, test, epochs = 50, batch_size = 32) {
    history <- model %>% fit(
      x = as.matrix(train[,1:input_dim]),
      y = as.matrix(train[,1:input_dim]),
      epochs = epochs,
      batch_size = batch_size,
      validation_data = list(as.matrix(test[,1:input_dim]), as.matrix(test[,1:input_dim])),
      workers = 4
    )
    score <- model %>% evaluate(
      x = as.matrix(test[,1:input_dim]),
      y = as.matrix(test[,1:input_dim])
    )
    return(score[[1]])
  }
  
  # Инициализация списка результатов
  results_df <- data.frame()
  
  # Обучение и оценка модели с каждой комбинацией гиперпараметров
  for (i in 1:nrow(param_grid)) {
    cat(sprintf("Training model %d of %d\n", i, nrow(param_grid)))
    params <- param_grid[i, ]
    model <- create_model(params$units1, params$units2, params$activation, params$dropout)
    score <- train_model(model, train, test)
    result <- c(params, score)
    names(result) <- c("units1", "units2", "activation", "dropout", "loss")
    results_df <- results_df %>% rbind(result)
  }
  
  #write.csv(results_df, file = "output/h20_results.csv")
  
  #Поиск наилучшей комбинации гиперпараметров
  best_params <- results_df %>% slice(which.min(loss))
  
  #Обучение модели с лучшей комбинацией гиперпараметров
  best_model <- create_model(best_params$units1, best_params$units2, best_params$activation, best_params$dropout)
  best_model %>% fit(
    x = as.matrix(df_h2o[,1:input_dim]),
    y = as.matrix(df_h2o[,1:input_dim]),
    epochs = 50,
    batch_size = 32
  )
  
  
  #Извлечение скрытых представлений для визуализации
  hidden_layer_model <- keras_model(inputs = best_model$input, outputs = best_model$layers[[3]]$output)
  hidden_layer <- predict(hidden_layer_model, as.matrix(df_h2o[,1:input_dim]))
  hidden_layer_df <- as.data.frame(hidden_layer)
  
  # tsne_model <- Rtsne(hidden_layer, dims = 2, perplexity = 45)
  # tsne_result <- as.data.frame(tsne_model$Y)
  
  results <- UMAP_n_thread(hidden_layer)
  results <- results[[1]]
  results <- as.data.frame(scale(results))
  results <- results[,1:2]
  h2o.shutdown()
  t_end <- Sys.time()
  exeTime <- difftime(t_end, t0, units='sec')
  return(list(results,exeTime))
}



#Bagging
Bagging <- function(data, Clusters){
  t0 <- Sys.time()
  bootstrap <- function(data, Clusters) {
    n <- nrow(data)/2
    df_random <- data %>% mutate(row_id = row_number(),
                                 Clusters = Clusters) %>% sample_n(n)
  }
  
  data_bootstraps <- replicate(20, bootstrap(data, Clusters), simplify = FALSE)
  
  
  # Define the UMAP reduction function
  umap_reduction <- function(data) {
    umap(data, metric = "cosine", n_threads = 8)
  }
  
  
  # Perform UMAP on each bootstrap sample
  data_bootstraps_umap <- map(data_bootstraps, function(x) {
    umap_reduction(x %>% dplyr::select(-c(row_id,Clusters)))$layout %>% 
      as_tibble() %>% 
      rename_all(~paste0("UMAP", seq_along(.))) %>% 
      bind_cols(label = x$Clusters) %>%
      bind_cols(row_id = x$row_id)
  })
  
  # Combine the results into a single data frame
  data_reduced <- bind_rows(data_bootstraps_umap)
  
  # Aggregate the results to obtain the final predictions
  data_reduced_agg <- as.data.frame(data_reduced %>% 
                                      group_by(row_id,labels) %>% 
                                      summarize(across(starts_with("UMAP"), mean)) %>% 
                                      ungroup() %>%
                                      arrange(row_id)) %>%
    dplyr::select(-row_id) %>%
    dplyr::select(UMAP1,UMAP2,labels)
  
  data_reduced_agg <- data_reduced_agg %>%
    mutate(UMAP1 = scale(UMAP1),
           UMAP2 = scale(UMAP2))
  
  t_end <- Sys.time()
  exeTime <- difftime(t_end, t0, units='sec')
  return(list(data_reduced_agg,exeTime))
  
}


