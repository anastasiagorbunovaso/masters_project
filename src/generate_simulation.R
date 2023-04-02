generate_expression_data <- function(n_samples, n_genes, n_clusters, noise_level) {
  
  # set parameters
  cluster_sizes <- sample(ceiling(n_samples/n_clusters):ceiling(n_samples/(n_clusters)), n_clusters-1, replace = TRUE)
  cluster_sizes[n_clusters] = n_samples - sum(cluster_sizes)
  
  # generate cluster centers
  cluster_centers <- matrix(rnorm(n_clusters * n_genes), ncol = n_genes)
  
  # generate samples from clusters
  samples <- matrix(0, nrow = n_samples, ncol = n_genes)
  start_idx <- 1
  for (i in 1:n_clusters) {
    end_idx <- start_idx + cluster_sizes[i] - 1
    samples[start_idx:end_idx, ] <- mvrnorm(cluster_sizes[i], mu = cluster_centers[i, ], Sigma = diag(n_genes))
    start_idx <- end_idx + 1
  }
  
  # add noise
  noise <- matrix(rnorm(n_samples * n_genes, sd = noise_level), nrow = n_samples, ncol = n_genes)
  samples <- samples + noise
  
  # create labels
  labels <- rep(1:n_clusters, times = cluster_sizes)
  
  # return data and labels
  return(list(data = samples, labels = labels))
}



#Установка параметров
n_samples <- 2000 # Количество образцов
n_genes <- 2000 # Количество генов
n_clusters <- 5 # Количество кластеров
noise_levels <- c(0.1, 4, 6)

# generate data
data_list <- list()
for (i in 1:length(noise_levels)) {
  noise_level <- noise_levels[i]
  data_list[[i]] <- generate_expression_data(n_samples, n_genes, n_clusters, noise_level)
  data_gen <- data_list[[i]]
  data_gen_data <- as.data.frame(data_gen$data)
  data_gen_labels <- as.data.frame(data_gen$labels)
  data_gen_data <- data_gen_data %>% 
    mutate_all(~ (.-min(.))/(max(.)-min(.)) * 100)
  # Сохранение данных в файл
  fwrite(data_gen_data, file = paste0(data_dir,"gene_expression_data_system_", i, ".csv"))
  fwrite(data_gen_labels, file = paste0(data_dir,"labels_system_", i, ".csv"))
}

