library(snow)

# data <- data %>%
#   mutate(label = iris[,5])

PCA_par_snow <- function(data) {
  
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
    Z <- as.data.frame(Z)
    Z <- cbind(Z,labels)
    # Return Z
    return(Z)
  })
  
  # Stop cluster
  stopCluster(cl)
  
  # Combine results
  full_projected_data <- do.call(rbind, Z)
  
  # Return projected data
  full_projected_data
  
}

num_cores <- 3
data <- fread('data/data_mRNA_BRCA.csv')
Clusters <- fread('data/labels_mRNA_BRCA.csv')
Clusters <- as.factor(Clusters$labels)

t0 <- proc.time()
Z <- PCA_par(data)
exeTimePCA <- proc.time() - t0
exeTimePCA
Z <- as.data.frame(Z)
ggplot(data = Z, aes(x = Z[,1], y = Z[,2],color=Clusters))+geom_point(size = 2)+
  labs(title = "PCA", 
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

t0 <- proc.time()
Z <- PCA(data)
exeTimePCA <- proc.time() - t0
exeTimePCA
Z <- as.data.frame(Z)
ggplot(data = Z, aes(x = Z[,1], y = Z[,2],color=Clusters))+geom_point(size = 2)+
  labs(title = "PCA", 
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