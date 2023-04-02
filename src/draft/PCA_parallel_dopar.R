library(foreach)
library(doParallel)


PCA_par_dopar <- function(data) {
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
  
  # Delta_Z = foreach(j=1:ncol(Z), .combine="c") %dopar% {
  #   sd(Z[,j])
  # }
  # 
  # Delta_norm_x = foreach(j=1:ncol(norm_x), .combine="c") %dopar% {
  #   sd(norm_x[,j])
  # }
  # 
  # Sum_Z = sum(Delta_Z^2)
  # Sum_norm_x = sum(Delta_norm_x^2)
  # Alpha = matrix(0,ncol(data),1)
  # 
  # 
  # summa <- foreach(i = 1:ncol(data), .combine = "+") %dopar% {
  #   sum(Delta[i, i])
  # }
  # 
  # 
  # Alpha = foreach(i = 1:ncol(data), .combine = 'c') %dopar% {
  #   Delta[i,i]/summa
  # }
  # 
  # sum(Alpha[1],Alpha[2])
  
  stopCluster(cl)
  
  Z
}


t0 <- proc.time()
Z <- PCA(data)
exeTimePCA <- proc.time() - t0
exeTimePCA


t0 <- proc.time()
Z_par <- PCA_par_dopar(data)
exeTimePCA <- proc.time() - t0
exeTimePCA


t0 <- proc.time()
Z_par_snow <- PCA_par_snow(data)
exeTimePCA <- proc.time() - t0
exeTimePCA



Z <- scale(Z)
Z_par <- scale(Z_par)
Z_par_snow <- scale(Z_par_snow)


ggplot(data = as.data.frame(Z), aes(x = Z[,1], y = Z[,2],color=Clusters))+geom_point(size = 2)+
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



ggplot(data = as.data.frame(Z_par), aes(x = Z_par[,1], y = Z_par[,2],color=Clusters))+geom_point(size = 2)+
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


ggplot(data = as.data.frame(Z_par_snow), aes(x = Z_par_snow[,1], y = Z_par_snow[,2],color=Clusters))+geom_point(size = 2)+
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
