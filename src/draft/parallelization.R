data <- as.data.frame(fread('data/data_miRNA_BRCA.csv'))
Clusters <- fread('data/labels_miRNA_BRCA.csv')
Clusters <- as.factor(Clusters$labels)


# data <- as.data.frame(fread('data/data_mRNA_BRCA.csv'))
# Clusters <- fread('data/labels_mRNA_BRCA.csv')
# Clusters <- as.factor(Clusters$labels)


# data <- fread('data/data_meDNA_BRCA.csv')
# Clusters <- fread('data/labels_meDNA_BRCA.csv')
# Clusters <- as.factor(Clusters$labels)

t0 <- proc.time()
umap_result <- umap(data, metric = 'cosine', n_threads = 8)
exeTimeUMAP <- proc.time() - t0
exeTimeUMAP

t0 <- proc.time()
umap_result <- umap(data,metric = 'cosine')
exeTimeUMAP <- proc.time() - t0
exeTimeUMAP
umap_ <- as.matrix(umap_result$layout)
umap_ <- scale(umap_)
umap_ <- as.data.frame(umap_)
ggplot(data = umap_, aes(x = umap_[,1], y = umap_[,2],color=Clusters))+geom_point(size = 2)+
  labs(title = "UMAP", 
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




# Загрузка данных
library(bigmemory)
bm <- read.big.matrix("data/data_miRNA_BRCA.csv", sep=",", header=TRUE)

# Подготовка кластера snow
library(snow)
cl <- makeCluster(8)
clusterEvalQ(cl, library(Rtsne))
clusterEvalQ(cl, memory.limit(10000))

# Разбивка на пакеты и распределение на ядра
n <- nrow(bm)
batch_size <- round(n/3)
batches <- seq(1, n, by=batch_size)
if (batches[length(batches)] != n) {
  batches <- c(batches, n)
}
num_batches <- length(batches) - 1
batch_indices <- lapply(1:num_batches, function(i) {batches[i]:(batches[i+1]-1)})
clusterExport(cl, list("bm", "batch_indices"))

# Запуск алгоритма t-SNE на каждом пакете данных
tsne_result <- vector("list", num_batches)
for (i in 1:num_batches) {
  tsne_result[[i]] <- clusterApply(cl, batch_indices[[i]], function(batch) {
    tsne_result_batch <- Rtsne(bm[batch,], dims=2, perplexity=40, verbose=TRUE)
    return(tsne_result_batch)
  })
}
tsne_result <- unlist(tsne_result)

# Сбор результатов
tsne_final <- do.call(rbind, tsne_result)
colnames(tsne_final) <- c("tsne1", "tsne2")

# Закрытие кластера
stopCluster(cl)