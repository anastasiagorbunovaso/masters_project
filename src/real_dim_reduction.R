if(grepl("simulation", file_to_work)){
  
  system_number <- tail(strsplit(file_to_work, "_")[[1]], n = 1)
  
  data <- read.csv(paste0(data_dir,"gene_expression_data_system_",system_number,".csv"))
  Clusters <- read.csv(paste0(data_dir,"labels_system_",system_number,".csv"))
  colnames(Clusters) <- 'labels'
  
}else{
  
  data <- read.csv(paste0(data_dir,"data_",file_to_work,"_preproc.csv"))
  Clusters <- read.csv(paste0(data_dir,"labels_",file_to_work,"_preproc.csv"))
  colnames(Clusters) <- 'labels'
  
}

output_table <- data.frame()

#init rf acc and time
rf_init <- rf_acc(data,Clusters)

#all algorithms
#PCA
alg_name <- 'PCA'
res_alg <- PCA(data)
output_table <- save_func(res_alg,alg_name,Clusters,output_table)

#PCA_dopar
alg_name <- 'PCA_dopar'
res_alg <- PCA_dopar(data)
output_table <- save_func(res_alg,alg_name,Clusters,output_table)

#PCA_snow
alg_name <- 'PCA_snow'
res_alg <- PCA_snow(data,Clusters)
output_table <- save_func(res_alg,alg_name,Clusters,output_table)

#ICA
alg_name <- 'ICA'
res_alg <- ICA(data)
output_table <- save_func(res_alg,alg_name,Clusters,output_table)

#tSNE
alg_name <- 'tSNE'
res_alg <- tSNE(data)
output_table <- save_func(res_alg,alg_name,Clusters,output_table)

#tSNE_dopar
alg_name <- 'tSNE_dopar'
res_alg <- tSNE_dopar(data)
output_table <- save_func(res_alg,alg_name,Clusters,output_table)

#tSNE_snow
alg_name <- 'tSNE_snow'
res_alg <- tSNE_snow(data)
output_table <- save_func(res_alg,alg_name,Clusters,output_table)

#UMAP
alg_name <- 'UMAP'
res_alg <- UMAP(data)
output_table <- save_func(res_alg,alg_name,Clusters,output_table)

#UMAP_n_thread
alg_name <- 'UMAP_n_thread'
res_alg <- UMAP_n_thread(data)
output_table <- save_func(res_alg,alg_name,Clusters,output_table)

#UMAP_snow
alg_name <- 'UMAP_snow'
res_alg <- UMAP_snow(data)
output_table <- save_func(res_alg,alg_name,Clusters,output_table)

#UMAP_dopar
alg_name <- 'UMAP_dopar'
res_alg <- UMAP_dopar(data)
output_table <- save_func(res_alg,alg_name,Clusters,output_table)

#MDS
alg_name <- 'MDS'
res_alg <- MDS(data)
output_table <- save_func(res_alg,alg_name,Clusters,output_table)

#NMF
alg_name <- 'NMF'
res_alg <- NMF(data)
output_table <- save_func(res_alg,alg_name,Clusters,output_table)

#H2O_autoencoder
alg_name <- 'H2O'
res_alg <- H2O(data)
output_table <- save_func(res_alg,alg_name,Clusters,output_table)

#Bagging
alg_name <- 'Bagging'
res_alg <- Bagging(data)
output_table <- save_func(res_alg,alg_name,Clusters,output_table)


output_table <- output_table %>% cbind(rf_init[[1]])
output_table <- output_table %>% cbind(as.double(rf_init[[2]]))
colnames(output_table) <- c('FILENAME','ALGORITHM','TIME','Q1','Q2','ACCURACY','TIME_RF','ACCURACY_INIT','TIME_RF_INIT')
  
write_csv(output_table,paste0(output_file_dir,file_to_work,'/',"output_table_",file_to_work,".csv"))
write.xlsx(output_table,paste0(output_file_dir,file_to_work,'/',"output_table_",file_to_work,".xlsx"))
