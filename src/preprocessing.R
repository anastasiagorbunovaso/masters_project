preproc_start <- Sys.time()

source("configs/get_config.R")
source("configs/check_load_file.R")

config <- get_config(data_dir)

cat(format(Sys.time(), "%b %d %Y, %X"), "|", 'Read input files...\n')
tables_to_load <- lapply(config$tables$input, check_load_file)

cat('All files checked!', "\n")
cat('-----------------------------------', "\n")

cat(format(Sys.time(), "%b %d %Y, %X"), "|", 'Start preproc...\n')


start <- Sys.time()

data <- as.data.frame(tables_to_load$data)
Clusters <- tables_to_load$labels

before_after_preproc <- data.frame()
before_after_preproc <- rbind(before_after_preproc, c(ncol(data),nrow(data),nrow(Clusters)))


if(file_to_work == 'meDNA_BRCA') {
  
  sd_data <- as.data.frame(apply(data, MARGIN = 2, FUN = sd))
  
  sd_data1 <- t(as.data.frame(sd_data$`apply(data, MARGIN = 2, FUN = sd)` >= quantile(sd_data$`apply(data, MARGIN = 2, FUN = sd)`, 0.8)))

  colnames(sd_data1) <- colnames(data)
  
  data <- rbind(data, sd_data1)
  data <- data[,(data[nrow(data),]>0)]
  data <- data[-nrow(data),]
  
}else{
  
  #remove features that are not expressed in at least 5% people
  i <- colSums(data > 0)
  is_expressed <- i > round(nrow(data)*0.05)
  data <- data[,is_expressed]
  
  #And people not expressing at least one UMI in at least 5% genes.
  i <- rowSums(data > 0)
  is_expressed <- i > round(ncol(data)*0.05)
  data <- data[is_expressed,]
  
  Clusters <- Clusters %>% cbind(is_expressed) %>%
    filter(is_expressed == T) %>% dplyr::select(-is_expressed)
  
}



before_after_preproc <- rbind(before_after_preproc, c(ncol(data),nrow(data),nrow(Clusters)))
colnames(before_after_preproc) <- c('ncol_data','nrow_data','nrow_labels')

fwrite(data,paste0(data_dir,'data_',file_to_work,'_preproc.csv'))
fwrite(Clusters,paste0(data_dir,'labels_',file_to_work,'_preproc.csv'))
write.xlsx(before_after_preproc,paste0(output_file_dir,file_to_work,'/before_after_preproc_',file_to_work,'.xlsx'))

cat(format(Sys.time(), "%b %d %Y, %X"), "|", 'End preproc.\n')

preproc_end <- Sys.time()

cat('Preprocessing: ', difftime(preproc_end, preproc_start, units='mins'), 'min', "\n")
time_list$Preprocessing <- difftime(preproc_end, preproc_start, units='mins')

