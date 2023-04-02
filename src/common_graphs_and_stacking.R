output_table <- fread(paste0(output_file_dir,file_to_work,'/',"output_table_",file_to_work,".csv"))

rm(data)
gc()

all_alg_names <- c('PCA','PCA_dopar','PCA_snow','ICA','tSNE','tSNE_dopar','tSNE_snow','UMAP',
                   'UMAP_n_thread','UMAP_snow','UMAP_dopar','MDS','NMF','H2O','Bagging')

stacking_info_table <- data.frame()
ALL_res <- data.frame()

for(nm in all_alg_names){
  print(nm)
  
  tmp <- fread(paste0(output_file_dir,file_to_work,'/',model_dir,nm,"_res_",file_to_work,".csv"))
  tmp <- tmp[,1:2]
  colnames(tmp) <- c(paste0(nm,'_1'),paste0(nm,'_2'))
  
  if(nrow(ALL_res) == 0){
    ALL_res <- ALL_res %>% rbind(tmp)
  }else{
    ALL_res <- ALL_res %>% cbind(tmp)
  }
  
}

ALL_res_cor <- round(cor(ALL_res),2)
corrp.mat <- cor_pmat(ALL_res_cor)

plot <- ggcorrplot(ALL_res_cor, hc.order =TRUE, type ="lower", 
                   p.mat = corrp.mat) + 
  ggtitle("Сorrelation Matrix") +
  theme(plot.title = element_text(hjust = 0.5, size = 35))
png(paste0(output_file_dir,file_to_work,'/',"cor_plot_fo_stacking.png"), width = 1000, height = 600, res = 96)
print(plot)
dev.off()


# Расчет корреляции между парами переменных
corr_matrix <- round(cor(ALL_res),2)
corr_matrix[lower.tri(corr_matrix)] <- 0
diag(corr_matrix) <- 0

# Удаляем коррелирующие столбцы
corr_threshold <- 0.8 # порог корреляции
for_del_cols <- rownames(corr_matrix)[apply(abs(corr_matrix) >= corr_threshold, 1, any)]

ALL_res_new <- ALL_res %>% dplyr::select(-for_del_cols)

stacking_info_table <- rbind(stacking_info_table, c(ncol(ALL_res),paste(colnames(ALL_res), collapse = ", "),
                                                  ncol(ALL_res_new),paste(colnames(ALL_res_new), collapse = ", "),
                                                  paste(setdiff(colnames(ALL_res), colnames(ALL_res_new)), collapse = ", ")))
colnames(stacking_info_table) <- c('cols_before','colnames_before','cols_after','colnames_after','DIFF')

write.xlsx(stacking_info_table,paste0(output_file_dir,file_to_work,'/',"stacking_info_table_",file_to_work,".xlsx"))



#Stacking
alg_name <- 'Stacking'
res_alg <- UMAP_n_thread(ALL_res_new)

output_table_tmp <- as.data.frame(output_table[,1:(ncol(output_table)-2)])
output_table_tmp <- save_func(res_alg,alg_name,Clusters,output_table_tmp)

output_table <- output_table_tmp %>%
  mutate(ACCURACY_INIT = output_table$ACCURACY_INIT[1],
         TIME_RF_INIT = output_table$TIME_RF_INIT[1])

write_csv(output_table,paste0(output_file_dir,file_to_work,'/',"output_table_",file_to_work,".csv"))
write.xlsx(output_table,paste0(output_file_dir,file_to_work,'/',"output_table_",file_to_work,".xlsx"))


output_table <- output_table %>%
  mutate(TIME = as.double(TIME),
         Q1 = as.double(Q1),
         Q2 = as.double(Q2),
         ACCURACY = as.double(ACCURACY),
         TIME_RF = as.double(TIME_RF))


###all common graphs

#Q1
Q1_res <- output_table %>%
  arrange(Q1)

res <- ggplot(Q1_res, aes(x = reorder(ALGORITHM, Q1), y = Q1)) +
  geom_col(fill = "darkblue") +
  labs(title = "Ratio of average intra-cluster and inter-cluster distances", 
       x = "Algorithm")+
  theme(plot.title=element_text(hjust = 0.5))+
  theme(axis.text.y = element_text(size = 20, face = "italic")) +
  theme(axis.text.x = element_text(size = 20, face = "italic"))+
  theme(plot.title=element_text(hjust = 0.5))+
  theme(legend.title = element_text(size = 15)) +
  theme(legend.text = element_text(size = 10)) +
  guides(color = guide_legend(override.aes = list(size = 5)))+
  theme(plot.title = element_text(size = 35),
        axis.title.x = element_text(size = 30),
        axis.title.y = element_text(size = 30))

png(paste0(output_file_dir,file_to_work,'/',"Q1_res.png"), width = 1600, height = 1000, res = 96)
print(res)
dev.off()


#Q2
Q2_res <- output_table %>%
  arrange(Q2)

res <- ggplot(Q2_res, aes(x = reorder(ALGORITHM, Q2), y = Q2)) +
  geom_col(fill = "darkblue") +
  labs(title = "Ratio of sum of sq. intra-cluster and inter-cluster distances", 
       x = "Algorithm")+
  theme(plot.title=element_text(hjust = 0.5))+
  theme(axis.text.y = element_text(size = 20, face = "italic")) +
  theme(axis.text.x = element_text(size = 20, face = "italic"))+
  theme(plot.title=element_text(hjust = 0.5))+
  theme(legend.title = element_text(size = 15)) +
  theme(legend.text = element_text(size = 10)) +
  guides(color = guide_legend(override.aes = list(size = 5)))+
  theme(plot.title = element_text(size = 35),
        axis.title.x = element_text(size = 30),
        axis.title.y = element_text(size = 30))

png(paste0(output_file_dir,file_to_work,'/',"Q2_res.png"), width = 1600, height = 1000, res = 96)
print(res)
dev.off()


#TIME
TIME <- output_table %>%
  arrange(TIME)

res <- ggplot(TIME, aes(x = reorder(ALGORITHM, TIME), y = TIME)) +
  geom_col(fill = "darkblue") +
  labs(title = "Algorithm execution time", 
       x = "Algorithm")+
  theme(plot.title=element_text(hjust = 0.5))+
  theme(axis.text.y = element_text(size = 20, face = "italic")) +
  theme(axis.text.x = element_text(size = 20, face = "italic"))+
  theme(plot.title=element_text(hjust = 0.5))+
  theme(legend.title = element_text(size = 15)) +
  theme(legend.text = element_text(size = 10)) +
  guides(color = guide_legend(override.aes = list(size = 5)))+
  theme(plot.title = element_text(size = 35),
        axis.title.x = element_text(size = 30),
        axis.title.y = element_text(size = 30))

png(paste0(output_file_dir,file_to_work,'/',"TIME.png"), width = 1600, height = 1000, res = 96)
print(res)
dev.off()


#ACC_RF
ACC_RF <- output_table %>%
  dplyr::select(ALGORITHM,ACCURACY) %>%
  rbind(c('INIT',output_table$ACCURACY_INIT[1])) %>%
  mutate(ACCURACY = as.double(ACCURACY)) %>%
  arrange(desc(ACCURACY)) %>%
  mutate(ACCURACY = ACCURACY*100)
  

res <- ggplot(ACC_RF, aes(x = reorder(ALGORITHM, desc(ACCURACY)), y = ACCURACY)) +
  geom_col(fill = "darkblue") +
  labs(title = "Random Forest Accuracy", 
       x = "Algorithm")+
  theme(plot.title=element_text(hjust = 0.5))+
  theme(axis.text.y = element_text(size = 20, face = "italic")) +
  theme(axis.text.x = element_text(size = 20, face = "italic"))+
  theme(plot.title=element_text(hjust = 0.5))+
  theme(legend.title = element_text(size = 15)) +
  theme(legend.text = element_text(size = 10)) +
  guides(color = guide_legend(override.aes = list(size = 5)))+
  theme(plot.title = element_text(size = 35),
        axis.title.x = element_text(size = 30),
        axis.title.y = element_text(size = 30))

png(paste0(output_file_dir,file_to_work,'/',"ACC_RF.png"), width = 1600, height = 1000, res = 96)
print(res)
dev.off()



#TIME_RF
TIME_RF <- output_table %>%
  dplyr::select(ALGORITHM,TIME_RF) %>%
  rbind(c('INIT',output_table$TIME_RF_INIT[1])) %>%
  mutate(TIME_RF = as.double(TIME_RF)) %>%
  arrange(TIME_RF)


res <- ggplot(TIME_RF, aes(x = reorder(ALGORITHM, TIME_RF), y = TIME_RF)) +
  geom_col(fill = "darkblue") +
  labs(title = "Random Forest Time", 
       x = "Algorithm")+
  theme(plot.title=element_text(hjust = 0.5))+
  theme(axis.text.y = element_text(size = 20, face = "italic")) +
  theme(axis.text.x = element_text(size = 20, face = "italic"))+
  theme(plot.title=element_text(hjust = 0.5))+
  theme(legend.title = element_text(size = 15)) +
  theme(legend.text = element_text(size = 10)) +
  guides(color = guide_legend(override.aes = list(size = 5)))+
  theme(plot.title = element_text(size = 35),
        axis.title.x = element_text(size = 30),
        axis.title.y = element_text(size = 30))

png(paste0(output_file_dir,file_to_work,'/',"TIME_RF.png"), width = 1600, height = 1000, res = 96)
print(res)
dev.off()
