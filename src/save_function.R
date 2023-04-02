save_func <- function(res_alg,alg_name,Clusters,output_table){
  
  for_save <- res_alg[[1]][,1:2] %>% cbind(Clusters)
  #save_res
  write_csv(for_save,paste0(output_file_dir,file_to_work,'/',model_dir,alg_name,"_res_",file_to_work,".csv"))
  write.xlsx(for_save,paste0(output_file_dir,file_to_work,'/',model_dir,alg_name,"_res_",file_to_work,".xlsx"))
  #visualization
  viz_dim_red(res_alg[[1]][,1:2],as.factor(Clusters$labels),alg_name,file_to_work)
  #metrics
  met <- metrics(res_alg[[1]][,1:2],Clusters)
  #rf_acc_res
  rf_res <- rf_acc(res_alg[[1]][,1:2],Clusters)
  #bind row to table
  output_table <- rbind(output_table,c(file_to_work,alg_name,res_alg[[2]],met[[1]],met[[2]],rf_res[[1]],rf_res[[2]]))
  
  return(output_table)
}
