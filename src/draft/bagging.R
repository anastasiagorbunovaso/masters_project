data <- read.csv("data/gene_expression_data_system_1.csv")
Clusters <- read.csv("data/labels_system_1.csv")
Clusters <- as.factor(Clusters$data_gen.labels)


bootstrap <- function(data) {
  n <- nrow(data)/2
  df_random <- data %>% mutate(row_id = row_number(),
                               Clusters = Clusters) %>% sample_n(n)
}

data_bootstraps <- replicate(20, bootstrap(data), simplify = FALSE)


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
  group_by(row_id,label) %>% 
  summarize(across(starts_with("UMAP"), mean)) %>% 
  ungroup() %>%
  arrange(row_id))

# # Print the reduced data
# print(data_reduced_agg)


# Plot the UMAP reduced data with bagging
ggplot(data_reduced_agg, 
       aes(x = data_reduced_agg[,3], y = data_reduced_agg[,4],color=Clusters))+geom_point(size = 2)+
  labs(title = "Bagging", 
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
