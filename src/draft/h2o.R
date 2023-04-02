# Загрузка библиотек
library(h2o)
library(keras)
library(tidyverse)
library(gridExtra)
library(Rtsne)

# Инициализация h2o
h2o.init(nthreads = 4)

# Загрузка датафрейма
df <- read.csv("data/gene_expression_data_system_1.csv")
Clusters <- read.csv("data/labels_system_1.csv")
Clusters <- as.factor(Clusters$data_gen.labels)

# Преобразование в h2o.frame
df_h2o <- as.h2o(df)

# Разделение на обучающий и тестовый наборы
split <- h2o.splitFrame(df_h2o, ratios = c(0.8), seed = 42)
train <- split[[1]]
test <- split[[2]]

# Определение размерности входных данных
input_dim <- ncol(df)

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

write.csv(results_df, file = "output/h20_results.csv")

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

tsne_model <- Rtsne(hidden_layer, dims = 2, perplexity = 45)
tsne_result <- as.data.frame(tsne_model$Y)

res_plot <- ggplot(data = tsne_result, aes(x = tsne_result[,1], y = tsne_result[,2],color=Clusters))+geom_point(size = 2)+
  labs(title = "H2O", 
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

res_plot

ggsave(paste0(output_file_dir,"h2o.jpg"), plot = res_plot, width = 10, height = 7, dpi = 600)
png(paste0(output_file_dir,"h2o.png"), width = 1000, height = 700)
ggplot(data = tsne_result, aes(x = tsne_result[,1], y = tsne_result[,2],color=Clusters))+geom_point(size = 2)+
  labs(title = "H2O", 
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
dev.off()

h2o.shutdown()
