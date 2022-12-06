library(ggplot2)

### set the working directory ----
rm(list=ls())
setwd("~/Documents/Universität/Erasmus_kurse/big_data/project_shared/Ohne_Titel")

results = readRDS("results/conf_matrices.rds")

models = c('lda', 'cart', 'knn', 'svm')
dataframes = c('AAC', 'DESC', 'AAC_DESC')

data_results = data.frame(matrix(ncol = 3))
colnames(data_results) = c('model', 'Accuracy', 'dataframe')

for (dataframe in dataframes) {
  for (model in models) {
    
    temp = results[[dataframe]][[model]]
    acc = temp$overall[['Accuracy']]
    acc_low = temp$overall[['AccuracyLower']]
    acc_upp = temp$overall[['AccuracyUpper']]
    
    data_results = rbind(data_results, c(model, acc, dataframe))
    data_results = rbind(data_results, c(model, acc_low, dataframe))
    data_results = rbind(data_results, c(model, acc_upp, dataframe))
  }
}
data_results = data_results[2:nrow(data_results),]
data_results$Accuracy = as.numeric(data_results$Accuracy)
data_results$Accuracy

p <- ggplot(data_results, aes(x=dataframe, y=Accuracy, fill=model)) + 
  geom_boxplot() +
  facet_wrap(~dataframe, scale="free")
p

