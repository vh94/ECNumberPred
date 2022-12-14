library(ggplot2)

### set the working directory ----
#rm(list=ls())
#setwd("~/Documents/Master_CompBioCoimbra/BigData/PredECMain/")

results = readRDS("results/conf_matrices.rds")

models = c('lda', 'knn', 'svm',"rf")
dataframes = c('AAC', 'DESC', 'AAC_DESC',"APAAC")

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
#data_results$Accuracy

p <- ggplot(data_results, aes(x=dataframe, y=Accuracy, fill=model)) + 
  geom_boxplot() +
  facet_grid(.~dataframe, scale="free")
p

