#########################################
### Author: Friederike Marie Moroff
### Date: 15/10/22
#########################################

### Loading any libraries ----
library(corrr)
library(usethis) 
library(dplyr)
library(caret)
library(ellipse)
library(MASS)
library(ISLR)
library(tidyverse)
library(ROSE)

### change memory limit ----
# usethis::edit_r_environ()
# -> save R_MAX_VSIZE=100Gb 


### set the working directory ----
rm(list=ls())
setwd("~/Documents/Universität/Erasmus_kurse/big_data/project_shared/ECNumberPred/data")


### load the data and set seed ----
set.seed(1)

data_counts = readRDS("TRAINING_EC.rds")            # dataframe with the aminoacid counts
data_features = readRDS("TRAINING_descriptors.rds") # dataframe with the features
rownames(data_features) = data_features$id          # set rownames
data_features = data_features[,-1]                  # delete id column
data_counts$EC <- as.factor(data_counts$EC)         # change type of EC variable to factor
data_features$EC <- as.factor(data_features$EC)

# merge both datasets
data_merged = merge(data_counts, data_features)


### inspect the data ----

# dimension of the data set
dim(data_counts)   # -> 14956  instances and 21 attributes
dim(data_features) # -> 14956  instances and 13 attributes
dim(data_merged)   

# list the levels for the class
levels(data_counts$EC)
# "Hydrolases"      "Isomerases"      "Ligases"         "Lyases"         
# "Oxidoreductases" "Transferases" 
# -> same for every dataframe

# class distribution
percentage <- prop.table(table(data_counts$EC)) * 100
cbind(freq=table(data_counts$EC), percentage=percentage)
#freq percentage
#Hydrolases      2497   16.69564
#Isomerases      2493   16.66890
#Ligases         2499   16.70901
#Lyases          2498   16.70233
#Oxidoreductases 2480   16.58197
#Transferases    2489   16.64215
# -> same for every dataframe
# -> data is balanced


### check for colinearity ----
corr_1 <- data_counts %>%
  correlate() %>%    # Create correlation data frame (cor_df)
  rearrange() %>%    # rearrange by correlations
  shave()            # Shave off the upper triangle for a clean result

fashion(corr_1)      # -> highest correlation between K and S with -.58
rplot(corr_1)

corr_2 <- data_features %>%
  correlate() %>%    # Create correlation data frame (cor_df)
  rearrange() %>%    # rearrange by correlations
  shave()            # Shave off the upper triangle for a clean result

fashion(corr_2)      # -> highest correlation between Non and Polar with -1.00,
rplot(corr_2)        # Charged and Acidic with -.83, Charged and Basic with -.79


### evaluate algorithms ----

control <- trainControl(method="cv", number=5)
metric <- "Accuracy"

# dataframe with aminoacids counts
fit.lda <- train(EC~., data=data_counts, method="lda", metric=metric, trControl=control)
fit.cart <- train(EC~., data=data_counts, method="rpart", metric=metric, trControl=control)
fit.knn <- train(EC~., data=data_counts, method="knn", metric=metric, trControl=control)
fit.svm <- train(EC~., data=data_counts, method="svmRadial", metric=metric, trControl=control)

results_counts <- resamples(list(lda=fit.lda, cart=fit.cart, knn=fit.knn, svm=fit.svm))
summary(results_counts)
#Accuracy 
#lda  0.3001337
#cart 0.2360415 
#knn  0.3580742 
#svm  0.4035440

# dataframe with features
fit.lda <- train(EC~., data=data_features, method="lda", metric=metric, trControl=control)
fit.cart <- train(EC~., data=data_features, method="rpart", metric=metric, trControl=control)
fit.knn <- train(EC~., data=data_features, method="knn", metric=metric, trControl=control)
fit.svm <- train(EC~., data=data_features, method="svmRadial", metric=metric, trControl=control)

results_features <- resamples(list(lda=fit.lda, cart=fit.cart, knn=fit.knn, svm=fit.svm))
summary(results_features)
#Accuracy 
#lda  0.2844251 
#cart 0.2515870 
#knn  0.2304348 
#svm  0.3620863

# both dataframes
fit.lda <- train(EC~., data=data_merged, method="lda", metric=metric, trControl=control)
fit.cart <- train(EC~., data=data_merged, method="rpart", metric=metric, trControl=control)
fit.knn <- train(EC~., data=data_merged, method="knn", metric=metric, trControl=control)
fit.svm <- train(EC~., data=data_merged, method="svmRadial", metric=metric, trControl=control)

results_merged <- resamples(list(lda=fit.lda, cart=fit.cart, knn=fit.knn, svm=fit.svm))
summary(results_merged)

###hsuwzs

