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
#setwd("~/Documents/Universität/Erasmus_kurse/big_data/project_shared/ECNumberPred/data")


### load the data and set seed ----
set.seed(1)

data_counts = readRDS("data/TRAINING_EC.rds")            # dataframe with the aminoacid counts
data_features = as.data.frame(readRDS("data/TRAINING_descriptors.rds"))# dataframe with the features
data_APAAC = readRDS("data/TRAINING_AA_DESC_APAAC.rds") #df with APAAC features


rownames(data_features) = data_features$id          # set rownames : 
data_features = data_features[,-1]                  # delete id column


# merge both datasets
data_merged = merge(data_counts, data_features, by =0) # by = 0 indicates by rowname

all(data_merged$EC.x == data_merged$EC.y) # check if EC is preserved persitantly

data_merged<-data_merged %>%   # remove doubled EC column
  dplyr::select(-EC.y) %>%
  dplyr::rename(EC=EC.x)%>% 
  column_to_rownames("Row.names") # re-set rownames


# set EC as factor
data_counts$EC <- as.factor(data_counts$EC)         
data_features$EC <- as.factor(data_features$EC)
data_APAAC$EC <- as.factor(data_APAAC$EC)
data_merged$EC <- as.factor(data_merged$EC)

### inspect the data ----

# dimension of the data set
dim(data_counts)   # -> 14956  instances and 21 attributes
dim(data_features) # -> 14956  instances and 13 attributes
dim(data_merged)   # -> 14956  ;  33
dim(data_APAAC)    # 14956    81
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


### Train ML algorithms ----
## NOTE: crossvalidation for support vector machine seems to be bugged,
## even though long runtime is expected since space complexity is O(n^2) ,still takes way to long ,
## see: https://stackoverflow.com/questions/30385347/r-caret-unusually-slow-when-tuning-svm-with-linear-kernel

## Params:
control <- trainControl(method="cv", number=5)
metric <- "Accuracy"

# 1: Aminoacids counts as only features

## Preprocessing - does not influence performance here
processing <- preProcess(data_counts, method = c("nzv","center","scale","corr")) 
data_counts <- predict(processing, data_counts) 

## training
fit.lda_AAC <- train(EC~., data=data_counts, method="lda", metric=metric, trControl=control)
fit.cart_AAC <- train(EC~., data=data_counts, method="rpart", metric=metric, trControl=control)
fit.knn_AAC <- train(EC~., data=data_counts, method="knn", metric=metric, trControl=control)
fit.svm_AAC <- train(EC~., data=data_counts, method="svmRadial", metric=metric, trControl=control)

models_AAC<-list(lda=fit.lda_AAC, cart=fit.cart_AAC, knn=fit.knn_AAC, svm=fit.svm_AAC)
saveRDS(models_AAC,"results/models_AAC.rds")
## evaluation
results_counts <- resamples(models_AAC)
s<-summary(results_counts)
s$statistics$Accuracy[,"Median"]
# lda      cart       knn       svm 
# 0.2938816 0.2360415 0.3711133 0.4035440 
saveRDS(results_counts,"results/AAC_counts")

# 2: simple AA features as model features

processing <- preProcess(data_features, method = c("nzv","center","scale","corr")) 
data_features <- predict(processing, data_features) 

fit.lda_features <- train(EC~., data=data_features, method="lda", metric=metric, trControl=control)
fit.cart_features <- train(EC~., data=data_features, method="rpart", metric=metric, trControl=control)
fit.knn_features <- train(EC~., data=data_features, method="knn", metric=metric, trControl=control)
fit.svm_features <- train(EC~., data=data_features, method="svmRadial", metric=metric,trControl=control)

models_features<-list(lda=fit.lda_features, cart=fit.cart_features, knn=fit.knn_features, svm=fit.svm_features)
saveRDS(models_features,"results/models_features.rds")

results_features <- resamples(models_features)
s<-summary(results_features)
s$statistics$Accuracy[,"Median"]

# lda      cart       knn       svm 
# 0.2782609 0.2573529 0.3395785 0.3636364 

saveRDS(results_features,"results/simple_features")

# 3: merged dataframe with AAC and simple features

processing <- preProcess(data_merged, method = c("nzv","center","scale","corr")) 
data_merged <- predict(processing, data_merged)

fit.lda_merged <- train(EC~., data=data_merged, method="lda", metric=metric, trControl=control)
fit.cart_merged <- train(EC~., data=data_merged, method="rpart", metric=metric, trControl=control)
fit.knn_merged <- train(EC~., data=data_merged, method="knn", metric=metric, trControl=control)
fit.svm_merged <- train(EC~., data=data_merged, method="svmRadial", metric=metric, trControl=control)

models_merged<-list(lda=fit.lda_merged, cart=fit.cart_merged, knn=fit.knn_merged,svm=fit.svm_merged)
saveRDS(models_merged,"results/models_merged.rds")
results_merged <- resamples(models_merged)
s<-summary(results_merged)
s$statistics$Accuracy[,"Median"]
# lda      cart       knn       svm 
# 0.3141519 0.2733957 0.4070856 0.4194519 
saveRDS(results_merged,"results/simple_features_and_AAC_counts")

## Preprocessing
processing <- preProcess(data_APAAC, method = c("nzv","center","scale","corr"),cutoff =  0.8) 
data_APAAC <- predict(processing, data_APAAC) 

# 4 dataframe with APAAC descriptors as model features
fit.lda_APAAC <- train(EC~., data=data_APAAC, method="lda", metric=metric, trControl=control)
fit.cart_APAAC <- train(EC~., data=data_APAAC, method="rpart", metric=metric, trControl=control)
fit.knn_APAAC <- train(EC~., data=data_APAAC, method="knn", metric=metric, trControl=control)
fit.svm_APAAC <- train(EC~., data=data_APAAC, method="svmRadial", metric=metric, trControl=control)

models_APAAC<-list(lda=fit.lda_APAAC, cart=fit.cart_APAAC, knn=fit.knn_APAAC,svm=fit.svm_APAAC)
saveRDS(models_APAAC,"results/models_APAAC.rds")

results_APAAC <- resamples(models_APAAC)
s<-summary(results_APAAC)

s$statistics$Accuracy[,"Median"]
# > s$statistics$Accuracy[,"Median"]
# lda      cart       knn       svm 
# 0.3259779 0.2530926 0.3823529 0.4122367 
saveRDS(results_APAAC,"results/APAAC_descriptors")

### create boxplot ----




