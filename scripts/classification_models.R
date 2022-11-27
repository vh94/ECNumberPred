#########################################
### Author: Friederike Marie Moroff
### Date: 15/10/22
#########################################

### Loading any libraries ----
install.packages("ROSE")

library(caret)
library(ellipse)
library(MASS)
library(ISLR)
library(tidyverse)
library(ROSE)


### set the working directory ----
rm(list=ls())
setwd("~/Documents/Universität/Erasmus_kurse/big_data/project")

### load the data and set seed ----

data = readRDS("comp_matrix_HSapiens.rds")
data_subset <- data[,2:22]

set.seed(1)

### inspect the data ----

# dimension of the data set
dim(data_subset)
# -> 4595 instances and 22 attributes

# list types for each attribute
sapply(data_subset, class)
#Ala       Cys       Asp       Glu       Phe       Gly       His       Ile       Lys       Leu 
#"integer" "integer" "integer" "integer" "integer" "integer" "integer" "integer" "integer" "integer" 
#Met       Asn       Pro       Gln       Arg       Ser       Thr       Val       Trp       Tyr 
#"integer" "integer" "integer" "integer" "integer" "integer" "integer" "integer" "integer" "integer" 
#EC 
#"factor" 

# list the levels for the class
levels(data_subset$EC)
# "EC 1" "EC 2" "EC 3" "EC 4" "EC 5" "EC 6" "EC 7"

# class distribution
percentage <- prop.table(table(data_subset$EC)) * 100
cbind(freq=table(data_subset$EC), percentage=percentage)
#freq percentage
#EC 1  551  11.991295
#EC 2 1849  40.239391
#EC 3 1677  36.496192
#EC 4  160   3.482046
#EC 5  129   2.807399
#EC 6  127   2.763874
#EC 7  102   2.219804

plot(data_subset$EC) # -> very unbalanced

### evaluate algorithms ----

control <- trainControl(method="cv", number=5)
metric <- "Accuracy"

fit.lda <- train(EC~., data=data_subset, method="lda", metric=metric, trControl=control)
fit.cart <- train(EC~., data=data_subset, method="rpart", metric=metric, trControl=control)
fit.knn <- train(EC~., data=data_subset, method="knn", metric=metric, trControl=control)
fit.svm <- train(EC~., data=data_subset, method="svmRadial", metric=metric, trControl=control)
fit.rf <- train(EC~., data=data_subset, method="rf", metric=metric, trControl=control)

results <- resamples(list(lda=fit.lda, cart=fit.cart, knn=fit.knn, svm=fit.svm, rf=fit.rf))
summary(results)


###undersampling because of class imbalance (EC 1,2 and 3)

df_EC_1 <- data_subset %>% filter(EC=="EC 1") 
df_EC_2 <- data_subset %>% filter(EC=="EC 2")
df_EC_3 <- data_subset %>% filter(EC=="EC 3")
df_EC_4 <- data_subset %>% filter(EC=="EC 4") 
df_EC_5 <- data_subset %>% filter(EC=="EC 5") 
df_EC_6 <- data_subset %>% filter(EC=="EC 6") 
df_EC_7 <- data_subset %>% filter(EC=="EC 7") 

smp_size <- 150

lines <- sample(1:nrow(df_EC_1), size = smp_size)
sample_EC_1 <- df_EC_1[lines,]

lines <- sample(1:nrow(df_EC_2), size = smp_size)
sample_EC_2 <- df_EC_2[lines,]

lines <- sample(1:nrow(df_EC_3), size = smp_size)
sample_EC_3 <- df_EC_3[lines,]

new_df <- rbind(sample_EC_1, sample_EC_2, sample_EC_3,
                df_EC_4, df_EC_5, df_EC_6, df_EC_7)


percentage <- prop.table(table(new_df$EC)) * 100
cbind(freq=table(new_df$EC), percentage=percentage)

### algorithm after undersampling ----

control <- trainControl(method="cv", number=5)
metric <- "Accuracy"

fit.lda <- train(EC~., data=new_df, method="lda", metric=metric, trControl=control)
fit.cart <- train(EC~., data=new_df, method="rpart", metric=metric, trControl=control)
fit.knn <- train(EC~., data=new_df, method="knn", metric=metric, trControl=control)
fit.svm <- train(EC~., data=new_df, method="svmRadial", metric=metric, trControl=control)
fit.rf <- train(EC~., data=new_df, method="rf", metric=metric, trControl=control)

results_undersampling <- resamples(list(lda=fit.lda, cart=fit.cart, knn=fit.knn, svm=fit.svm, rf=fit.rf))
summary(results_undersampling)


### oversampling
## TPFN und confusion matrix

### big data 
### neue daten durch die analyse
# 1. aa comp
# 2. features
# 3. beides 

### applied omics
### welche strains für die analyse?
### datenbanken abricate was ist der output
### punkt mutation, transfer -> impakt auf resistenzen
### geografisch oder haupt/nebenplasmid

## report - overleaf 
## r markdwon

