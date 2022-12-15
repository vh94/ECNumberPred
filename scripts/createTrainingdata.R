setwd("~/Documents/Master_CompBioCoimbra/BigData/PredECMain/")
set.seed(1)
SAMPLESIZE = 2500
## load librarys:
library(protr)
library(purrr)
library(tibble)
library(tidyr)
library(magrittr)# base R pipe does not support "." operator yet

## source functions from utils:
invisible(sapply(paste0("utils/",list.files("utils")),source))

## create file names to aminoacid identifier files:
path_train ="ECPred/ECPred Datasets/EC_MainClass_PositiveTrain/"
filenames_train<-paste0(path_train,list.files(path_train))
path_test ="ECPred/ECPred Datasets/EC_MainClass_PositiveValidation/"
filenames_test<-paste0(path_test,list.files(path_test))


## Get Aminoacid sequence list sample
AAsamples<-sapply(filenames_train, getAAsample,n=SAMPLESIZE)
## replace filepath with Enzymeclass- names
regex<-"ECPred/ECPred Datasets/EC_MainClass_PositiveTrain/\\s*(.*?)\\s*_Positive_Train.txt"
names(AAsamples)<- regmatches(names(AAsamples),regexec(regex,names(AAsamples))) %>% map_chr(2)

## create aminoacid content matrices
mats <- sapply(AAsamples, createAACmatrix,simplify = FALSE)
## combine matrices to dataframe
AAC_df <- purrr::map_df(mats, ~as.data.frame(.x), .id="EC")

## create basic Descriptors (seqinr):
AA_DESC<-sapply(AAsamples,createDescriptors,simplify = FALSE)
AA_DESC_df <- purrr::map_df(AA_DESC, ~as.data.frame(.x), .id="EC")

## create APAAC
APAAC<-sapply(AAsamples, createDescriptors_APAAC,simplify = FALSE)
AA_DESC_APAAC <- purrr::map_df(APAAC, ~as.data.frame(.x), .id="EC")

## compress and save:
saveRDS(AAsamples,"data/TRAINING_EC_SEQ.rds")
saveRDS(AAC_df,"data/TRAINING_EC.rds")
saveRDS(AA_DESC,"data/TRAINING_descriptors.rds")
saveRDS(AA_DESC_APAAC,"data/TRAINING_AA_DESC_APAAC.rds")

## read in:
AAC_d<-readRDS("data/TRAINING_EC.rds")
AAsamples<- readRDS("data/TRAINING_EC_SEQ.rds")
AA_DESC_df<-readRDS("data/TRAINING_descriptors.rds")
