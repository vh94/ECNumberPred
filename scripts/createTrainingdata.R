setwd("~/Documents/Master_CompBioCoimbra/BigData/PredECMain/")
set.seed(1)
SAMPLESIZE = 1000
## load librarys:
library(protr)
library(purrr)
library(magrittr)# base R pipe does not support "." operator yet

## source functions from utils:
invisible(sapply(paste0("utils/",list.files("utils")),source))
## create file names to aminoacid identifier files:
Enzymes<-c("Hydrolases","Isomerases","Ligases","Oxidoreductases","Transferases")
path ="ECPred/ECPred Datasets/EC_MainClass_PositiveTrain/"
filenames<-paste0(path,Enzymes,"_Positive_Train.txt")
## Get Aminoacid sequence list sample
AAsamples<-sapply(filenames, getAAsample,n=SAMPLESIZE)
names(AAsamples)<-Enzymes
## create AminoAcidContent matrices
mats <- sapply(AAsamples, createAACmatrix,simplify = FALSE)
## combine matrices to dataframe
AAC_df <- purrr::map_df(mats, ~as.data.frame(.x), .id="EC")
## compress and save:
saveRDS(AAC_df,"data/TRAINING_EC.rds")

