# Evaluate models using Human Proteom Data:


## Read in human proteom data:
HS_prot_AAC   <- readRDS("data/HSapiens_proteom_AAC_df.rds")
HS_prot_DESC  <- readRDS("data/HSapiens_proteom_AA_DESC_df.rds")
HS_prot_APAAC <- readRDS("data/HSapiens_proteom_APAAC_df.rds")


## Merge DESC und AAC
HS_prot_merged <- merge(HS_prot_AAC, HS_prot_DESC, by =0) # by = 0 indicates by rowname
all(HS_prot_merged$EC.x == HS_prot_merged$EC.y) # check if EC is preserved persitantly

HS_prot_merged <- HS_prot_merged %>%   # remove doubled EC column
  dplyr::select(-EC.y) %>%
  dplyr::rename(EC=EC.x) %>% 
  column_to_rownames("Row.names")


# set EC as factor
HS_prot_AAC$EC    <- as.factor(HS_prot_AAC$EC)         
HS_prot_DESC$EC   <- as.factor(HS_prot_DESC$EC)
HS_prot_APAAC$EC  <- as.factor(HS_prot_APAAC$EC)
HS_prot_merged$EC <- as.factor(HS_prot_merged$EC)


## Read in trained Models:

models_AAC      <- readRDS("results/models/models_AAC.rds")
models_DESC     <- readRDS("results/models/models_features.rds")
models_merged   <- readRDS("results/models/models_merged.rds")
models_APAAC    <- readRDS("results/models/models_APAAC.rds")


All_data<-list(
  "AAC"      = HS_prot_AAC,
  "DESC"     = HS_prot_DESC,
  "AAC_DESC" = HS_prot_merged,
  "APAAC"    = HS_prot_APAAC
  )
All_data<-lapply(All_data,\(x) predict(preProcess(x, method = c("center","scale")),x))

All_models<-list(models_AAC,models_DESC,models_merged,models_APAAC)

# i want to evaluate like this:
#|    DATA        |   Models    |
#|----------------|-------------|
#|  HS_prot_AAC   | models_AAC  |
#|  HS_prot_DESC  | models_DESC |
#|  HS_prot_merged| models_merged|
#|  HS_prot_APAAC | models_APAAC |
#
# where each of the 4 models is a list with the 
# trained "lda" ,"cart" ,"knn"  &"svm" models



createConfusionMatrices<-function(data,models){
  Matrices<-empty.dump()
  
  for (d in 1:length(All_data)) {
    message(paste("Data" ,names(All_data)[[d]],":"))
    
    cmlist<-empty.dump()
    for (mod in names(All_models[[d]])) {
      message(mod)
      p     <- predict(All_models[[d]][[mod]], newdata=All_data[[d]][-1])
      cmlist[[mod]]<- caret::confusionMatrix(data = p, reference=All_data[[d]]$EC)
    }
   
    Matrices[[d]] <- cmlist
  }
  names(Matrices)<-names(All_data)
  return(Matrices)
    
}

conf_Matrices<-createConfusionMatrices(All_data,All_models)
saveRDS(conf_Matrices,"results/conf_matrices.rds")
m<-readRDS("results/conf_matrices.rds")
