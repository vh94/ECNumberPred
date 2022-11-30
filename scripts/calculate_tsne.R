library(Rtsne)
library(ggplot2)
library(protr)
library(caret)
library(dplyr)
library(MASS)
APAAC<-readRDS("data/TRAINING_AA_DESC_APAAC.rds")
APAAC$EC<-as.factor(APAAC$EC)
## identifing near zero variablitiy
# nearZeroVar(APAAC,freqCut = 95/1 ) # no nzv 
# ## checking for correlations:
# descrCOR<-cor(APAAC[-1])
# ## find higly correlated 
# highlyCorDescr <- findCorrelation(descrCOR, cutoff = .75)
# corrplot::corrplot(descrCOR)
# APAAC_filt<-APAAC[-1][,-highlyCorDescr]
# 
# descrCOR_filt<-cor(APAAC_filt)
# corrplot::corrplot(descrCOR_filt)


###### PCA
processing <- preProcess(APAAC, method = c("center", "scale", "pca","nzv")) # thresh  to change the PCA cut-off
transformed <- predict(processing, APAAC) 

ggplot(transformed,aes(x=PC1,y=PC2,color=APAAC$EC))+geom_point()
transformed_log<-as.data.frame(scale(log10(transformed)))
ggplot(transformed_log,aes(x=PC1,y=PC2,color=APAAC$EC))+geom_point()


ggplot(APAAC,aes(log10(Pc1.A),log10(Pc2.Hydrophobicity.1) ,color=EC))+geom_point()

#### TSNE::

### not very good and RTNSE performs scaling and PCA under the hood anyways ..
## and the default LR(200) is way to low
learning_rate<-NROW(APAAC)/12
transformed<-normalize_input(as.matrix(APAAC[-1]))

Rtsne_fit<-Rtsne(transformed,
                 pca=FALSE,
                 pca_center = FALSE,
                 pca_scale = FALSE,
                 normalize = FALSE,
                 is_distance=FALSE,
                 max_iter=500,
                 eta=learning_rate,
                 theta=0.02,
                 num_threads=0)

plot1<-Rtsne_fit$Y %>% 
  as.data.frame() %>%
  mutate(EC=APAAC$EC) %>% 
  ggplot(aes(x=V1,y=V2,color=EC))+
  geom_point(alpha=0.5)+theme_classic()
plot1

# perplexity?? perp 3*30<NROW(APAAC)-1
processing <- preProcess(APAAC, method = c("nzv","BoxCox","pca")) #PCA is perfor
transformed <- predict(processing, APAAC) 

transformend_norm<-normalize_input(as.matrix(transformed[-1]))

Rtsne_fit2<-Rtsne(transformend_norm,
                  pca=FALSE,
                  pca_center = FALSE,
                  pca_scale = FALSE,
                  normalize = FALSE,
                  is_distance=FALSE,
                  max_iter=3500,
                  eta=learning_rate,
                  theta=0,
                  num_threads=0) # use all avialable threads


plot2<-Rtsne_fit2$Y %>% 
  as.data.frame() %>%
  mutate(EC=APAAC$EC) %>% 
  ggplot(aes(x=V1,y=V2,color=EC))+
  geom_point(alpha=0.5)+theme_classic()
plot2

##nnusing distances


transformend_norm<-normalize_input(as.matrix(transformed[-1]))
dists<-dist(transformend_norm,method = "manhattan")
rm(APAAC)
rm(transformed)
rm(processing)
rm(transformend_norm)
gc()
Rtsne_dist<-Rtsne(dists, is_distance=TRUE) 





#NMDS


fit <- cmdscale(dists,eig=TRUE, k=2)

fit <- isoMDS(dists,k=2)
