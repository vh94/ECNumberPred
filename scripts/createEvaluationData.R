#Create evaluation Data Homo Sapiens:
library(protr)
library(purr)
source("utils/createAAPAAC.R")

## These are the enzyme commission numbers and their names:
EC_codes<-c(
  "Oxidoreductases"=1,
  "Transferases"=2,
  "Hydrolases"=3,
  "Lyases"=4,
  "Isomerases"=5,
  "Ligases"=6
)

## This time we don't know the uniprot IDs, only the EC we want and the taxon
## So I create a set of six URLS and six files for each EC number and download all
## HSampiens Protein fasta files from uniprot using wget:

for( i in EC_codes){
  url <- paste0("https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28%28taxonomy_id%3A9606%29%20AND%20%28ec%3A",i,"%29%29%20AND%20%28reviewed%3Atrue%29")
  dest <- paste0("data/raw/uniprot/EC_",i,"_HomoSapiens")
  download.file(url,method = 'wget',destfile = dest)
}

## Now i read the six fasta files into the session and remove problematic sequences
## and bind them into a nested list: Human_prots
readfna<-function(file){
  protr::readFASTA(file) %>% 
    lapply(.,protr::removeGaps) %>%  # remove gaps in seq 
    subset(lapply(.,protr::protcheck)==TRUE)
}

Human_prots<-empty.dump()  ## Empty list
for(i in EC_codes){
  file<-paste0("data/raw/uniprot/EC_",i,"_HomoSapiens") # same as dest
  Human_prots[[names(EC_codes)[i]]]<-readfna(file)
}
## Human_prots is a nested list with 6 sublists containing thes sequences
# > attributes(Human_prots)
# $names
# [1] "Oxidoreductases" "Transferases" "Hydrolases" "Lyases" "Isomerases" "Ligases"   

## Create a dataframe from Human_prots list:
Human_prots_df <-data.frame()
for (i in EC_codes) {
  df<-data.frame(
    "EC"  = names(EC_codes)[i],
    "Sequencee" = do.call(rbind,Human_prots[[names(EC_codes)[i]]])
  )
  Human_prots_df<-rbind(Human_prots_df,df)
}

# > glimpse(Human_prots_df)
# Rows: 4,477
# Columns: 2
# $ EC        <chr> "Oxidoreductases", "Oxidoreductases", "Oxidoreductases", ..
# $ Sequencee <chr> "MGLEALVPLAMIVAIFLLLVDLMHRHQRWAARYPPGPLPLPGLGNLLHVDFQNT ...
saveRDS(Human_prots_df,"data/HSapiens_proteom_df.rds")

## Create APAAC Descriptors for Human Proteom:
Human_APAAC <- sapply(Human_prots, createDescriptors_APAAC,simplify = FALSE)
Human_APAAC_df <- purrr::map_df(Human_APAAC, ~as.data.frame(.x), .id="EC")
saveRDS(Human_APAAC_df,"data/HSapiens_proteom_APAAC_df.rds")









