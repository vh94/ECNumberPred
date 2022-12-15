#Create evaluation Data Homo Sapiens:
library(protr)
library(purr)
## source functions from utils:
invisible(sapply(paste0("utils/",list.files("utils")),source))

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

Human_prots<-empty.dump()  ## Empty list
for( i in EC_codes){
  url <- paste0("https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28%28taxonomy_id%3A9606%29%20AND%20%28ec%3A",i,"%29%29%20AND%20%28reviewed%3Atrue%29")
  Human_prots[[names(EC_codes)[i]]]<-  protr::readFASTA(url,seqonly = TRUE) %>% 
    lapply(.,protr::removeGaps) %>%  # remove gaps in seq 
    subset(lapply(.,protr::protcheck)==TRUE)
}
saveRDS(Human_prots,"data/HSapiens_proteom.rds")


Human_prots<-readRDS("data/HSapiens_proteom.rds")
## create aminoacid content matrices For Human Proteom
Human_AAC     <- sapply(Human_prots, createAACmatrix,simplify = FALSE)
Human_AAC_df  <- purrr::map_df(Human_AAC, ~as.data.frame(.x), .id="EC")
saveRDS(Human_AAC_df,"data/HSapiens_proteom_AAC_df.rds")

## create basic Descriptors (seqinr):
Human_AA_DESC <- sapply(Human_prots,createDescriptors,simplify = FALSE)
Human_AA_DESC_df <- purrr::map_df(Human_AA_DESC, ~as.data.frame(.x), .id="EC")
saveRDS(Human_AA_DESC_df,"data/HSapiens_proteom_AA_DESC_df.rds")

## Create APAAC Descriptors for Human Proteom:
Human_APAAC   <- sapply(Human_prots, createDescriptors_APAAC,simplify = FALSE)
Human_APAAC_df <- purrr::map_df(Human_APAAC, ~as.data.frame(.x), .id="EC")
saveRDS(Human_APAAC_df,"data/HSapiens_proteom_APAAC_df.rds")









