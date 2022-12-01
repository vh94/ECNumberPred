#Create evaluation Data Homo Sapiens:

## Download Data from Uniprot:
for( i in 1:7){
    url <- paste0("https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28%28taxonomy_id%3A9606%29%20AND%20%28ec%3A",i,"%29%29%20AND%20%28reviewed%3Atrue%29")
    dest <- paste0("data/raw/uniprot/EC_",i,"_HomoSapiens")
    download.file(url,method = 'wget',destfile = dest)
}

## Read files and put into list  EC_HSapiens_list
EC_HSapiens_list<- empty.dump()
for(i in 1:7){
  EC_HSapiens_list[[i]]<-seqinr::read.fasta(
    paste0("data/raw/uniprot/EC_",i,"_HomoSapiens"),
    seqtype = "AA",
    seqonly = FALSE)
}

## 