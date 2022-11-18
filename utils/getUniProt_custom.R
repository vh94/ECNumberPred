getUniProt_custom<-function(id){
  id <- as.character(id)
  n <- length(id)
  proteins<-empty.dump()
  for (i in 1:n) {
    tryCatch(
      {
        temp<-readFASTA(paste("https://www.uniprot.org/uniprot/", 
                                     id[i], ".fasta", sep = ""),seqonly = TRUE)[[1]]
        names(temp)<-id[i]
        proteins<-append(proteins,temp)
     },
      error=function(error_message){
        message(paste("Skip sequence: ",id[i]))
        message(error_message)
      }
      )
  }
  proteins
}
