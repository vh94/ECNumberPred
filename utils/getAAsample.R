#' Function to query amioacid list from Uniprot
#' 
#' This function querys aminoacids from Uniprot by calling protr routines.
#' It randomly samples n identifiers, querys them to uniprot recieves the 
#' Aminoacid seqences, removes Gaps in the sequences, and checks for problems
#' 
#' @param path path to file containing aminoacid identifiers
#' @param n Samplesize 
#' @return A named list with Aminoacid sequences
getAAsample <- function(filename,n) {
  IDs <- scan(filename,skip = 3,character(),quote="") #read file
  stopifnot("Samplesize too large"=n<=length(IDs)) 
  message(paste("Sample",n,"out of",filename))
  Samp_IDs<-IDs %>% sample(n)  # Sample IDs
  remove(IDs) # clear cache
  # Query uniprot and clean:
  message("Query Data...")
  Samp_IDs %>%
    ### <b>
    getUniProt_custom() %>% # query uniprot
    ### </b>
    lapply(.,protr::removeGaps) %>%  # remove gaps in seq 
    subset(unlist(lapply(., protr::protcheck)))  #check sequences
}
