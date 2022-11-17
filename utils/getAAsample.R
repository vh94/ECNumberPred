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
  # Read File:
  IDs <- scan(filename,skip = 3,character(),quote="") #read file
  stopifnot("Samplesize too large"=n<=length(IDs))
  # Sample IDs
  Samp_IDs<-IDs %>% sample(n) # sample
  remove(IDs) # clear cache
  # Query uniprot and clean:
  Samp_IDs %>%
    protr::getUniProt() %>% # query uniprot
    protr::removeGaps() %>%  # remove gaps in seq
    setNames(nm=Samp_IDs) %>% # set names
    subset(unlist(lapply(., protr::protcheck))) %>% # remove problematic sequences
    list()

}