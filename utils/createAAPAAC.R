#' Function to createType II PseAAC descriptor form list of Aminoacid-sequences.
#' 
#' @param AAlist List with aminoacid sequences
#' @return list of matrices with Type II PseAAC descriptors
#' 
createDescriptors_APAAC <- function(AAlist) {
  AAlist %>% 
    ### <b>
    sapply(protr::extractAPAAC,simplify = TRUE) %>%
    ### </b>
    as.data.frame() %>% 
    t()
}




