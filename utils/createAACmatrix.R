#' Function to create Aminacid contents matrix from Aminoacids sequence list.
#' 
#' @param AAlist List with aminoacid sequences
#' @return list of matrices with AAC, columns are one letter AminAcid IDs
#' 
createAACmatrix <- function(AAlist) {
  AAlist %>% 
    ### <b>
    sapply(protr::extractAAC,simplify = TRUE ) %>% 
    ### </b>
    as.data.frame() %>% 
    t()
}