#' Function to create Aminacid contents matrix from Aminoacis sequence list.
#' 
#' @param AAlist List with aminoacid sequences
#' @return matrix with aminoacid contents, columns are one letter Aminoacids
#' 
createAACmatrix <- function(AAlist) {
  AAlist %>% 
    sapply(protr::extractAAC ) %>% # extract AAC
    as.data.frame() %>% 
    t()
}