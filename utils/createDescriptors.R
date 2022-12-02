#' Creates Descriptors from seqinr package
#' 

library(seqinr)
library(dplyr)

createDescriptors_sub <- function(sequence) {
  
  s <- seqinr::s2c(sequence[[1]])
  p <- seqinr::AAstat(s,plot = FALSE)
  
  c(unlist(p$Prop),PI=p$Pi,MW=seqinr::pmw(s))
  
}


createDescriptors <- function(AAlist) {
  AAlist %>% 
    sapply(createDescriptors_sub,simplify = TRUE) %>%
    unlist() %>% 
    as.data.frame() %>% t()
}



