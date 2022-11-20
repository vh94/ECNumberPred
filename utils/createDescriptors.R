#' Creates Descriptors from seqinr package
#' 

library(seqinr)
library(dplyr)

createDescriptors_sub <- function(sequence) {
  
  s <- seqinr::s2c(sequence[[1]])
  
  seqinr::AAstat(s,plot = FALSE)$Prop %>%
    as.data.frame() %>%
    dplyr::mutate(PI=seqinr::computePI(s),MW=seqinr::pmw(s)) %>% 
    `rownames<-`(names(sequence))
  
  
}

createDescriptors <- function(AAlist) {
  AAlist %>% 
    sapply(createDescriptors_sub,simplify = FALSE) %>%
    as.data.frame() %>% 
    t()
}



