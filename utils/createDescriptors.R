#' Creates Descriptors from seqinr package
#' 

library(seqinr)
library(dplyr)

createDescriptors_sub <- function(sequence) {
  
  s <- seqinr::s2c(sequence[[1]])
  
  c(seqinr::AAstat(s,plot = FALSE)$Prop,
    PI=seqinr::computePI(s),
    MW=seqinr::pmw(s))
  
}


createDescriptors <- function(AAlist) {
  AAlist %>% 
    sapply(createDescriptors_sub,simplify = TRUE) %>%
    as.data.frame() %>% 
    t()
}



