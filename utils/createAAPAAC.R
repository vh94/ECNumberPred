#' Creates Descriptors from seqinr package
#' 

library(seqinr)
library(dplyr)

createDescriptors_APAAC <- function(AAlist) {
  AAlist %>% 
    sapply(protr::extractAPAAC,simplify = TRUE) %>%
    as.data.frame() %>% 
    t()
}




