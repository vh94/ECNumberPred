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


# testlist<-list(a=list(bla=AAsamples[[1]][2],blub=AAsamples[[3]][2]),b=list(bub=AAsamples[[1]][4],
#                foo=AAsamples[[4]][2]))
#  out<-sapply(testlist, createDescriptors,simplify = FALSE)
# df <-purrr::map_df(out, ~as.data.frame(.x), .id="EC")
# 
# 
# 
# out<-sapply(AAsamples, createDescriptors)
# purrr::map_df(out, ~as.data.frame(.x), .id="EC")


