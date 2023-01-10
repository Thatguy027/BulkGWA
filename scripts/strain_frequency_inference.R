library(tidyverse)
library(RcppML)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

load("temp/prunedGenotypes_FlipCommonVariants.RData")

#################################################### 
#################################################### 
####################################################  - remove sites with any NA

gt_with_na <- apply(flipped_bootstrap_input[[1]], 1, function(x){
  sum(is.na(x))
})

a_ct <- flipped_bootstrap_input[[2]][which(gt_with_na==0),]
gt <- flipped_bootstrap_input[[1]][which(gt_with_na==0),]

fullGy=crossprod(gt,a_ct)
fullGGp=crossprod(gt)

predictions=apply(fullGy,2,function(x) {
  as.vector(RcppML::nnls(fullGGp,matrix(x), fast_nnls=T)) #fast_nnls=F, cd_maxit=10000, cd_tol=1e-10))
  
})

predictions=apply(predictions, 2, function(x) x/sum(x))

predictions_df <- data.frame(predictions) %>%
  dplyr::mutate(strain = colnames(fullGGp)) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(strain = strsplit(strain,split = "_")[[1]][1]) %>%
  tidyr::gather(sample, frq, -strain)

write.table(predictions_df, file = "temp/strain_frequencies.tsv", quote = F, row.names = F, col.names = T, sep = "\t")
