library(BEDMatrix)
library(RcppML)
library(extraDistr)
library(tidyverse)
library(ggpubr)
library(GGally)

# g is a matrix (n x s , n markers by s strains)
# t1 is a vector alternate allele counts  (n x 1 , n markers)

setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# load in genotypes
mats=list()
for(chr in c('I','II','III', 'IV', 'V', 'X')){
  print(chr)
  g=as.matrix(BEDMatrix(paste0('data/CeNDR20210121_Plink/', chr, '.bed') ))
  g=t(g)
  # g[is.na(g)]=0
  #to simplify for now
  g[g==2]=1
  mats[[chr]]=g
}

g=rbind(mats[[1]], mats[[2]], mats[[3]], mats[[4]], mats[[5]], mats[[6]])

rm(mats)

# clean up strain names
isotype_lookup <- data.table::fread("data/strain_isotype_lookup.tsv") 

strains <- data.table::fread("data/strains_used.txt", col.names = c("strain"), header = F,sep = " ") %>%
  dplyr::left_join(., isotype_lookup, by = "strain") %>%
  dplyr::select(strain) %>%
  dplyr::distinct()

f_strains <- paste(strains$strain,strains$strain,sep = "_") %>% sort()

# load bulk sample ALT calls
alt_dir <- "data/aser/"

for(alt_files in grep("table",list.files(alt_dir), value = T)){
  
  sample = gsub(".table.gz","",alt_files)
  print(sample)
  afdf <- data.table::fread(glue::glue("{alt_dir}{alt_files}"))%>%
    # dplyr::filter(improperPairs==0) %>%
    dplyr::select(chrom = contig, pos = position, ref = refAllele, alt = altAllele, alt_ct = altCount, dp = totalCount)%>%
    dplyr::mutate(af = alt_ct/dp,
                  sample = sample,
                  marker = paste0(chrom,":",pos,"_",alt)) 
  
  if(!exists("alt_df_bias")){
    alt_df_bias <- afdf
  } else {
    alt_df_bias <- dplyr::bind_rows(alt_df_bias, afdf)
  }
}

# number of samples in experiment
n_samples <- length(unique(alt_df_bias$sample))

# find low coverage variants across all samples
dp_cut_df <- alt_df_bias %>%
  dplyr::select(marker, sample, dp, ref, alt) %>%
  tidyr::spread(sample, dp)

dp_per_site <- apply(dp_cut_df[,4:ncol(dp_cut_df)], 1, function(x){sum(x, na.rm = T)})

dp_cut_df$sum_dp = dp_per_site

# ggplot(dp_cut_df) +
#   aes(x = sum_dp/n_samples)+
#   geom_histogram()+
#   xlim(0,200)

good_dp_marker <- dp_cut_df %>% # 2919290
  dplyr::filter(sum_dp/n_samples > 10) %>% # 2898656
  dplyr::filter(sum_dp/n_samples < 150) %>% # 2897148
  na.omit() %>% # 2896940
  dplyr::pull(marker)

# filter sites to contain those with good depth representation
pr_alt_df_bias <- alt_df_bias %>%
  dplyr::filter(marker %in% good_dp_marker)

uniq_mrkr <- unique(pr_alt_df_bias$marker)

# filter genotype matrix to contain the same markers as the ASER 
g_pruned <- g[row.names(g) %in% uniq_mrkr,colnames(g) %in% f_strains]

# g_pruned is the baseline for all other genotype processing
# pr_alt_df_bias is the baseline for all other allele count processing

save(g_pruned, pr_alt_df_bias, file="temp/test_gPruned_and_prASER.Rdata")

# # # # # # generate bootstrapping input for original genotype matrix
# 
# t1 <- pr_alt_df_bias %>%
#   dplyr::select(marker, sample, alt_ct) %>%
#   tidyr::spread(sample, alt_ct) 
# 
# t1_sorted <- t1 %>%
#   tidyr::separate(marker, into = c("chrom","pos","alt"), convert = T, remove = F) %>%
#   dplyr::arrange(chrom, pos) 
# 
# t1_sorted_pr <- t1_sorted %>%
#   dplyr::filter(marker %in% row.names(g_pruned))
# 
# identical(t1_sorted_pr$marker, row.names(g_pruned))
# 
# t1m <- dplyr::select(t1_sorted_pr, -marker, -chrom, -pos, -alt) %>%
#   as.matrix()
# 
# save(g_pruned, t1m, file="cluster_data/20220908_BulkL1_Bootstrap_Input_originalGenotypes.RData")
