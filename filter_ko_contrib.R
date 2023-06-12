library(cowplot)
library(magrittr)
library(qiime2R)
library(tidyverse)

## putting the above commands into a function so that it's easier to use 
## assumes that the input file is a tsv 
## only inputs are the input file path and the output file path 

contrib_red <- function(in_fp, out_fp){
  kometa_contrib <- read_tsv(file = in_fp)
  names(kometa_contrib)[names(kometa_contrib) == 'function'] <- 'ko'
  wanted_kos <- c('K00929', 'K01034','K15873', 'K15874')
  kometa_contrib %>% 
    filter(ko %in% wanted_kos) -> min_kometa_contrib
  write_tsv(min_kometa_contrib, out_fp)
}

## testing new function 
ko_in <- '~/gut_microbiome_metabolomics/lacto_sum_scaled/butyrate_bile_pred_outputs/new_picrust2_out_pipeline/KO_metagenome_out/pred_metagenome_contrib.tsv.gz'
ko_out <- '~/gut_microbiome_metabolomics/lacto_sum_scaled/butyrate_bile_pred_outputs/TEST_meta_contrib.tsv'

contrib_red(ko_in, ko_out)


## using new function 

ko_in <- 
  '~/gut_microbiome_metabolomics/total_sum_scaled/butyrate_bile_pred_outputs/picrust2_out_pipeline_tss3/KO_metagenome_out/pred_metagenome_contrib.tsv.gz'
ko_out <- '~/gut_microbiome_metabolomics/total_sum_scaled/butyrate_bile_pred_outputs/picrust2_out_pipeline_tss3/tss3_meta_contrib.tsv'

contrib_red(ko_in, ko_out)
