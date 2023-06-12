## 5-23-23
## fixing my butyrate/bile acid functions 
## so that my log10 scaled axes don't break 

library(cowplot)
library(magrittr)
library(qiime2R)
library(tidyverse)
library(viridis)
library(microshades)

## metadata 
buk_meta_fp <- 'merged_metadata1.tsv'

ko_metadata_prep <- function(metadata_fp) {
  tmpMeta <- read_tsv(metadata_fp, n_max = 2)
  mycols <- colnames(tmpMeta)
  metadata <- read_tsv(metadata_fp, skip = 2, col_names = mycols)
  names(metadata)[names(metadata) == '#SampleID'] <- 'sample'
  metadata %>% 
    select(sample, diet, study, day_post_inf) -> metadata
  return(metadata)
} 

metadata <- ko_metadata_prep(buk_meta_fp)

## taxonomy
buk_taxonomy_fp <- '~/gut_microbiome_metabolomics/lacto_sum_scaled/qiime/taxonomy.qza'

ko_tax_prep <- function(tax_fp){
  taxonomy <- read_qza(file = tax_fp)$data %>% 
    parse_taxonomy() %>% 
    as_tibble(rownames = 'taxon') %>% 
    rename_all(tolower)
  return(taxonomy)
}

taxonomy <- ko_tax_prep(buk_taxonomy_fp)

## ko contrib file 
## and overall biom table construction 
buk_contrib_fp <- 
  '~/gut_microbiome_metabolomics/lacto_sum_scaled/butyrate_bile_pred_outputs/new_picrust2_out_pipeline/lacto_meta_contrib.tsv'

ko_contrib_prep <- function(ko_contrib_fp, 
                            metadata, 
                            taxonomy){
  kometa_contrib <- read_tsv(file = ko_contrib_fp)
  names(kometa_contrib)[names(kometa_contrib) == 'function'] <- 'ko'
  kometa_contrib %>% 
    left_join(taxonomy, by = 'taxon') %>% 
    left_join(metadata, by = 'sample') -> kometa_contrib_big
  return(kometa_contrib_big)
}

ko_biom <- ko_contrib_prep(buk_contrib_fp,
                           metadata = metadata,
                           taxonomy = taxonomy)


## adding a filter step and list of diets since we want to filter out the chow diet 
## it can't be compared to the other diets since it didn't have lactococcus contamination 
## this function filters our ko biom table for the wanted taxonomic level, ko(s), and diet(s)
buk_ko <- 'K00929'
lacto_diets <- c('HF/HF', 'HF/LF', 'LF/HF', 'LF/LF')
buk_tax_level <- 'class'

ko_diet_filter <- function(kometa_contrib_big, 
                      ko_list, 
                      taxonomy_level, 
                      diet_list) { 
  kometa_contrib_big %>% 
    select(sample, ko, taxon_function_abun, study, diet, day_post_inf, 
           any_of(taxonomy_level)) %>% 
    filter(ko %in% ko_list) %>% 
    filter(diet %in% diet_list) -> filtered_biom
  return(filtered_biom)
}

buk_biom_filt <- ko_diet_filter(ko_biom,
                                buk_ko,
                                buk_tax_level,
                                lacto_diets)


## this function takes our filtered ko biom and sums the functional abundances per taxonomic level per sample 
## (so we don't have repeat points on our ggplot for the same sample)

func_abun_sum <- function(filtered_biom, 
                          taxonomy_level) {
  filtered_biom %>% 
    group_by(sample, ko, study, diet, day_post_inf, .data[[taxonomy_level]]) %>% 
    summarize(taxon_function_abun = sum(taxon_function_abun)) %>% 
    ungroup() -> filtered_biom_sum
  return(filtered_biom_sum)
}

buk_biom_sum <- func_abun_sum(buk_biom_filt,
                              buk_tax_level)

## the functions in question 
## need to be able to add sigma to our filled values so that log10 doesn't freak out 
## I'm dumb and can't figure out what sigma means besides standard deviation (which is really high for this dataset)
## so I added pi to everything. 
## this may work out horribly 

buk_threshold <- 100 
THRESHOLD = buk_threshold

buk_biom_sum %>% 
  group_by(sample, day_post_inf, ko, diet, study) %>% 
  summarize(taxon_function_abun = sum(taxon_function_abun)) %>% 
  mutate(class = 'Total') %>% 
  filter(!is.na(day_post_inf)) -> buk_stats

values <- buk_stats$taxon_function_abun
sd(values)

## these commands are creating the 'total' facet section without over-writing the 
## actual class of the bacteria 
## also putting low abundances in the 'other' category 

buk_biom_sum %>% 
  group_by(sample, day_post_inf, ko, diet, study) %>% 
  summarize(taxon_function_abun = sum(taxon_function_abun)) %>% 
  mutate(class = 'Total') %>% 
  filter(!is.na(day_post_inf)) %>% 
  spread(day_post_inf, taxon_function_abun, fill = (0 + pi)) %>% 
  gather(-sample, -ko, -diet, -study, -class,
         key = day_post_inf, value = taxon_function_abun) %>%
  mutate(day_post_inf = as.numeric(day_post_inf)) %>%
  group_by(class) %>%
  mutate(other_col = mean(taxon_function_abun),
         class = if_else(other_col < THRESHOLD, 'Other', class)) %>%
  arrange(other_col) -> buk_sample_abun



## I need to be able to rbind buk_sample_abun to the output of these commands
## so the tables need to be exactly the same so that I can put them together
## these commands are creating the main 'other' facet without over-writing the 
## actual taxonomic class 

buk_threshold <- 100 
THRESHOLD = buk_threshold

buk_biom_sum %>% 
  filter(!is.na(day_post_inf)) %>%
  spread(day_post_inf, taxon_function_abun, fill = (0 + pi)) %>% 
  gather(-sample, -ko, -diet, -study, -class,
         key = day_post_inf, value = taxon_function_abun) %>% 
  mutate(day_post_inf = as.numeric(day_post_inf)) %>% 
  group_by(class) %>% 
  mutate(other_col = mean(taxon_function_abun),
         class = if_else(other_col < THRESHOLD, 'Other', class)) %>% 
  arrange(other_col) -> buk_sum_other

## the rbinding in question 
buk_sum_other %>% 
  rbind(buk_sample_abun) -> buk_all_filt


## ggplot?
## trend line still isn't accurate but there aren't anymore half-domes ?
## pls help 
## may need to adjust trend line based off of mean? or the thing that isn't mean?

buk_all_filt %>% 
  filter(!is.na(diet)) %>% 
  ggplot(aes(x = day_post_inf, y = taxon_function_abun)) +
  scale_y_continuous(trans = 'log10') +
  geom_boxplot(aes(group = day_post_inf), outlier.shape = NA) +
  geom_smooth(se = FALSE, size = 0.5) +
  geom_vline(xintercept = -3, linetype = 'dashed', color = 'red', size = 0.2) +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'purple', size = 0.2) +
  geom_jitter(width = 0.1, height = 0, 
              alpha = 0.4) +
  scale_x_continuous(breaks = c(-15, -8, -3, 0, 3)) +
  facet_grid(class~diet, scales = 'free_x') +
  theme_bw(base_size = 14) +
  theme(legend.text = element_text(size = 8.5),
        strip.text.y = element_text(angle = 0)) +
  guides(color = guide_legend(override.aes = list(size = 0.9))) +
  xlab('Days Relative to Infection') +
  ylab('KO Counts') -> test_buk


ggsave("test_buk1.pdf",
       plot = test_buk,
       width = 15,
       height = 10, 
       path = '~/gut_microbiome_metabolomics/lacto_sum_scaled/butyrate_bile_pred_outputs')


## the functions for above (creating different facet groups)
ko_sample_abun <- function(filtered_biom_sum, threshold){
  THRESHOLD = threshold
  filtered_biom_sum %>% 
    group_by(sample, day_post_inf, ko, diet, study) %>% 
    summarize(taxon_function_abun = sum(taxon_function_abun)) %>% 
    mutate(class = 'Total') %>% 
    filter(!is.na(day_post_inf)) %>%
    spread(day_post_inf, taxon_function_abun, fill = 0) %>% 
    gather(-sample, -ko, -diet, -study, -class,
           key = day_post_inf, value = taxon_function_abun) %>% 
    mutate(day_post_inf = as.numeric(day_post_inf)) %>% 
    group_by(class) %>% 
    mutate(other_col = mean(taxon_function_abun),
           class = if_else(other_col < THRESHOLD, 'Other', class)) %>% 
    arrange(other_col) -> filtered_sample_abun
  return(filtered_sample_abun)
}

ko_sample_facet <- function(filtered_biom_sum, filtered_sample_abun, threshold){
  THRESHOLD = threshold
  filtered_biom_sum %>% 
    filter(!is.na(day_post_inf)) %>%
    spread(day_post_inf, taxon_function_abun, fill = 0) %>% 
    gather(-sample, -ko, -diet, -study, -class,
           key = day_post_inf, value = taxon_function_abun) %>% 
    mutate(day_post_inf = as.numeric(day_post_inf)) %>% 
    group_by(class) %>% 
    mutate(other_col = mean(taxon_function_abun),
           class = if_else(other_col < THRESHOLD, 'Other', class)) %>% 
    arrange(other_col) -> filtered_sum_other
  filtered_sum_other %>% 
    rbind(filtered_sample_abun) -> filtered_sum_other
  return(filtered_sum_other)
}