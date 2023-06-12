## 6-12-23
## putting together functions for Qiime2 core diversity analysis plots 

## needed libraries
library(qiime2R)
library(tidyverse)
library(cowplot)
library(magrittr)
library(vegan)
library(viridis)
library(microshades)
library(phyloseq)
library(ggh4x)
library(broom)

## metadata file editor 
## assumes the file is a .tsv 
## can be personalized (lines 12-15)
metadata_fixer <- function(metadata_fp) {
  tmpMeta <- read_tsv(metadata_fp, n_max = 2)
  mycols <- colnames(tmpMeta)
  metadata <- read_tsv(metadata_fp, skip = 2, col_names = mycols)
  names(metadata)[names(metadata) == '#SampleID'] <- 'sampleid'
  metadata %>% 
    filter(!is.na(diet)) %>% 
    mutate(day_post_inf = if_else(day_post_inf == 2, 3, day_post_inf)) %>% 
    mutate(diet = as.factor(diet)) -> metadata
  return(metadata)
}

## similar function to above for Qiime2 files that need the second row taken out
## assumes the input is a .qza file 
file_fixer <- function(metadata_fp) {
  tmpMeta <- read_tsv(metadata_fp, n_max = 2)
  mycols <- colnames(tmpMeta)
  metadata <- read_tsv(metadata_fp, skip = 2, col_names = mycols)
  names(metadata)[names(metadata) == '#SampleID'] <- 'sampleid'
  return(metadata)
}

## this function will pull out the percent variations from a specified column so you can add it to your pcoa plots 
## you need to have the Proportion Explained pulled out of the weighted/unweighted unifrac pcoa results file as a string
## col_name should be "PC1, PC2, etc"
pcoa_ax_lab <- function(unifrac_var, col_name){
  uni_lab <- as.character(round(unifrac_var[col_name] * 100, 2))
  uni_lab <- paste0(col_name, ' - ', uni_lab, '%')
  return(uni_lab)
}

## unweighted/weighted unifrac pcoa result, faith's pd, and shannon entropy file prep 
## going to attempt to return multiple outputs so I can just have one function for file prep
biom_table_prep <- function(unweighted_fp,
                            weighted_fp,
                            faith_fp,
                            shannon_fp,
                            metadata_file){
  ## unweighted pcoa
  unweighted <- read_qza(unweighted_fp)$data
  unweighted_var <- unweighted$ProportionExplained
  unweighted_pcoa <- unweighted$Vectors ##used for pcoa plot
  names(unweighted_pcoa)[names(unweighted_pcoa) == 'SampleID'] <- 'sampleid'
  ## weighted pcoa
  weighted <- read_qza(weighted_fp)$data
  weighted_var <- weighted$ProportionExplained
  weighted_pcoa <- weighted$Vectors
  names(weighted_pcoa)[names(weighted_pcoa) == 'SampleID'] <- 'sampleid'
  ## faith's 
  faith <- read_tsv(faith_fp)
  names(faith)[names(faith) == '#SampleID'] <- 'sampleid'
  ## shannon 
  shannon <- read_tsv(shannon_fp)
  names(shannon)[names(shannon) == '...1'] <- 'sampleid'
  ## unweighted biom 
  unweighted_pcoa %>% 
    left_join(metadata_file, by = 'sampleid') %>% 
    left_join(faith, by = 'sampleid') %>% 
    left_join(shannon, by = 'sampleid') -> unweighted_biom
  ## weighted biom
  weighted_pcoa %>% 
    left_join(metadata_file, by = 'sampleid') %>% 
    left_join(faith, by = 'sampleid') %>% 
    left_join(shannon, by = 'sampleid') -> weighted_biom
  ## creating a list to return multiple outputs 
  my_list <- list(UnweightedVar = unweighted_var, 
                  WeightedVar = weighted_var,
                  UnweightedBiom = unweighted_biom,
                  WeightedBiom = weighted_biom)
  return(my_list)
}

## pcoa plot function
## xlab and ylab are outputs from pcoa_ax_lab function
pcoa_plot <- function(biom_file,
                      labels,
                      names_labels,
                      xlab,
                      ylab,
                      title){
  ## what you want the grid labels to be (list)
  labs <- (labels) 
  ## what grid labels currently are (list)
  names(labs) <- (names_labels) 
  biom_file %>% 
    filter(!is.na(diet)) %>% 
    ggplot(aes(x = PC1, y = PC2)) +
    geom_point(aes(fill = faith_pd), pch = 21, alpha = 0.7) +
    theme_bw(base_size = 14) +
    scale_fill_distiller(palette = 'Spectral', name = "Faith's PD") +
    facet_grid(day_post_inf~diet, 
               labeller = labeller(diet = labs)) +
    theme(legend.text = element_text(size = 8.5),
          strip.text.y = element_text(angle = 0)) +
    ggtitle(title) +
    labs(x = xlab, y = ylab) -> pcoa
  return(pcoa)
}

## faith's pd plot 
## assumes that the files is a .tsv
faith_pd_plot <- function(faith_fp,
                          metadata_file,
                          labels,
                          names_labels,
                          title){
  faith <- read_tsv(faith_fp)
  names(faith)[names(faith) == '#SampleID'] <- 'sampleid'
  metadata_file %>% 
    left_join(faith, by = 'sampleid') -> faith_pd
  ## what you want the grid labels to be (list)
  labs <- (labels) 
  ## what grid labels currently are (list)
  names(labs) <- (names_labels) 
  faith_pd %>%
    filter(!is.na(diet)) %>% 
    ggplot(aes(x = day_post_inf, y = faith_pd)) +
    geom_boxplot(aes(group = day_post_inf), outlier.shape = NA) +
    geom_jitter(width = 0.1, height = 0, alpha = 0.4) +
    scale_x_continuous(breaks = c(-15, -8, -3, 0, 3)) +
    geom_smooth(se = FALSE) +
    theme_bw(base_size = 16) +
    facet_grid(~diet, labeller = labeller(diet = labs) ) +
    ggtitle(title) +
    xlab('Days Relative to Infection') +
    ylab("Faith's PD") -> faith_plot
  return(faith_plot)
}

## shannon entropy plot
## assumes that the file is a .tsv
shannon_plot <- function(shannon_fp,
                          metadata_file,
                          labels,
                          names_labels,
                          title){
  shannon <- read_tsv(shannon_fp)
  names(shannon)[names(shannon) == '...1'] <- 'sampleid'
  metadata_file %>% 
    left_join(shannon, by = 'sampleid') -> shannon_entropy
  ## what you want the grid labels to be (list)
  labs <- (labels) 
  ## what grid labels currently are (list)
  names(labs) <- (names_labels) 
  shannon_entropy %>%
    filter(!is.na(diet)) %>% 
    ggplot(aes(x = day_post_inf, y = shannon_entropy)) +
    geom_boxplot(aes(group = day_post_inf), outlier.shape = NA) +
    geom_jitter(width = 0.1, height = 0, alpha = 0.4) +
    scale_x_continuous(breaks = c(-15, -8, -3, 0, 3)) +
    geom_smooth(se = FALSE) +
    theme_bw(base_size = 16) +
    facet_grid(~diet, labeller = labeller(diet = labs) ) +
    ggtitle(title) +
    xlab('Days Relative to Infection') +
    ylab("Shannon Entropy") -> shannon_plot
  return(shannon_plot)
}

## taxa barplot can be constructed using the microshades functions


## testing out my functions to make sure they work properly 

## metadata file prep
metadata_FP <- '~/gut_microbiome_metabolomics/total_sum_scaled/tss_R/merged_metadata1.tsv'
metadata <- metadata_fixer(metadata_FP)

## all other files (hopefully)
## to pull individual dataframes out of the result of this function, use
## core_files$(UnweightedVar, WeightedVar, UnweightedBiom, WeightedBiom)
tss_unweighted <- '~/gut_microbiome_metabolomics/total_sum_scaled/tss3/tss_core_metrics3/unweighted_unifrac_pcoa_results.qza'
tss_weighted <- '~/gut_microbiome_metabolomics/total_sum_scaled/tss3/tss_core_metrics3/weighted_unifrac_pcoa_results.qza'
tss_faith <- '~/gut_microbiome_metabolomics/total_sum_scaled/tss3/tss_core_metrics3/faith_pd/data/alpha-diversity.tsv'
tss_shannon <- '~/gut_microbiome_metabolomics/total_sum_scaled/tss3/tss_core_metrics3/shannon_entropy/data/alpha-diversity.tsv'

core_files <- biom_table_prep(tss_unweighted,
                              tss_weighted,
                              tss_faith,
                              tss_shannon,
                              metadata)
## pulling the individual dataframes out of the list to use for plots in the future 
uw_var <- core_files$UnweightedVar
w_var <- core_files$WeightedVar
unweighted_biom <- core_files$UnweightedBiom
weighted_biom <- core_files$WeightedBiom

## pcoa plot construction
## it works!!
uw_uni_xlab <- pcoa_ax_lab(uw_var, "PC1")
uw_uni_ylab <- pcoa_ax_lab(uw_var, "PC2")

diet_labs <- c('Chow', 
    'High Fat / High Fiber', 
    'High Fat / Low Fiber', 
    'Low Fat / High Fiber', 
    'Low Fat / Low Fiber')

diet_names_labels <- c('Chow', 'HF/HF', 'HF/LF', 'LF/HF', 'LF/LF')

uw_title <- 'Total Sum Scaled Unweighted UniFrac'

unweighted <- pcoa_plot(unweighted_biom,
                        diet_labs,
                        diet_names_labels,
                        uw_uni_xlab,
                        uw_uni_ylab,
                        uw_title)

## faith's pd plot construction
## it works!!
tss_faith <- '~/gut_microbiome_metabolomics/total_sum_scaled/tss3/tss_core_metrics3/faith_pd/data/alpha-diversity.tsv'
diet_labs <- c('Chow', 
               'High Fat / High Fiber', 
               'High Fat / Low Fiber', 
               'Low Fat / High Fiber', 
               'Low Fat / Low Fiber')
diet_names_labels <- c('Chow', 'HF/HF', 'HF/LF', 'LF/HF', 'LF/LF')
faith_title <- "Total Sum Scaled Faith's Phylogenic Diversity"

faith_plot <- faith_pd_plot(tss_faith,
                            metadata,
                            diet_labs,
                            diet_names_labels,
                            faith_title)

## shannon entropy plot construction
## it works!!
tss_shannon <- '~/gut_microbiome_metabolomics/total_sum_scaled/tss3/tss_core_metrics3/shannon_entropy/data/alpha-diversity.tsv'
diet_labs <- c('Chow', 
               'High Fat / High Fiber', 
               'High Fat / Low Fiber', 
               'Low Fat / High Fiber', 
               'Low Fat / Low Fiber')
diet_names_labels <- c('Chow', 'HF/HF', 'HF/LF', 'LF/HF', 'LF/LF')
shannon_title <- "Total Sum Scaled Shannon Entropy"

shannon <- shannon_plot(tss_shannon,
                        metadata,
                        diet_labs,
                        diet_names_labels,
                        shannon_title)
