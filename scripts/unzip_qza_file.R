## 7-27-23
## creating an R function that can read in zipped files (.qzas) and unzip them

library(tidyverse)
library(magrittr)
library(broom)
library(qiime2R)

## experimental file, file path
unweighted_pcoa  <- './data/qiime/core_outputs/unweighted_unifrac_pcoa_results.qza'

## experimental code
zipdir <- tempfile()
dir.create(zipdir)

unzip(unweighted_pcoa, exdir=zipdir)

artifact.dir <- list.files(zipdir)
full.path <- file.path(zipdir, artifact.dir)

files <- list.files(file.path(full.path, "data"))
test <- read_tsv(file.path(full.path, "data", files), skip = 9, col_names = FALSE)

temp_col <- colnames(test)
new_cols <- gsub(pattern = "X", replacement = "PC", x=temp_col)

new_test <- test <- read_tsv(file.path(full.path, "data", files), skip = 9, col_names = new_cols)
names(new_test)[1:ncol(new_test)+1] <- names(new_test)[2:ncol(new_test)]


## using the read_qza function 
actual <- read_qza(unweighted_pcoa)$data
vectors <- actual$Vectors
