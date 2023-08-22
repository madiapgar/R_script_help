## file is going to accept a biom table and a metadata file 
## outputs one thing - figure one 

library(argparse)
library(tidyverse)

parser <- ArgumentParser()

parser$add_argument("-b",
                    "--biom",
                    dest = "biom_fp",
                    help = "Filepath to biom table in X format.")
parser$add_argument("-m",
                    "--meta",
                    dest = "meta_fp",
                    help = "Filepath to metadata in X format.")
parser$add_argument("-o",
                    "--output",
                    dest = "fig1_fp",
                    help = "Filepath to location for figure one to be stored.")

args <- parser$parse_args()
print("Here I am at line 23!")
meta <- read_tsv(args$meta_fp)

write.table(meta, args$fig1_fp)

#print(paste("biom fp", args$biom_fp))
#print(paste("meta fp", args$meta_fp))
#print(paste("figure 1 fp", args$fig1_fp))