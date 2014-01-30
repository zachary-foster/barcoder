## ---- imports ----
library(png)
library(grid)
library(igraph)
library(ggplot2)
library(spider)
library(tools)
library(phyloch)
library(Hmisc)
library(knitr)
library(xtable)
library(plyr)
library(plotrix)

source("functions.R")
source("constants.R")

## ---- argument_parsing ----
#arguments <- commandArgs(TRUE)
arguments <- c("/home/local/USDA-ARS/fosterz/ITS_analysis/PR2_database/refseq_gb-191_head1000.fasta",
               "/home/local/USDA-ARS/fosterz/ITS_analysis/SSU_primers.fasta",
               "/home/local/USDA-ARS/fosterz/ITS_analysis",
               "0")
names(arguments) <- c("database", "primers", "output", "mismatch")
output_path <- file.path(arguments["output"], basename(file_path_sans_ext(arguments["primers"])))
unlink(output_path, recursive=TRUE)
dir.create(output_path)

