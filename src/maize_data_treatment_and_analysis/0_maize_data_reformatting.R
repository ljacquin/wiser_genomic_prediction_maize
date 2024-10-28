# script meant to reformat data for genomic prediction
# note: text is formatted from Addins using Style active file from styler package

# clear memory and source libraries
rm(list = ls())
library(devtools)
install_other_requirements <- F
if (install_other_requirements) {
  install.packages("BiocManager")
  library(BiocManager)
  BiocManager::install("snpStats")
  BiocManager::install("mixOmicsTeam/mixOmics")
  py_install("umap-learn", pip = T, pip_ignore_installed = T)
}
library(mixOmics)
library(data.table)
library(plotly)
library(ggplot2)
library(umap)
library(dplyr)
library(Matrix)
library(graphics)
library(htmlwidgets)
library(rstudioapi)
library(stringr)
library(tidyr)
library(dplyr)
library(lsmeans)

# detect and set script path automatically, and source functions
setwd(dirname(getActiveDocumentContext()$path))
source("../functions.R")

# set options to increase memory
options(expressions = 5e5)
options(warn = -1)
emm_options(rg.limit = 10e6)

# set path for genomic data and phenotype data
genom_dir_path <- "../../data/genomic_data/"
pheno_dir_path <- "../../data/phenotype_data/"

# set output result path for genomic graphics
output_genom_graphics_path <- "../../results/genomic_prediction_graphics/"

# get genomic data
geno_df <- as.data.frame(fread(paste0(
  genom_dir_path,
  "7a-Genotyping_50K_41722.csv"
)))
row_names_ <- geno_df[, 1]
rownames(geno_df) <- row_names_
geno_df <- geno_df[, -1]
dim(geno_df)

# get phenotype data
pheno_df <- as.data.frame(fread(paste0(
  pheno_dir_path,
  "2a-GrainYield_components_Plot_level.csv"
)))

# rename columns associated to genotypes and environments
colnames(pheno_df)[
  match("Variety_ID", colnames(pheno_df))
] <- "Genotype"
pheno_df$Envir <- paste0(
  pheno_df$Experiment,
  "_block_", pheno_df$block
)

# replace management types by "R" and "W"
pheno_df <- pheno_df %>%
  mutate(Management = ifelse(treatment == "rainfed", "R", "W"))

# apply function to generate row and position columns by environment
pheno_df <- generate_row_column_variables_by_environment(pheno_df)

# get common genotypes between genomic and phenotype data
pheno_df <- match_indices(pheno_df, geno_df)
geno_df <- geno_df[rownames(geno_df) %in% pheno_df$Genotype, ]
dim(geno_df)
length(unique(pheno_df$Genotype))

# write reformatted datasets
fwrite(pheno_df,
  file = paste0(pheno_dir_path, "phenotype_data.csv")
)
fwrite(as.data.frame(geno_df),
  file = paste0(genom_dir_path, "genomic_data.csv"),
  row.names = T
)
