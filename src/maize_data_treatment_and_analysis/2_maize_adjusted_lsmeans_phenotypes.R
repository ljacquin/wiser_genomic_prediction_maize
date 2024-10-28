# script meant to compute lsmeans for all traits
# note: text is formatted from Addins using Style active file from styler package

# clear memory and source libraries
rm(list = ls())
library(reticulate)
library(devtools)
if ("maize_env" %in% conda_list()$name) {
  use_condaenv("maize_env")
}
library(tidyverse)
library(tidyr)
library(data.table)
library(lubridate)
library(plotly)
library(htmlwidgets)
library(emmeans)
library(SpATS)
library(stringr)
library(lme4)
library(anytime)
library(foreach)
library(parallel)
library(doParallel)

# define computation mode, i.e. "local" or "cluster"
computation_mode <- "cluster"

# if comutations are local in rstudio, detect and set script path
# automatically using rstudioapi
if (identical(computation_mode, "local")) {
  library(rstudioapi)
  setwd(dirname(getActiveDocumentContext()$path))
}

# source functions
source("../functions.R")

# set options to increase memory and suppress warnings
options(expressions = 1e5)
emm_options(rg.limit = 5e5)
options(warn = -1)

# set paths
pheno_dir_path_ <- "../../data/phenotype_data/"
pheno_file_path_ <- paste0(
  pheno_dir_path_,
  "phenotype_data.csv"
)
spats_adj_pheno_path <- paste0(
  pheno_dir_path_,
  "spats_per_env_adjusted_phenotypes/"
)

# threshold for removing columns with too much na
col_na_thresh_ <- 0.3

# set vars to keep for lsmeans computation
vars_to_keep_ <- c("Envir", "Genotype")

# get file names for spats adjusted phenotypes and replace pattern
# "_spats_adjusted_.*" with "" to get trait names
files_names_spats_adj_pheno <- list.files(spats_adj_pheno_path)
trait_names_ <- str_replace_all(files_names_spats_adj_pheno,
  "_spats_adjusted_.*",
  replacement = ""
)
# initialize list for lsmeans associated to all traits
list_ls_means_adj_pheno_per_geno <- vector("list", length(trait_names_))
names(list_ls_means_adj_pheno_per_geno) <- trait_names_

# compute lsmeans across all traits
for (file_ in files_names_spats_adj_pheno) {
  print(paste0("computation for file : ", file_))

  df_ <- as.data.frame(fread(paste0(spats_adj_pheno_path, file_)))

  df_ <- df_[, c(vars_to_keep_, colnames(df_)[str_detect(
    colnames(df_),
    "spats_adj_pheno"
  )])]

  Y <- colnames(df_)[str_detect(colnames(df_), "spats_adj_pheno")]

  if (length(unique(df_$Envir)) > 1) {
    # compute adjusted ls-means for genotypes across environments
    lm_model <- lm(formula(paste0(Y, "~ Genotype + Envir")), data = df_)

    ls_means <- as.data.frame(
      lsmeans(lm_model, ~Genotype)
    )[, c("Genotype", "lsmean")]

    colnames(ls_means)[match("lsmean", colnames(ls_means))] <-
      paste0(
        str_replace_all(file_, "_spats_adjusted_.*", replacement = ""),
        "_lsmean"
      )

    list_ls_means_adj_pheno_per_geno[[
      str_replace_all(file_, "_spats_adjusted_.*",
        replacement = ""
      )
    ]] <- ls_means
  } else {
    # compute adjusted ls-means for genotypes for unique environment
    lm_model <- lm(formula(paste0(Y, "~ Genotype")), data = df_)
    ls_means <- as.data.frame(
      lsmeans(lm_model, ~Genotype)
    )[, c("Genotype", "lsmean")]

    colnames(ls_means)[match("lsmean", colnames(ls_means))] <-
      paste0(
        str_replace_all(file_, "_spats_adjusted_.*", replacement = ""),
        "_lsmean"
      )

    list_ls_means_adj_pheno_per_geno[[
      str_replace_all(file_, "_spats_adjusted_.*",
        replacement = ""
      )
    ]] <- ls_means
  }
}

#  merge list of ls_means into a single data frame for genotypes
pheno_df <- Reduce(
  function(x, y) {
    merge(x, y, by = "Genotype", all = T)
  },
  list_ls_means_adj_pheno_per_geno
)

# convert merge object to data.frame
pheno_df <- as.data.frame(pheno_df)
na_count <- colSums(is.na(pheno_df))
idx_col_to_drop <- which(na_count / nrow(pheno_df) > col_na_thresh_)
if (length(idx_col_to_drop) > 0) {
  pheno_df <- pheno_df[, -idx_col_to_drop]
}
colnames(pheno_df)[str_detect(colnames(pheno_df), "_lsmean")] <-
  str_replace_all(
    colnames(pheno_df)[str_detect(colnames(pheno_df), "_lsmean")],
    pattern = "_lsmean", replacement = ""
  )
fwrite(pheno_df,
  file = paste0(
    pheno_dir_path_,
    "adjusted_ls_mean_phenotypes.csv"
  )
)
