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
library(MASS)

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

# set maximum number of principal components to be tested using akaike
# information criterion
max_n_comp_ <- 10

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
genom_dir_path <- "../../data/genomic_data/"

# set vars to keep for lsmeans computation
vars_to_keep_ <- c("Envir", "Genotype")

# get file names for spats adjusted phenotypes and replace pattern
# "_spats_adjusted_.*" with "" to get trait names
files_names_spats_adj_pheno <- list.files(spats_adj_pheno_path)
trait_names_ <- str_replace_all(files_names_spats_adj_pheno,
  "_spats_adjusted_.*",
  replacement = ""
)
# initialize lists for lsmeans and blups associated to all traits
list_ls_means_adj_pheno_per_geno <- vector("list", length(trait_names_))
names(list_ls_means_adj_pheno_per_geno) <- trait_names_
blup_list_ <- vector("list", length(trait_names_))
names(blup_list_) <- trait_names_

# initialize vector for aic values
aic_ <- rep(0, max_n_comp_)

# compute principal components to be used as fixed covariates for population
# structure correction
geno_df <- as.data.frame(fread(paste0(genom_dir_path, "genomic_data.csv")))
geno_names_ <- geno_df[, 1]
geno_df <- geno_df[, -1]
geno_pca <- mixOmics::pca(apply(geno_df, 2, as.numeric), ncomp = max_n_comp_)
pc_coord_df_ <- as.data.frame(geno_pca$variates)[, 1:max_n_comp_]
pc_var_names_ <- colnames(pc_coord_df_)
pc_coord_df_$Genotype <- geno_names_

# for each trait, compute blups and lsmeans across all environments
# for each genotype
for (file_ in files_names_spats_adj_pheno) {
  print(paste0("computation for file : ", file_))

  df_ <- as.data.frame(fread(paste0(spats_adj_pheno_path, file_)))
  df_ <- df_[, c(vars_to_keep_, colnames(df_)[str_detect(
    colnames(df_),
    "spats_adj_pheno"
  )])]
  Y <- colnames(df_)[str_detect(colnames(df_), "spats_adj_pheno")]
  # for isolated genotypes with single phenotypic values (i.e. non repeated),
  # repeat their unique phenotypic values one time to make blup computation
  # possible
  idx_non_repeat_geno <- which(table(df_$Genotype) < 2)
  if (length(idx_non_repeat_geno) > 0) {
    df_ <- rbind(df_, df_[idx_non_repeat_geno, ])
  }
  # merge principal components and adjusted phenotypes based on genotype key
  df_ <- merge(df_, pc_coord_df_,
    by = "Genotype", all = TRUE
  )
  # remove any NA for trait adjusted phenotype before blup and lsmeans computation
  df_ <- df_[!is.na(df_[, Y]), ]

  if (length(unique(df_$Envir)) > 1) {
    # compute blups for genotypes using a linear mixed model (LMM) which fits
    # pc as fixed effects inorder to account for population structure

    # compute aic values in order to select number of pcs
    for (n_comp_ in 1:max_n_comp_) {
      lmer_model_ <- lmer(
        as.formula(paste0(
          Y,
          " ~ 1 + Envir + ", paste(pc_var_names_[1:n_comp_],
            collapse = " + "
          ),
          " + (1 | Genotype)"
        )),
        data = df_
      )
      aic_[n_comp_] <- AIC(lmer_model_)
    }
    n_opt_comp_aic_ <- which.min(aic_)
    print(paste0("number of pc selected: ", n_opt_comp_aic_))

    # estimate model based on selected number of pcs which minimize aic
    lmer_model_ <- lmer(
      as.formula(paste0(
        Y,
        " ~ 1 + Envir + ", paste(pc_var_names_[1:n_opt_comp_aic_],
          collapse = " + "
        ),
        " + (1 | Genotype)"
      )),
      data = df_
    )
    blup_list_[[str_replace_all(file_, "_spats_adjusted_.*",
      replacement = ""
    )]] <- data.frame(
      "Genotype" = rownames(ranef(lmer_model_)$Genotype),
      "blup" = as.numeric(unlist(ranef(lmer_model_)$Genotype))
    )

    # compute adjusted ls-means for genotypes across environments
    lm_model <- lm(formula(paste0(Y, "~ 1 + Genotype + Envir")), data = df_)
    ls_means <- as.data.frame(
      lsmeans(lm_model, ~Genotype)
    )[, c("Genotype", "lsmean")]

    # if the model coefficients are not estimable due to collinearity, use lm_()
    # to compute the coefficients with a pseudo-inverse approach. Then, use
    # lsmeans_() for calculating least squares means, which is based on group_by()
    # to replace the original lsmeans() function, as the latter requires a standard lm object.
    if (sum(is.na(lm_model$coefficients)) > 0) {
      lm_model <- lm_(formula(paste0(Y, "~ 1 + Genotype + Envir")), data = df_)
      ls_means <- as.data.frame(
        lsmeans_(lm_model, df_)
      )
    }
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
    # compute blups for genotypes using a linear mixed model (LMM) which fits
    # pc as fixed effects inorder to account for population structure

    # compute aic values in order to select number of pcs
    for (n_comp_ in 1:max_n_comp_) {
      lmer_model_ <- lmer(
        as.formula(paste0(
          Y, " ~ 1 + ", paste(pc_var_names_[1:n_comp_],
            collapse = " + "
          ),
          " + (1 | Genotype)"
        )),
        data = df_
      )
      aic_[n_comp_] <- AIC(lmer_model_)
    }
    n_opt_comp_aic_ <- which.min(aic_)
    print(paste0("number of pc selected: ", n_opt_comp_aic_))

    # estimate model based on selected number of pcs which minimize aic
    lmer_model_ <- lmer(
      as.formula(paste0(
        Y, " ~ 1 + ", paste(pc_var_names_[1:n_opt_comp_aic_],
          collapse = " + "
        ),
        " + (1 | Genotype)"
      )),
      data = df_
    )
    blup_list_[[str_replace_all(file_, "_spats_adjusted_.*",
      replacement = ""
    )]] <- data.frame(
      "Genotype" = rownames(ranef(lmer_model_)$Genotype),
      "blup" = as.numeric(unlist(ranef(lmer_model_)$Genotype))
    )

    # compute adjusted ls-means for genotypes for unique environment
    lm_model <- lm(formula(paste0(Y, "~ 1 + Genotype")), data = df_)
    ls_means <- as.data.frame(
      lsmeans(lm_model, ~Genotype)
    )[, c("Genotype", "lsmean")]

    # if the model coefficients are not estimable due to collinearity, use lm_()
    # to compute the coefficients with a pseudo-inverse approach. Then, utilize
    # lsmeans_() for calculating least squares means, which is based on group_by()
    # to replace the original lsmeans() function, as the latter requires a standard lm object.
    if (sum(is.na(lm_model$coefficients)) > 0) {
      lm_model <- lm_(formula(paste0(Y, "~ 1 + Genotype")), data = df_)
      ls_means <- as.data.frame(
        lsmeans_(lm_model, df_)
      )
    }
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

# reduce blup list
blup_df <- Reduce(
  function(x, y) {
    merge(x, y,
      by = "Genotype",
      all = T
    )
  },
  blup_list_
)
colnames(blup_df) <- c("Genotype", trait_names_)

# write blups
fwrite(blup_df, file = paste0(
  pheno_dir_path_,
  "blup_phenotypes.csv"
))

# reduce lsmeans list
lsmean_df <- Reduce(
  function(x, y) {
    merge(x, y, by = "Genotype", all = T)
  },
  list_ls_means_adj_pheno_per_geno
)

# write lsmeans
colnames(lsmean_df) <- c("Genotype", trait_names_)

fwrite(lsmean_df,
  file = paste0(
    pheno_dir_path_,
    "adjusted_ls_mean_phenotypes.csv"
  )
)
