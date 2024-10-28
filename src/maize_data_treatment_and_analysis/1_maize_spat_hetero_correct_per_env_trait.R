# script meant to correct for spatial heterogenity for all traits
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
options(expressions = 5e5)
options(warn = -1)

# set paths
# input and output data paths
pheno_dir_path_ <- "../../data/phenotype_data/"
pheno_file_path_ <- paste0(
  pheno_dir_path_,
  "phenotype_data.csv"
)
output_spats_file_path <- paste0(
  pheno_dir_path_,
  "spats_per_env_adjusted_phenotypes/"
)

# define function(s) and package(s) to export for parallelization
func_to_export_ <- c("fread")
pkgs_to_export_ <- c(
  "data.table", "stringr", "SpATS", "lme4",
  "lubridate", "emmeans", "plotly", "tidyr", "htmlwidgets"
)

# define selected_traits_ and vars_to_keep_ for output
selected_traits_ <- c(
  "plant.height", "tassel.height", "ear.height", "anthesis",
  "silking", "anthesis.silking.interval", "grain.number",
  "grain.yield", "grain.weight"
)
vars_to_keep_ <- c(
  "Envir", "Row", "Column", "Genotype"
)

# colors for boxplots
blue_gradient <- c("#3D9BC5", "#005AB5", "#00407A")
yellow_orange_gradient <- colorRampPalette(c("#149414", "#3b5534"))(3)
bp_colors_ <- c(blue_gradient, yellow_orange_gradient)

# define parameters for computations
min_obs_lmer_ <- 5 # cannot fit lmer if less than that.. Note  5 is pretty small
# and doesn't necessarily make sense either, its somewhat arbitrary

# get pheno_df and detect attributes, e.g. number of modalities or levels for specific variables
pheno_df_ <- as.data.frame(fread(pheno_file_path_))
management_types <- unique(pheno_df_$Management)
n_management <- length(management_types)

# define a list for singular models associated to SpATS
singular_model_list_ <<- vector("list", length(selected_traits_))
names(singular_model_list_) <- selected_traits_

# parallelize treatments for each trait_, for sequential treatment replace %dopar% by %do%
# save and return errors in singular_model_out_vect_
cl <- makeCluster(detectCores())
registerDoParallel(cl)

singular_model_out_vect_ <-
  foreach(
    trait_ = selected_traits_,
    .export = func_to_export_,
    .packages = pkgs_to_export_,
    .combine = c
  ) %dopar% {
    print(paste0("performing computation for ", trait_))

    # keep variables of interest
    df_ <- pheno_df_[, c(vars_to_keep_, trait_)]

    # get unique environments
    env_list_ <- unique(df_$Envir)

    # initialize list for spatial heterogeneity correction for each environment
    list_spats_envir_ <- vector("list", length(env_list_))
    names(list_spats_envir_) <- env_list_

    # perform a spatial heterogeneity correction for each environment for trait_
    for (env_ in env_list_)
    {
      list_spats_env_ <- spat_hetero_env_correct_trait(
        trait_,
        env_,
        df_
      )
      list_spats_envir_[[env_]] <- list_spats_env_$df_envir_

      # save environment for which an error occured
      if (list_spats_env_$message_ != "no error") {
        singular_model_list_[[trait_]] <- c(
          singular_model_list_[[trait_]],
          paste0(
            "Error for ", env_, " during SpATs: ",
            list_spats_env_$message_
          )
        )
      }
    }

    # concatenate list elements for spatial heterogeneity correction into a single df_
    df_ <- do.call(rbind, list_spats_envir_)
    df_ <- drop_na(df_)

    # sort list
    singular_model_list_[[trait_]] <-
      sort(singular_model_list_[[trait_]], decreasing = T)

    # define rename exceptions
    exception_cols <- c(
      "Genotype", "Envir",
      "Row", "Column",
      "R", "C", trait_
    )
    # rename columns excluding the exception columns
    new_names <- colnames(df_)
    new_names[!(new_names %in% exception_cols)] <- paste0(
      trait_, "_", new_names[
        !(new_names %in% exception_cols)
      ]
    )

    # replace the existing column names with the new names
    colnames(df_) <- new_names

    # write adjusted phenotype, from spatial heterogeneity correction, to long format
    fwrite(df_, paste0(
      output_spats_file_path, trait_,
      "_spats_adjusted_phenotypes_long_format.csv"
    ))

    return(paste0(
      trait_, ": ",
      singular_model_list_[[trait_]]
    ))
  }

# stop cluster
stopCluster(cl)

# write errors saved in singular_model_out_vect_
writeLines(singular_model_out_vect_, paste0(
  pheno_dir_path_,
  "envir_per_trait_with_miss_data_or_singular_spats_model.csv"
))
