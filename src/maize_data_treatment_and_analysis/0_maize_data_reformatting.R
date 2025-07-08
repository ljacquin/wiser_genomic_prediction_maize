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
library(lme4)
library(pheatmap)

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

# set maximum number of principal components to be tested using akaike
# information criterion
max_n_comp_ <- 10

# umap parameters, most sensitive ones
random_state_umap_ <- 15
n_neighbors_umap_ <- 15
min_dist_ <- 0.1

# define traits_
traits_ <- c(
  "plant.height", "tassel.height", "ear.height",
  "anthesis", "silking", "anthesis.silking.interval",
  "grain.number", "grain.yield", "grain.weight"
)

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

# convert geno_df to numeric matrix
geno_mat <- scale(apply(geno_df, 2, as.numeric), center = T, scale = F)
k_mat <- crossprod(t(geno_mat))

# create a heatmap for genomic covariance matrix
png(
  paste0(
    output_genom_graphics_path,
    "maize_genomic_covariance_heatmap.png"
  ),
  width = 1200, height = 1200, res = 150
)
pheatmap(k_mat,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  main = "Maize dataset genomic covariance matrix heatmap",
  color = colorRampPalette(c("blue", "red"))(100),
  show_rownames = T,
  show_colnames = T
)
dev.off()

# compute pca for genomic data
geno_pca <- mixOmics::pca(apply(geno_df, 2, as.numeric), 2,
  center = T,
  scale = F,
  ncomp = max_n_comp_
)
pc_coord_df_ <- as.data.frame(geno_pca$variates)[, 1:max_n_comp_]
colnames(pc_coord_df_) <- str_replace_all(
  colnames(pc_coord_df_),
  pattern = "X.", replacement = ""
)
pc_var_names_ <- colnames(pc_coord_df_)
pc_coord_df_$Genotype <- rownames(geno_df)
geno_pca_exp_var_ <- geno_pca$explained_variance

# create pca plot
pca_fig_x_y <- plot_ly(
  data = pc_coord_df_,
  x = ~PC1, y = ~PC2,
  type = "scatter", mode = "markers"
) %>%
  layout(
    plot_bgcolor = "#e5ecf6",
    title = "PCA 2D plot for maize dataset genomic data (41722 SNP)",
    xaxis = list(title = paste0(
      names(geno_pca_exp_var_)[1], ": ",
      signif(100 * as.numeric(geno_pca_exp_var_)[1], 2), "%"
    )),
    yaxis = list(title = paste0(
      names(geno_pca_exp_var_)[2], ": ",
      signif(100 * as.numeric(geno_pca_exp_var_)[2], 2), "%"
    ))
  )

# save pca graphic
saveWidget(pca_fig_x_y, file = paste0(
  output_genom_graphics_path, "maize_pca_genomic_data.html"
))

# compute umap in 2d
geno_umap_2d_model <- umap(
  apply(geno_df, 2, as.numeric),
  n_components = 2,
  random_state = random_state_umap_,
  n_neighbors = n_neighbors_umap_,
  min_dist = min_dist_
)
geno_umap_2d <- data.frame(geno_umap_2d_model[["layout"]])

# create umap plot
umap_fig_x_y <- plot_ly(
  data = geno_umap_2d,
  x = ~X1, y = ~X2,
  type = "scatter", mode = "markers", color = "orange"
) %>%
  layout(
    plot_bgcolor = "#e5ecf6",
    title = "UMAP 2D plot for maize dataset genomic data (41722 SNP)",
    xaxis = list(title = "First component"),
    yaxis = list(title = "Second component")
  )

# save umap graphic
saveWidget(umap_fig_x_y, file = paste0(
  output_genom_graphics_path, "maize_umap_genomic_data.html"
))
