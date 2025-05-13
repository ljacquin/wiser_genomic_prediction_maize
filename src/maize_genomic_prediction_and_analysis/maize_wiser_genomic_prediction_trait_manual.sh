#!/bin/sh
#=================================================================#
#  script for launching genomic prediction for a trait and kernel #
#=================================================================#
### Requirements
#SBATCH --partition=p01
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=90
#SBATCH --cpus-per-task=12
kernel_num_par=1
trait_num_par=4
Rscript maize_wiser_genomic_prediction_trait.R $kernel_num_par $trait_num_par
