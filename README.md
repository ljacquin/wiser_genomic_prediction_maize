[<img src="img/maize.jpg" width="600"/>]()

# genomic prediction for WISER, LS-means and BLUP phenotypes associated to maize traits

### 🎯 Objective

This repository contains R scripts designed for reproducible data analysis and results, aligned with the FAIR principles. The scripts perform data reformatting and phenotypic estimation using WISER, LS-means, and BLUP. The BLUP specifically integrate principal component coordinates of genotypes, derived from genomic data, as fixed effects to account for population structure.

### 💻 Instructions

Download the ```wiser_genomic_prediction_maize``` repository in the current user's directory on a computing cluster or personal computer using one of the following commands :

  *  ```git clone git@github.com:ljacquin/wiser_genomic_prediction_maize.git``` <p> </p>
    or
  * ```git clone https://github.com/ljacquin/wiser_genomic_prediction_maize.git``` 
  <p> </p>
  
  ⚠️ Make sure``` git``` is installed beforehand; if not, install it with ```sudo apt install git```.
  <p> </p>

* Given that ```R ≥ 4.1.2``` is already installed, within the ```wiser_genomic_prediction_maize``` folder use the following command to install and test ```wiser_genomic_prediction_maize``` required ```R``` libraries : 

  * ```R -q --vanilla < src/requirements.R```
  * ```R -q --vanilla < src/test_requirements.R```
  <p> </p>
  
* The ```R``` scripts ```0_maize_data_reformatting.R```, ```1_maize_spat_hetero_correct_per_env_trait.R``` and ```2_maize_adjusted_blups_lsmeans_phenotypes``` in the ```src/maize_data_treatment_and_analysis/``` folder perform 0) data reformatting, 1) spatial heterogeneity correction based on rows and columns within each environment (defined as site, year, management type and block), and 2) phenotype estimation using LS-means and BLUP, respectively. The BLUP incorporate principal component coordinates of genotypes, derived from genomic data, as fixed effects to account for population structure. These scripts must be executed sequentially from 0) to 2) as the inputs of the next script are the output of the previous one. Note that ```data/phenotype_data/phenotype_data.zip``` must be decompressed before executing these scripts. 

* The ```R``` script ```src/maize_genomic_prediction_and_analysis/maize_wiser_genomic_prediction_trait.R``` performs, for each trait, the genomic prediction tasks and analyses for the phenotypes estimated using WISER, LS-means, and BLUP. Note that this script also computes WISER's phenotypic estimates prior to the genomic prediction tasks.

* For genomic prediction tasks and analyses, execute the following commands to make scripts and programs executable :

  *  ```chmod u+rwx src/maize_genomic_prediction_and_analysis/*.sh```
  <p> </p>

* Finally, execute one of the following commands for executing the genomic prediction tasks and analyses :

  *  ```sbatch src/execute_maize_wiser_genomic_prediction_all_traits.sh```<p> </p>
    or
  * ```./src/execute_maize_wiser_genomic_prediction_all_traits.sh``` (i.e., interactive execution)
  <p> </p>

⚠️ The tasks and analyses performed by the ```R``` scripts in the ```wiser_genomic_prediction_maize``` repository can be run in either ```Unix/Linux``` or ```Windows``` environments, as long as ```R``` and the necessary libraries are installed. For local computations in ```RStudio```, ensure that the ```computation_mode``` variable is set to "local" in the ```R``` scripts located in ```src/```.

