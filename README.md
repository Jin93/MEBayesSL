# MEBayesSL

MEBayesSL is an R-based command line tool for implementing ME-Bayes SL, a powerful method for developing ancestry-specific polygenic risk score (PRS) that integrates information from GWAS summary statistics and external LD reference data from multiple populations (ancestry groups). ME-Bayes SL infers SNP effect sizes via a Bayesian model with an induced prior correlation structure across populations followed by an ensemble learning step with the [Super Learner](https://www.degruyter.com/document/doi/10.2202/1544-6115.1309/html).

The preprint will be put on bioRxiv soon and is available upon request. Please contact Jin Jin (Jin.Jin@Pennmedicine.upenn.edu) for details.

## Installation & Data Preparation ([Warning: tutorial incomplete])

- Download the source files from https://github.com/Jin93/MEBayesSL/. From now on we call the folder /MEBayesSL/ for simplicity.

Download the ref_bim.txt from [this link](https://www.dropbox.com/s/58uzwqewxv34wal/ref_bim.txt?dl=0) and save it under /MEBayesSL/.

- Download and decompress the LD reference data that contains: (1) raw LD reference genotype data and the (2) accompanying LD information and save them in (1) ${path_ref} and (2) ${path_LD}, respectively. 

The LD reference data contains SNP information and LD estimates by LD block for genetic variants that are in the [HapMap3](https://www.broadinstitute.org/medical-and-population-genetics/hapmap-3) plus [MEGA Chip](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5156387/) list. Note: in some scenarios, the training GWAS sample for each population consists of multiple ancestry groups, and ideally a customized LD reference dataset should be created for each population with matched ancestral composition. Code for constructing such LD reference dataset can be requested (Jin.Jin@Pennmedicine.upenn.edu).

Two options for lD reference panel are provided: 1000 Genomes Project phase 3 samples, or UK Biobank samples, both under genome build 37 (GRCh37). Each reference data contains two folders: 

(1) `./LD/`: raw LD reference genotype data, which are input files for estimating LD matrices in LDpred2 and for an intermediate step of summarizing LD information in ME-Bayes. Save the decompressed folders in ${path_ref}.

(2) `./raw/`: Precalculated LD matrices and SNP information by LD block, which are input files in ME-Bayes. Save the decompressed folders in ${path_LD}.

#### 1. LD reference data constructed based on the 1000 Genomes Project phase 3 samples (498 EUR, 659 AFR, 347 AMR, 503 EAS, 487 SAS): 

[EUR reference data](https://www.dropbox.com/s/wvxh4yqthm8m7uf/EUR.zip?dl=0) (~6.73G): `tar -zxvf EUR.tar.gz`

[AFR reference data](https://www.dropbox.com/s/iwqg65uieevfzj2/AFR.zip?dl=0) (~7.69G): `tar -zxvf AFR.tar.gz`

[AMR reference data](https://www.dropbox.com/s/mev5zyf4x6m076q/AMR.zip?dl=0) (~8.80G): `tar -zxvf AMR.tar.gz`

[EAS reference data](https://www.dropbox.com/s/o28mlovtakv5n7v/EAS.zip?dl=0) (~5.63G): `tar -zxvf EAS.tar.gz`

[SAS reference data](https://www.dropbox.com/s/idp02rgl8xv379b/SAS.zip?dl=0) (~2.60G): `tar -zxvf SAS.tar.gz`


#### 2. LD reference data constructed based on UK Biobank samples (10,000 EUR, 4,585 AFR, 687 AMR, 1,010 EAS, 5,427 SAS):

[EUR reference data]() (~13.15G): `tar -zxvf EUR.tar.gz`

[AFR reference data]() (~11.59G): `tar -zxvf AFR.tar.gz`

[AMR reference data](https://www.dropbox.com/s/2ba4tsbhz03rg83/AMR.zip?dl=0) (~4.88G): `tar -zxvf AMR.tar.gz`

[EAS reference data](https://www.dropbox.com/s/uofu788707dp4xv/EAS.zip?dl=0) (~4.27G): `tar -zxvf EAS.tar.gz`

[SAS reference data]() (~11.44G): `tar -zxvf SAS.tar.gz`


Launch R and install required libraries:

``` r
install.packages(c('optparse','bigreadr','bigsnpr','bigparallelr', 'bigmemory','stringr','caret','Rcpp', 'RcppArmadillo','RcppTN','inline','doMC','foreach','doParallel','data.table','readr','MASS','reshape','parallel',
'devtools','genio','dplyr','pryr','Matrix','lavaan','xtable','SuperLearner'))
```

Note: there are several command lines that need to be customized by users:
1. The command line "module load conda_R" in `LDpred2_jobs.R` and `MEBayes_jobs.R` may need to be modified.
2. The command lines on lines 121 - 122 in `LDpred2_jobs.R` and lines 148 - 149 in `MEBayes_jobs.R`: "qsub -cwd -l mem_free=23G,h_vmem=23G,h_fsize=100g", may need to be modified.


## MEBayesSL Manual

MEBayesSL workflow: 
<p align="center">
<img
  src="/img/MEBayesSL_Workflow.png"
  title="ME-Bayes SL Workflow"
  width=70% 
  height=70%>
</p>

#### [Step 0] 

Obtain tuned causal SNP proportion  ($p_k, k=1,2,\ldots,K$) and heritability ($h^2_k, k=1,2,\ldots,K$) for each training ancestry group from LDpred2. These parameters will be used to specify the prior causal SNP proportions and heritability parameters in ME-Bayes.

LDpred2_jobs.R: submit LDpred2 jobs by chromosome.
```r
LDpred2_jobs.R --PATH_package --PATH_data --PATH_LDref --PATH_out --FILE_sst --pop --chrom --bfile_tuning --NCORES
```

LDpred2_tuning.R: obtain tuned LDpred2 parameters.
```r
LDpred2_tuning.R --PATH_package --PATH_out --PATH_plink --FILE_sst --pop --chrom 1-22 --bfile_tuning --pheno_tuning --bfile_testing --pheno_testing --testing --NCORES
```

#### [Step 1] 

ME-Bayes: a Bayesian model that jointly models the GWAS summary data across all training populations to obtain a total of $L \times K$ PRS models under L different tuning parameter settings for $Pr⁡(δ_{1j},…,δ_{Kj})$ (functions of $p_k$s) and $\rho_{k_1,k_2}$s across all K training populations.

```r
MEBayes_jobs.R --PATH_package --PATH_data --PATH_LDref --PATH_out --FILE_sst --pop --LDpred2_params --chrom --bfile_tuning --NCORES
```

#### [Step 2] 

For the target population, apply the super learning (SL) algorithm (default base learners: elastic net regression, ridge regression, and linear regression) to train an “optimal” linear combination of the ($L \time K$) PRS models, which we call the ME-Bayes SL PRS model, based on the tuning set of the target population. Optional: the prediction performance of the final ME-Bayes SL PRS model can be reported on an independent testing set, if the testing set is provided as an input.

```r
MEBayesSL.R --PATH_package --PATH_out --PATH_plink --FILE_sst --pop --chrom --bfile_tuning --pheno_tuning --bfile_testing --pheno_testing --testing --NCORES
```

- PATH_package (required): path to the directory where the downloaded files (decompressed) are saved.

- PATH_data (required): path to the directory where the training data by ancestry group are saved.

- PATH_LDref (required): path to the directory where the LD reference data by ancestry group and chromosome are saved.

- PATH_out (required): path to the output directory where the results are saved.

- PATH_plink (required): path to plink2.

- FILE_sst (required): paths followed by file names of the population-specific GWAS summary statistics, separated by comma. Required columns: chr, rsid, pos, a0, a1, beta, beta_se, n_eff.

- pop (required: populations of the GWAS samples, separated by comma.

- chrom (required): the chromosome on which the model is fitted, input in the format of 1-22 or 1,2,3. Default: 1-22

- p: candidate values for tuning parameter p (causal SNP proportion). Default: 1e-04,0.00018,0.00032,0.00056,0.001,0.0018,0.0032,0.0056,0.01,0.018,0.032,0.056,0.1,0.18,0.32,0.56,1.

- H2: candidate values for tuning parameter H2 (heritability = H2 * h2_est from LDSC). Default: 0.7, 1, 1.4.

- sparse: whether to consider a sparse model: 0, 1, or 0,1. Default: 0.

- bfile_tuning (required): path to PLINK binary input file prefix (excluding ".bed"/".bim"/".fam"") for tuning, save by chromosome.

- pheno_tuning (optional): path to phenotype file (PLINK format) for tuning.

- covar_tuning (optional): path to quantitative covariates (PLINK format) for tuning.

- testing (required): whether to perform testing in seperate dataset. Default: F.

- bfile_testing (optional): path to PLINK binary input file prefix (.bed/.bim/.fam) for testing, save by chromosome.

- pheno_testing (optional): path to phenotype file (PLINK format) for testing.

- covar_testing (optional): path to quantitative covariates (PLINK format) for testing.

- verbose: how much chatter to print: 0=nothing; 1=minimal; 2=all. Default: 1.

- cleanup: cleanup temporary files or not. Default: T.

- NCORES: how many cores to use. (Default: 17 for LDpred2_jobs.R, 5 for MEBayes_jobs.R, and 1 for LDpred2_tuning.R and MEBayesSL.R)

- LDpred2_params (required): path to the directory where the tuned LDpred2 parameters (population-specific causal SNP proportions, heritability and whether or not a sparse model is used) are saved, separated by comma.

- cors_additional (optional): additional candidate values for tuning parameter: genetic correlation across ancestry groups, example: 3 groups with label 1,2,3, want to add two additional settings: cor_setting1(1,2),cor_setting1(1,3),cor_setting1(2,3);cor_setting2(1,2),cor_setting2(1,3),cor_setting2(2,3).

- ps_additional (optional): typically not necessary. Additional candidate values for tuning parameter: ancestry-specific causal SNP proportions, example: 3 groups with label 1,2,3, want to add two additional settings: p1_setting1,p2_setting1,p3_setting1;p1_setting2,p2_setting2,p3_setting2.

- SL_library (optional): the base learners implemented in SuperLearner, separated by comma. Default: SL.glmnet,SL.ridge,SL.lm.

- linear_score (optional): whether the trained linear models will be saved. If not, only the Super Learner model will be saved. Note: some models in SL_library are non-linear. In this case, linear score file cannot be generated.

- target_pop (required): Target population (used to save output).



## Example
Download [example data](https://www.dropbox.com/s/xxw3t17k66il3k5/example.tar.gz?dl=0), decompress it by `tar -zxvf example.tar.gz` and save the files under the directory ${path_example}. Download the 1000 Genomes reference data and save the decompressed files in ${path_LDref}. Create a new folder `path_out` (e.g., in this example, `/dcs04/nilanjan/data/jjin/mebayessl/test`) to save the output. Run the example code below with your own data directories and check if the results are consistent with the results here: [example results]().




``` r
package='/dcs04/nilanjan/data/jjin/MEBayesSL'
path_data='/dcs04/nilanjan/data/jjin/example'
path_LDref='/dcs04/nilanjan/data/jjin/LD_1kg'
path_out='/dcs04/nilanjan/data/jjin/mebayessl/test'
path_plink='/dcl01/chatterj/data/jin/software/plink2'
target_pop='AFR'
```

### Step 1: Run LDpred2 by chromosome (by submitting 22 jobs simultaneously, each for one chromosome). In each job, the algorithm will run under different tuning parameter settings in parallel.

``` r
Rscript ${package}/R/LDpred2_jobs.R \
--PATH_package ${package} \
--PATH_data ${path_data} \
--PATH_LDref ${path_LDref} \
--PATH_out ${path_out} \
--FILE_sst ${path_data}/summdata/EUR.txt,${path_data}/summdata/AFR.txt \
--pop EUR,AFR \
--chrom 1-22 \
--bfile_tuning ${path_data}/sample_data/EUR/tuning_geno,${path_data}/sample_data/AFR/tuning_geno \
--NCORES 17

```



### Step 2: Wait until all LDpred2 jobs are completed. The next step is to obtain tuned LDpred2 parameters. Tuned LDpred2 effect size estimates and the optimal tuning parameters are saved in ${path_out}.

``` r
Rscript ${package}/R/LDpred2_tuning.R \
--PATH_package ${package} \
--PATH_out ${path_out} \
--PATH_plink ${path_plink} \
--FILE_sst ${path_data}/summdata/EUR.txt,${path_data}/summdata/AFR.txt \
--pop EUR,AFR \
--chrom 1-22 \
--bfile_tuning ${path_data}/sample_data/EUR/tuning_geno,${path_data}/sample_data/AFR/tuning_geno \
--pheno_tuning ${path_data}/sample_data/EUR/pheno.fam,${path_data}/sample_data/AFR/pheno.fam \
--bfile_testing ${path_data}/sample_data/EUR/testing_geno,${path_data}/sample_data/AFR/testing_geno \
--pheno_testing ${path_data}/sample_data/EUR/pheno.fam,${path_data}/sample_data/AFR/pheno.fam \
--testing TRUE \
--NCORES 1

```

### Step 3: Run ME-Bayes by chromosome (by submitting 22 jobs simultaneously, each for one chromosome). In each job, the algorithm will run under different tuning parameter settings in parallel.

``` r
Rscript ${package}/R/MEBayes_jobs.R \
--PATH_package ${package} \
--PATH_data ${path_data} \
--PATH_LDref ${path_LDref} \
--PATH_out ${path_out} \
--FILE_sst ${path_data}/summdata/EUR.txt,${path_data}/summdata/AFR.txt \
--pop EUR,AFR \
--LDpred2_params ${path_out}/EUR/optim_params.txt,${path_out}/AFR/optim_params.txt \
--chrom 1-22 \
--bfile_tuning ${path_data}/sample_data/EUR/tuning_geno,${path_data}/sample_data/AFR/tuning_geno \
--NCORES 5

```

### Step 4: Combine PRS models generated under different parameter settings with a Super Learner (SL) algorithm to obtain the final ensembled ME-Bayes SL PRS model. Here, with the testing dataset provided, the prediction $R^2$ of the final ME-Bayes SL PRS model is reported on the testing set.

``` r
Rscript ${package}/R/MEBayesSL.R \
--PATH_package ${package} \
--PATH_out ${path_out} \
--PATH_plink ${path_plink} \
--FILE_sst ${path_data}/summdata/EUR.txt,${path_data}/summdata/AFR.txt \
--pop EUR,AFR \
--chrom 1-22 \
--bfile_tuning ${path_data}/sample_data/EUR/tuning_geno,${path_data}/sample_data/AFR/tuning_geno \
--pheno_tuning ${path_data}/sample_data/EUR/pheno.fam,${path_data}/sample_data/AFR/pheno.fam \
--bfile_testing ${path_data}/sample_data/EUR/testing_geno,${path_data}/sample_data/AFR/testing_geno \
--pheno_testing ${path_data}/sample_data/EUR/pheno.fam,${path_data}/sample_data/AFR/pheno.fam \
--testing TRUE \--NCORES 1

```
