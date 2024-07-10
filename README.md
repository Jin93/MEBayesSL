# MUSSEL

MUSSEL is an R-based command line tool for implementing MUSSEL, a method for developing ancestry-specific polygenic risk score (PRS) that integrates information from GWAS summary statistics and external LD reference data from multiple populations (ancestry groups). MUSSEL infers SNP effect sizes via a Bayesian model with an induced prior correlation structure across populations followed by an ensemble learning step with the Super Learner. As intermediate products, LDpred2 PRS models trained separately on GWAS data for each population will also be generated (please see [Example](#example) for details).

To use the tool, please follow the instructions in [Getting Started](#gettingstarted) to download the required code and data, then try our example code in [Example](#example). Please refer to the [paper](https://www.cell.com/cell-genomics/fulltext/S2666-979X(24)00095-8) or contact Jin Jin (Jin.Jin@Pennmedicine.upenn.edu) for details.
</br>



## Version History
- [ ] __April 10, 2024:__  Updated MUSSEL to incorporate the most recent version of [LDpred2](https://privefl.github.io/bigsnpr/articles/LDpred2.html) (June 8, 2023 version): now LDpred2_jobs.R will only submit one job to the server instead of 22 (by chromosome) in the previous version; fixed bugs with covariate adjustment; now allow for calculating AUC of PRS for binary traits.

- [ ] __Feb 3, 2023:__  The first version of MUSSEL was made available on Github.
</br>




## Getting Started

- Download or clone the Github repository by `git clone https://github.com/Jin93/MUSSEL.git`. From now on, we will refer to the folder as `/MUSSEL/` for simplicity.

- Download the `ref_bim.txt` from [this link](https://www.dropbox.com/s/58uzwqewxv34wal/ref_bim.txt?dl=0) and save it under `/MUSSEL/`.


- Download the LD reference information and save the unzipped folder in ${path_LDref}.

The (pre-computed) LD information contains SNP information and LD matrix estimates for:
  1. ~ 2.18 million [HapMap3](https://www.broadinstitute.org/medical-and-population-genetics/hapmap-3) and [MEGA ](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5156387/) SNPs.
  2. For five ancestry groups: European (EUR), African/African America (AFR), Hispanic/Latino/Admixed American (AMR), East Asian (EAS), and South Asian (SAS).
  3. Based on either the UK Biobank samples (recommended when training GWAS sample sizes are relatively large, e.g., N<sub>GWAS</sub> > 50K for at least one ancestry), or the 1000 Genomes samples (recommended when GWAS training sample sizes are relatively small, e.g., N<sub>GWAS</sub> < 50K for all ancestry groups).


Note: in some scenarios, the training GWAS sample for a population consists of multiple ancestry groups. In this case, ideally, a customized LD reference dataset should be created for this population with matched ancestral composition. The code for constructing such LD reference dataset can be accessed [here](R/Generate_LD_info_by_LDblock.R).


Please choose one of the two LD reference panels according to the training GWAS sample sizes (please see detailed description below). Each reference data contains two folders: 

`./precalLD/`: Pre-computed LD matrices and SNP information by LD block for implementing LDpred2 ([Step 0](#step-0:-run-ldpred2-by-ancestry)).

`./LD/`: Pre-computed LD matrices and SNP information by LD block, which are input files in MUSS.

`./map/`: SNP information (SNP ID, alleles) for mapping alleles between LD reference panel and GWAS summary data in LDpred2 and MUSS.

`./raw/`: raw LD reference genotype data (.bim) by chromosome, which will be used in MUSS ([Step 1](#step-1:-muss)).



### 1. UK Biobank LD reference data (currently been updated)

- 10,000 EUR, 4,585 AFR, 687 AMR, 1,010 EAS, 5,427 SAS.
- __Recommended when training GWAS sample sizes are relatively large, e.g., N<sub>GWAS</sub> > 50K for at least one ancestry groups__.


### 2. LD reference data constructed based on the 1000 Genomes Project phase 3 samples 

- 498 EUR, 659 AFR, 347 AMR, 503 EAS, 487 SAS.
- __Recommended when GWAS training sample sizes are relatively small, e.g., N<sub>GWAS</sub> < 50K for all ancestry groups__.

[EUR reference data](https://drive.google.com/file/d/1448w1cMdFuBmhnKEpamUBkGWX93WcaI-/view?usp=drive_link) (~26.33G): Google Drive File ID: `1448w1cMdFuBmhnKEpamUBkGWX93WcaI-`; decompress by `tar -zxvf EUR.tar.gz`

[AFR reference data](https://drive.google.com/file/d/12j9y0brD004HNakC7LoEot8RHFhToVAq/view?usp=drive_link) (~36.81G): Google Drive File ID: `12j9y0brD004HNakC7LoEot8RHFhToVAq`; decompress by `tar -zxvf AFR.tar.gz`

[AMR reference data](https://drive.google.com/file/d/1ItDxm8mllEAoaUyddd_wz8UTp0clNmLW/view?usp=drive_link) (~32.82G): Google Drive File ID: `1ItDxm8mllEAoaUyddd_wz8UTp0clNmLW`; decompress by `tar -zxvf AMR.tar.gz`

[EAS reference data](https://drive.google.com/file/d/1z4bztGLssmXmOiA4AHq4jzRQ5u9A-w7-/view?usp=drive_link) (~19.90G): Google Drive File ID: `1z4bztGLssmXmOiA4AHq4jzRQ5u9A-w7-`; decompress by `tar -zxvf EAS.tar.gz`

[SAS reference data](https://drive.google.com/file/d/1PCiy_rMDGxPpBZpQzRSVNocc9ZwoHt0-/view?usp=drive_link) (~21.61G): Google Drive File ID: `1PCiy_rMDGxPpBZpQzRSVNocc9ZwoHt0-`; decompress by `tar -zxvf SAS.tar.gz`


- Install [PLINK1.9](https://www.cog-genomics.org/plink/) and [PLINK2](https://www.cog-genomics.org/plink/2.0/).

- Launch R and install required libraries:

``` r
install.packages(c('RISCA','optparse','bigreadr','bigsnpr','bigparallelr', 'bigmemory','stringr','caret','Rcpp', 'RcppArmadillo','RcppTN','inline','doMC','foreach','doParallel','data.table','readr','MASS','reshape','parallel',
'devtools','genio','dplyr','pryr','Matrix','lavaan','xtable','SuperLearner'))
```

- Prepare input data files:

Please use the `example_data` in [Example](#example) as a reference to prepare  the input data files. Please note that the SNPs are matched across the training GWAS data, LD reference data, and validation data by RSID.

An example of the GWAS summary data format:
```
    rsid          chr      a1     a0     beta        beta_se      n_eff
    rs3131972	  1	   G	  A	 -0.005679   0.006108     100000
    rs3131969	  1	   G	  A	 -0.006049   0.006711     100000
    rs1048488	  1        T	  C	 -0.007662   0.006125     100000
```
The following columns are required for the GWAS summary data files:

 1. rsid: SNP ID, in the format of rsXXXX. If only the position information is available (in the training GWAS data and/or validation data), please impute rsid (e.g., using LD reference data with the same genome build).
 2. chr: shromosome, 1, 2, ..., or 22.
 3. beta: SNP effect. For binary traits, beta is the log of odds ratio (logOR) from logistic regressions.
 4. beta_se: standard error of beta.
 5. a1: effective (reference) allele (counted allele in regression), the allele which beta corresponds to.
 6. a0: other (alternative) allele. Note that sometimes it is referred to as "a2".
 7. n_eff: GWAS sample size by SNP. For binary traits, it is the effective sample sizes = 4 / (1 / N_control + 1 / N_case); and for continuous traits, it is the same as the total sample size.

Before running the MUSSEL pipeline, please consider applying the following quality control steps on the GWAS summary data files, if possible:

 1. Only keep the biallelic HapMap3 + MEGA SNPs (SNP IDs can be found in the second column of `ref_bim.txt`) to avoid troubles caused by reading huge files (e.g., > 8 million SNPs) in R.
 2. Remove SNPs with minor allele frequencies lower than 1% in all populations.
 3. Remove SNPs with small GWAS sample sizes (90% of median per SNP GWAS sample size).


__Note:__ 

there are several command lines that may need to be customized by users because of discrepancies in server:

- The linux commands for submitting jobs on server (lines 88 - 94 and line 110: "sbatch --mem=30G" in `LDpred2_jobs.R`, and lines 120 - 126 and lines 145 - 146 (e.g., "sbatch --mem=35G"), in `MUSS_jobs.R`), may need to be modified according to the server used (e.g., "module load conda_R", "module load R/4.3", etc.).
- For the commands in line 110 in `LDpred2_jobs.R`: "sbatch --mem=25G", the memory is required for running LDpred2 sequentially on two ancestry groups with NCORES=11 cores. If more ancestries are included, a larger memory (e.g., 40G for $K=5$) may need to requested (mainly for loading the LD information). Similarly, for the commands in lines 145 - 146 in `MUSS_jobs.R`: the required memory should be customized according to the number of training ancestry groups ($K$). The default memory requirement in MUSS_jobs.R is for jointly modeling 2 ancestry groups. For modeling more ancestry groups, the requested memory can be estimated as a linear function of $K$.


</br>


## MUSSEL Manual

MUSSEL workflow: 
<p align="center">
<img
  src="/img/MUSSEL_Workflow.png"
  title="MUSSEL Workflow"
  width=85% 
  height=85%>
</p>

### Step 0: run LDpred2 by ancestry

This step is to obtain the estimated causal SNP proportion ($p_k, k=1,2,\ldots,K$) and heritability ($h^2_k, k=1,2,\ldots,K$) parameters in LDpred2 for each of $K$ training ancestry groups. These parameters will be used to specify the prior causal SNP proportions and heritability parameters in [MUSS](#step-1:-muss).


LDpred2_jobs.R: submit LDpred2 jobs by chromosome.
```r
LDpred2_jobs.R --PATH_package --PATH_data --PATH_LDref --PATH_out --FILE_sst --pop --bfile_tuning --NCORES
```

LDpred2_tuning.R: obtain estimated LDpred2 parameters.
```r
LDpred2_tuning.R --PATH_package --PATH_out --PATH_plink --FILE_sst --pop --chrom 1-22 --bfile_tuning --pheno_tuning --bfile_testing --pheno_testing --testing --NCORES
```

### Step 1: MUSS

MUSS: a Bayesian model that jointly models the GWAS summary data across all training populations to obtain a total of $L \times K$ PRS models under $L$ different tuning parameter settings for $Pr⁡(δ_{1j},…,δ_{Kj})$ (functions of $p_k$s) and $\rho_{k_1,k_2}$s across all $K$ training populations.

```r
MUSS_jobs.R --PATH_package --PATH_data --PATH_LDref --PATH_out --FILE_sst --pop --LDpred2_params --chrom --bfile_tuning --NCORES
```

### Step 2: MUSSEL

For each target population, apply the Super Learning (SL) algorithm (default base learners: elastic net regression, ridge regression, and linear regression) to train an “optimal” linear combination of the ($L \times K$) PRS models, which we call the MUSSEL PRS model, based on the tuning set of the target population. 

Optional: the prediction performance of the final MUSSEL PRS model can be reported on an independent testing set, if the testing set is provided as an input.

```r
MUSSEL.R --PATH_package --PATH_out --PATH_plink --FILE_sst --pop --chrom --bfile_tuning --pheno_tuning --bfile_testing --pheno_testing --testing --NCORES
```

- PATH_package (required): path to the directory where the downloaded files (decompressed) are saved.

- PATH_data (required): path to the directory where the training data by ancestry group are saved.

- PATH_LDref (required): path to the directory where the LD reference data by ancestry group and chromosome are saved.

- PATH_out (required): path to the output directory where the results are saved.

- PATH_plink (required): path to plink2.

- FILE_sst (required): paths followed by file names of the population-specific GWAS summary statistics, separated by comma. Required columns: chr, rsid, pos, a0, a1, beta, beta_se, n_eff.

- pop (required: populations of the GWAS samples, separated by comma.

- chrom (required): the chromosome on which the model is fitted, input in the format of 1-22 or 1,2,3. Default: 1-22

- p: candidate values for tuning parameter p (causal SNP proportion). Default:
1.0e-05, 1.8e-05, 3.2e-05, 5.6e-05, 1.0e-04, 1.8e-04, 3.2e-04, 5.6e-04, 1.0e-03, 1.8e-03, 3.2e-03, 5.6e-03, 1.0e-02, 1.8e-02, 3.2e-02, 5.6e-02, 1.0e-01, 1.8e-01, 3.2e-01, 5.6e-01, 1.0e+00.

- H2: candidate values for tuning parameter H2 (heritability = H2 * h2_est from LDSC). Default: 0.3, 0.7, 1, 1.4.

- sparse: whether to consider a sparse model: 0, 1, or 0,1. Default: 0.

- bfile_tuning (required): path to PLINK binary input file prefix (excluding ".bed"/".bim"/".fam"") for tuning, save by chromosome.

- pheno_tuning (optional): path to phenotype file (PLINK format) for tuning.

- covar_tuning (optional): path to quantitative covariates (PLINK format) for tuning.

- testing (required): whether to perform testing in seperate dataset. Default: F.

- trait_type (required): Type of phenotype, continuous or binary. Default: 'continuous'.

- bfile_testing (optional): path to PLINK binary input file prefix (.bed/.bim/.fam) for testing, save by chromosome.

- pheno_testing (optional): path to phenotype file (PLINK format) for testing.

- covar_testing (optional): path to quantitative covariates (PLINK format) for testing.

- verbose: how much chatter to print: 0=nothing; 1=minimal; 2=all. Default: 1.

- cleanup: cleanup temporary files or not. Default: T.

- NCORES: how many cores to use. (Default: 13 for LDpred2_jobs.R, 5 for MUSS_jobs.R, and 1 for LDpred2_tuning.R and MUSSEL.R)

- LDpred2_params (required): path to the directory where the tuned LDpred2 parameters (population-specific causal SNP proportions, heritability and whether or not a sparse model is used) are saved, separated by comma.

- cors_additional (optional): additional candidate values for tuning parameter: genetic correlation across ancestry groups, example: 3 groups with label 1,2,3, want to add two additional settings: cor_setting1(1,2),cor_setting1(1,3),cor_setting1(2,3);cor_setting2(1,2),cor_setting2(1,3),cor_setting2(2,3).

- ps_additional (optional): typically not necessary. Additional candidate values for tuning parameter: ancestry-specific causal SNP proportions, example: 3 groups with label 1,2,3, want to add two additional settings: p1_setting1,p2_setting1,p3_setting1;p1_setting2,p2_setting2,p3_setting2.

- SL_library (optional): the base learners implemented in SuperLearner, separated by comma. Default: SL.glmnet,SL.ridge,SL.lm.

- linear_score (optional): whether the trained linear models will be saved. If not, only the Super Learner model will be saved. Note: some models in SL_library are non-linear. In this case, linear score file cannot be generated.

- target_pop (required): Target population (used to save output).

</br>


## Example
Download [example data](https://www.dropbox.com/scl/fi/bne781g2qsq67p0r9y9cl/example.zip?rlkey=6aw5tnfpbnjc3ieq1ee4do2ay&dl=0), decompress it by `tar -zxvf example.tar.gz` and save the files under the directory ${path_example}. Download the 1000 Genomes reference data and save the decompressed files in ${path_LDref}. Create a new folder `path_out` (e.g., in this example, `/dcs04/nilanjan/data/jjin/MUSSEL/test`) to save the output. Run the example code below with your own data directories and check if the results/outputs (saved in ${path_out}) are consistent with the results/outputs here: [example results](https://www.dropbox.com/scl/fi/44z9rkku5nsc45yvdm6hk/MEBayesSL_example_data_results.zip?rlkey=rrmgjwu4at0nbr4spjj6836l5&dl=0).

```r 
module load R

package='/dcs04/nilanjan/data/jjin/MUSSEL'
path_data='/dcs04/nilanjan/data/jjin/example'
path_LDref='/dcs04/nilanjan/data/jjin/LD_1kg'
path_out='/dcs04/nilanjan/data/jjin/MUSSEL/test'
path_plink='/dcl01/chatterj/data/jin/software/plink2'
target_pop='EUR,AFR'
trait_type='continuous'
```
Note: load the R version for which the required R packages were installed, in this example, R Version 4.2.2 Patched (2023-03-01 r83924).


### Step 1: Run LDpred2 by chromosome (22 jobs, each for one chromosome). 

In each job, the algorithm will run under different tuning parameter settings in parallel.
Note: as side products, $K$ LDpred2 PRS models trained based on GWAS data for each ancestry group will be generated by this step.

``` r
Rscript ${package}/R/LDpred2_jobs.R \
--PATH_package ${package} \
--PATH_data ${path_data} \
--PATH_LDref ${path_LDref} \
--PATH_out ${path_out} \
--FILE_sst ${path_data}/summdata/EUR.txt,${path_data}/summdata/AFR.txt \
--pop EUR,AFR \
--bfile_tuning ${path_data}/sample_data/EUR/tuning_geno,${path_data}/sample_data/AFR/tuning_geno \
--NCORES 11

```




### Step 2: Obtain tuned parameters from LDpred2

Wait until all LDpred2 jobs are completed to run this step. Tuned LDpred2 effect size estimates and the optimal tuning parameters are saved in ${path_out}. This step also generates the LDpred2 PRS models for each ancestry group as by products, which will be saved in `${path_out}/LDpred2/{race}_LDpred2_beta.txt`.

``` r
Rscript ${package}/R/LDpred2_tuning.R \
--PATH_package ${package} \
--PATH_out ${path_out} \
--PATH_plink ${path_plink} \
--FILE_sst ${path_data}/summdata/EUR.txt,${path_data}/summdata/AFR.txt \
--pop EUR,AFR \
--chrom 1-22 \
--bfile_tuning ${path_data}/sample_data/EUR/tuning_geno,${path_data}/sample_data/AFR/tuning_geno \
--pheno_tuning ${path_data}/sample_data/EUR/pheno.txt,${path_data}/sample_data/AFR/pheno.txt \
--covar_tuning ${path_data}/sample_data/EUR/covar.txt,${path_data}/sample_data/AFR/covar.txt \
--bfile_testing ${path_data}/sample_data/EUR/testing_geno,${path_data}/sample_data/AFR/testing_geno \
--pheno_testing ${path_data}/sample_data/EUR/pheno.txt,${path_data}/sample_data/AFR/pheno.txt \
--covar_testing ${path_data}/sample_data/EUR/covar.txt,${path_data}/sample_data/AFR/covar.txt \
--trait_type ${trait_type} \
--testing TRUE \
--NCORES 1

```
Note: 
 1. the `pheno.txt` dataset contains phenotype information for both tuning and testing individuals. 
 2. For continuous traits, instead of using the `--covar_tuning` and `--covar_testing` option to adjust for covariates, an alternative approach is to first regress the phenotype on covariates, then save the residual (adjusted phenotype data) in `pheno.txt`.

### Step 3: Run MUSS by chromosome (submit 22 jobs simultaneously, each for one chromosome). 

``` r
Rscript ${package}/R/MUSS_jobs.R \
--PATH_package ${package} \
--PATH_data ${path_data} \
--PATH_LDref ${path_LDref} \
--PATH_out ${path_out} \
--FILE_sst ${path_data}/summdata/EUR.txt,${path_data}/summdata/AFR.txt \
--pop EUR,AFR \
--LDpred2_params ${path_out}/LDpred2/EUR_optim_params.txt,${path_out}/LDpred2/AFR_optim_params.txt \
--chrom 1-22 \
--bfile_tuning ${path_data}/sample_data/EUR/tuning_geno,${path_data}/sample_data/AFR/tuning_geno \
--NCORES 5

```

### Step 4: Combine PRS models generated under different parameter settings with a Super Learner (SL) algorithm to obtain the final ensembled MUSSEL PRS model. 

If a testing dataset is provided, the prediction $R^2$ or $AUC$ of the final MUSSEL PRS model will be reported on the testing set.

``` r
Rscript ${package}/R/MUSSEL.R \
--PATH_package ${package} \
--PATH_out ${path_out} \
--PATH_plink ${path_plink} \
--pop EUR,AFR \
--target_pop ${target_pop} \
--chrom 1-22 \
--bfile_tuning ${path_data}/sample_data/EUR/tuning_geno,${path_data}/sample_data/AFR/tuning_geno \
--pheno_tuning ${path_data}/sample_data/EUR/pheno.txt,${path_data}/sample_data/AFR/pheno.txt \
--covar_tuning ${path_data}/sample_data/EUR/covar.txt,${path_data}/sample_data/AFR/covar.txt \
--bfile_testing ${path_data}/sample_data/EUR/testing_geno,${path_data}/sample_data/AFR/testing_geno \
--pheno_testing ${path_data}/sample_data/EUR/pheno.txt,${path_data}/sample_data/AFR/pheno.txt \
--covar_testing ${path_data}/sample_data/EUR/covar.txt,${path_data}/sample_data/AFR/covar.txt \
--trait_type ${trait_type} \
--testing TRUE \--NCORES 1

```


## Questions

Please report any issues on the Issues page, I will respond as soon as possible. For a quicker response, please contact Jin.Jin@Pennmedicine.upenn.edu.


## Citation

__Jin, J.__, Zhan, J., Zhang, J., Zhao, R., O’Connell J, Jiang, Y., Buyske, S., Gignoux, C., Haiman, C., Kenny, E.E., Kooperberg, C., et al. MUSSEL: Enhanced Bayesian polygenic risk prediction leveraging information across multiple ancestry groups. Cell Genomics 4(4), 100539, 2024. [Link](https://www.cell.com/cell-genomics/fulltext/S2666-979X(24)00095-8)


