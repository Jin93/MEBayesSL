# MUSSEL

MUSSEL is an R-based command line tool for implementing MUSSEL, a method for developing ancestry-specific polygenic risk score (PRS) that integrates information from GWAS summary statistics and external LD reference data from multiple populations (ancestry groups). MUSSEL infers SNP effect sizes via a Bayesian model with an induced prior correlation structure across populations followed by an ensemble learning step with the Super Learner. As intermediate products, LDpred2 PRS models trained separately on GWAS data for each population will also be generated.

Please refer to the [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10120638/) or contact Jin Jin (Jin.Jin@Pennmedicine.upenn.edu) for details.
</br>



## Version History
- [ ] __April 8, 2024__ Updated source code to incorporate updates for the [LDpred2](https://privefl.github.io/bigsnpr/articles/LDpred2.html) pipeline: now LDpred2_jobs.R will only submit one job per ancestry to the server instead of 22 (by chromosome) in the previous version; added more details regarding recommendation for LD reference panel.

- [ ] __Feb 3, 2023__ First version of MUSSEL was made publicly available on Github.
</br>


## Getting Started

- Download or clone the Github repository by `git clone https://github.com/Jin93/MUSSEL.git`. From now on, we will refer to the folder as `/MUSSEL/` for simplicity.

- Download the `ref_bim.txt` from [this link](https://www.dropbox.com/s/58uzwqewxv34wal/ref_bim.txt?dl=0) and save it under `/MUSSEL/`.


- Download and decompress the LD reference data and save the decompressed folders in ${path_LDref}.

The LD reference data contains SNP information and LD estimates by LD block for genetic variants that are in the [HapMap3](https://www.broadinstitute.org/medical-and-population-genetics/hapmap-3) plus [MEGA Chip](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5156387/) list. Note: in some scenarios, the training GWAS sample for a population consists of multiple ancestry groups, in this case, ideally, a customized LD reference dataset should be created for this population with matched ancestral composition. The code for constructing such LD reference dataset can be requested (Jin.Jin@Pennmedicine.upenn.edu).

Please choose one of the two LD reference panels according to the training GWAS sample sizes (please see detailed description for each LD panel). Each reference data contains two folders: 

`./raw/`: raw LD reference genotype data (.bim, .bed,. fam), which are input files for estimating LD matrices in LDpred2 and for an intermediate step of summarizing LD information in the [MUSS](#mussel-manual) step.

`./LD/`: Precalculated LD matrices and SNP information by LD block, which are input files in MUSS.


#### 1. UK Biobank LD reference data

- 10,000 EUR, 4,585 AFR, 687 AMR, 1,010 EAS, 5,427 SAS.
- Recommended when training GWAS sample sizes are relatively large, e.g., N<sub>GWAS</sub> > 50K for at least two ancestry groups.

[EUR reference data](https://www.dropbox.com/scl/fi/09yd12dest1tqxkt8p8ch/EUR.zip?rlkey=774vb1e5d6hfnyucilx160cyo&dl=0) (~13.15G): `tar -zxvf EUR.tar.gz`

[AFR reference data](https://www.dropbox.com/scl/fi/jfymih83anr2vuevmfqok/AFR.zip?rlkey=r1lxpn1fnbk98ssf8f8ji4xkk&dl=0) (~11.59G): `tar -zxvf AFR.tar.gz`

[AMR reference data](https://www.dropbox.com/s/2ba4tsbhz03rg83/AMR.zip?dl=0) (~4.88G): `tar -zxvf AMR.tar.gz`

[EAS reference data](https://www.dropbox.com/s/uofu788707dp4xv/EAS.zip?dl=0) (~4.27G): `tar -zxvf EAS.tar.gz`

[SAS reference data](https://www.dropbox.com/scl/fi/o635c86ylthbl3omfetbu/SAS.zip?rlkey=ot396toxl0phaiae15cnpbeyn&dl=0) (~11.44G): `tar -zxvf SAS.tar.gz`

Note: PRS trained using the larger UKBB LD reference data is usually more powerful than PRS trained using the 1000G reference data, especially with a sufficiently large discovery GWAS sample size (e.g., >100K).

#### 2. LD reference data constructed based on the 1000 Genomes Project phase 3 samples 

- 498 EUR, 659 AFR, 347 AMR, 503 EAS, 487 SAS.
- Recommended when GWAS training sample sizes are relatively small, e.g., N<sub>GWAS</sub> < 50K for all ancestry groups.

[EUR reference data](https://www.dropbox.com/s/wvxh4yqthm8m7uf/EUR.zip?dl=0) (~6.73G): `tar -zxvf EUR.tar.gz`

[AFR reference data](https://www.dropbox.com/s/iwqg65uieevfzj2/AFR.zip?dl=0) (~7.69G): `tar -zxvf AFR.tar.gz`

[AMR reference data](https://www.dropbox.com/s/mev5zyf4x6m076q/AMR.zip?dl=0) (~8.80G): `tar -zxvf AMR.tar.gz`

[EAS reference data](https://www.dropbox.com/s/o28mlovtakv5n7v/EAS.zip?dl=0) (~5.63G): `tar -zxvf EAS.tar.gz`

[SAS reference data](https://www.dropbox.com/s/idp02rgl8xv379b/SAS.zip?dl=0) (~2.60G): `tar -zxvf SAS.tar.gz`


- Install [PLINK1.9](https://www.cog-genomics.org/plink/) and [PLINK2](https://www.cog-genomics.org/plink/2.0/).

- Launch R and install required libraries:

``` r
install.packages(c('optparse','bigreadr','bigsnpr','bigparallelr', 'bigmemory','stringr','caret','Rcpp', 'RcppArmadillo','RcppTN','inline','doMC','foreach','doParallel','data.table','readr','MASS','reshape','parallel',
'devtools','genio','dplyr','pryr','Matrix','lavaan','xtable','SuperLearner'))
```

- Prepare input data files:

Please use the `example_data` in [Example](#example) as a reference to prepare  the input data files. Please note that the SNPs are matched across the training GWAS data, LD reference data, and validation data by RSID. If only the position information is available for SNPs in the training GWAS data and validation data, please impute the RSID (e.g., using LD reference data with the same genome build) before running the MUSSEL pipeline.

An example of the GWAS summary data format:
```
    rsid           chr     a1     a0     beta        beta_se      n_eff
    rs3131972	  1	   G	  A	 -0.005679   0.006108     100000
    rs3131969	  1	   G	  A	 -0.006049   0.006711     100000
    rs1048488	  1        T	  C	 -0.007662   0.006125     100000
```
The following columns are required for the GWAS summary data files:

 1. rsid: SNP ID in the format of rsXXXX.
 2. chr: chromosome, 1,2,...,22.
beta: SNP effect. Note that for binary traits, beta is the coefficient in logistic regression, i.e. log(OR).
beta_se: Standard error of beta.
a1: effective allele (counted allele in regression).
a0: alternative allele (non-A1 allele).
n_eff: Sample size per variant. Note that for binary traits, it is effective sample sizes = 4 / (1 / N_control + 1 / N_case); and for continuous traits, it is simply the sample size.
Note that the summary statistics files are suggested to be cleaned as follows before using:

Keep variants in reference panels to avoid troubles caused by reading huge files. The rsid of variants in reference panels can be found in ref_bim.txt.
Alleles are suggested to match to reference panels to avoid flipping strands. The alleles of variants in reference panels can be found in ref_bim.txt.
For each population, remove variants whose allele frequencies are less than 1%. If your summary statistics does not contain information of population-specific allele frequencies, we recommanded to use that in 1000G, which can be found in ref_af.txt.
Remove variants with small sample size (90% of median sample size per variant).

__Note:__ 

there are several command lines that need to be customized by users because of discrepancies in server:

- The command line, "module load R" in `LDpred2_jobs.R` and `MUSS_jobs.R`, for loading R on server, may need to be modified (e.g., "module load conda_R", "module load R/4.3", etc.)
- The command lines on lines 119 - 120 in `LDpred2_jobs.R` and lines 143 - 144 in `MUSS_jobs.R`, e.g., "sbatch --mem=23G", may need to be modified. Note: the memory required for MUSS_jobs.R should be customized according to the number of training ancestry groups. The default memory requirement in MUSS_jobs.R is for jointly modeling 5 ancestry groups. For modeling two ancestry groups, the memory requested can be set to about half of the default values.


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

### Step 0 

Run LDpred2 by ancestry, obtain the estimated causal SNP proportion ($p_k, k=1,2,\ldots,K$) and heritability ($h^2_k, k=1,2,\ldots,K$) parameters in LDpred2 for each of $K$ training ancestry groups. These parameters will be used to specify the prior causal SNP proportions and heritability parameters in [MUSS](#step-1:-muss).


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
Download [example data](https://www.dropbox.com/s/xxw3t17k66il3k5/example.tar.gz?dl=0), decompress it by `tar -zxvf example.tar.gz` and save the files under the directory ${path_example}. Download the 1000 Genomes reference data and save the decompressed files in ${path_LDref}. Create a new folder `path_out` (e.g., in this example, `/dcs04/nilanjan/data/jjin/MUSSEL/test`) to save the output. Run the example code below with your own data directories and check if the results are consistent with the results here: [example results](https://www.dropbox.com/s/hjmqghn2jva0950/MUSSEL_example_data_results.zip?dl=0).

```r 
module load R

package='/dcs04/nilanjan/data/jjin/MUSSEL'
path_data='/dcs04/nilanjan/data/jjin/example'
path_LDref='/dcs04/nilanjan/data/jjin/LD_1kg'
path_out='/dcs04/nilanjan/data/jjin/MUSSEL/test'
path_plink='/dcl01/chatterj/data/jin/software/plink2'
target_pop='AFR'
```
Note: load the R version for which the required R packages were installed, in this example, R Version 4.2.2 Patched (2023-03-01 r83924).


### Step 1: Run LDpred2 by chromosome (by submitting 22 jobs simultaneously, each for one chromosome). 

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

### Step 3: Run MUSS by chromosome (by submitting 22 jobs simultaneously, each for one chromosome). In each job, the algorithm will run under different tuning parameter settings in parallel.

``` r
Rscript ${package}/R/MUSS_jobs.R \
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

### Step 4: Combine PRS models generated under different parameter settings with a Super Learner (SL) algorithm to obtain the final ensembled MUSSEL PRS model. Here, with the testing dataset provided, the prediction $R^2$ of the final MUSSEL PRS model is reported on the testing set.

``` r
Rscript ${package}/R/MUSSEL.R \
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


## Questions

Please report any issues on the Issues page, I will respond as soon as possible. For a quicker response, please contact Jin.Jin@Pennmedicine.upenn.edu.


## Citation

Jin, Jin, et al. "MUSSEL: enhanced Bayesian polygenic risk prediction leveraging information across multiple ancestry groups". bioRxiv. [Preprint](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10120638/)


