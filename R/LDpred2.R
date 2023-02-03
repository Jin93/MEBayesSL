rm(list=ls())
suppressMessages(library("optparse"))
suppressMessages(library("bigsnpr"))
suppressMessages(library("bigreadr"))
suppressMessages(library("readr"))
suppressMessages(library("stringr"))
suppressMessages(library("caret"))
suppressMessages(library("Rcpp"))
suppressMessages(library("RcppArmadillo"))
suppressMessages(library("inline"))
suppressMessages(library("doMC"))
suppressMessages(library("foreach"))

option_list = list(
  make_option("--PATH_package", action="store", default=NA, type='character',
              help="Path to the directory where the downloaded files (decompressed) are saved [required]"),
  make_option("--PATH_ref", action="store", default=NA, type='character',
              help="Path to the directory where the LD reference data by ancestry group and chromosome are saved [required]"),
  make_option("--PATH_out", action="store", default=NA, type='character',
              help="Path to the output directory where the results are saved [required]"),
  make_option("--FILE_sst", action="store", default=NA, type='character',
              help="Paths followed by file names of the population-specific GWAS summary statistics, separated by comma [required] [required columns: chr, rsid, pos, a0, a1, beta, beta_se, n_eff]"),
  make_option("--pop", action="store", default=NA, type='character',
              help="Populations of the GWAS samples, separated by comma [required]"),
  make_option("--chrom", action="store", default=NA, type='character',
              help="The chromosome on which the model is fitted [required]"),
  
  make_option("--p", action="store", default=paste(signif(seq_log(1e-4, 1, length.out = 17), 2), collapse = ','), type='character',
              help="Candidate values for tuning parameter p (causal SNP proportion) [default: %default]"),
  make_option("--H2", action="store", default=paste(c(0.7, 1, 1.4), collapse = ','), type='character',
              help="Candidate values for tuning parameter H2 (heritability = H2 * h2_est from LDSC) [default: %default]"),
  make_option("--sparse", action="store", default='0', type='character',
              help="Whether to consider a sparse model: 0, 1, or 0,1 [default: %default]"),
  make_option("--bfile_tuning", action="store", default=NA, type='character',
              help="Path to PLINK binary input file prefix (minus .bed/.bim/.fam) for tuning, save by chromosome [required]"),

  make_option("--verbose", action="store", default=1, type="integer",
              help="How much chatter to print: 0=nothing; 1=minimal; 2=all [default: %default]"),
  make_option("--cleanup", action="store", default=T, type="logical",
              help="Cleanup temporary files or not [default: %default]"),
  make_option("--NCORES", action="store", default=17, type="integer",
              help="How many cores to use [default: %default]")
)
opt = parse_args(OptionParser(option_list=option_list))

if ( opt$verbose == 2 ) { SYS_PRINT = F } else { SYS_PRINT = T }

suppressWarnings(dir.create(opt$PATH_out))

races = str_split(opt$pop,",")[[1]]; K <- length(races)
sumdata_paths = str_split(opt$FILE_sst,",")[[1]]
ref_paths <- paste0(opt$PATH_ref,"/",races)
out_paths <- paste0(opt$PATH_out,"/",races)

opt$chrom <- gsub("-",":",opt$chrom)
eval(parse(text=paste0("chrom = c(",opt$chrom,")")))

bfile_tuning_vec <- str_split(opt$bfile_tuning,",")[[1]]

# Perform i/o check:
files <- NULL
files <- c(files, sumdata_paths)
for(mmm in 1:K){ files <- c(files, paste(bfile_tuning_vec[mmm],c(".bed",".bim",".fam"),sep='')) }

for ( f in files ) {
  if ( !file.exists(f) ){
    cat( "ERROR: ", f , " files for tuning data do not exist\n" , sep='', file=stderr() )
    q()
  }
}
rm(list="files")

NCORES <- opt$NCORES

source(paste0(opt$PATH_package,"/R/source-functions.R"))
ldr = 3/1000 # default ld radius in LDpred2

ref <- fread2(paste0(opt$PATH_package,"/ref_bim.txt"))
chr <- as.numeric(opt$chrom[1])
p_seq <- as.numeric(strsplit(opt$p, split = ',')[[1]])
H2s = as.numeric(strsplit(opt$H2, split = ',')[[1]])
sparse = as.numeric(strsplit(opt$sparse, split = ',')[[1]])



for(mmm in 1:K){
  race <- races[mmm]
  sumdata_path <- sumdata_paths[mmm]
  ref_path <- ref_paths[mmm]
  out_path <- out_paths[mmm]
  
  bfile_tuning <- bfile_tuning_vec[mmm]

  suppressWarnings(dir.create(paste0(out_path)))
  suppressWarnings(dir.create(paste0(out_path, "/tmp")))
  suppressWarnings(dir.create(paste0(out_path, "/tmp/ref_files")))
  suppressWarnings(dir.create(paste0(out_path, "/tmp/beta_files")))
  suppressWarnings(dir.create(paste0(out_path, "/tmp/beta_files/beta_in_all_settings_bychrom")))

  if ( opt$verbose >= 1 ) {
    cat(paste0("\n*************************"))
    cat(paste0("\n**** LDpred2 on ",race," ****"))
    cat(paste0("\n*************************\n"))
  }
    
  ########################################################################
  ########################################################################
  if ( opt$verbose >= 1 ) cat(paste0('\n** Step 1: data preparation **\n'))
  
  ############
  ## Load GWAS summary data, only keep SNPs that are present in both LD reference samples and tuning samples (according to rsid)
  sum.raw <-bigreadr::fread2(sumdata_path)
  sum.raw <- sum.raw[sum.raw$rsid %in% ref$V2,]
  tun <- bigreadr::fread2(paste0(bfile_tuning,'.bim'))
  sum.raw <- sum.raw[sum.raw$rsid %in% tun$V2,]
  ref_tmp <- ref[match(sum.raw$rsid, ref$V2),]
  
  temfile = paste0(ref_path, '/chr',chr,'.bk')
  system(paste0('rm -rf ',temfile))
  temfile = paste0(ref_path, '/chr',chr,'.rds')
  system(paste0('rm -rf ',temfile))
  temfile = paste0(ref_path, '/chr',chr,'.bk')
  if (!file.exists(temfile)){
    snp_readBed(paste0(ref_path, '/chr',chr,'.bed'))
  }
  obj.bigSNP <- snp_attach(paste0(ref_path, '/chr',chr,'.rds'))
  
  G   <- obj.bigSNP$genotypes
  CHR <- obj.bigSNP$map$chromosome
  POS <- obj.bigSNP$map$physical.pos
  if (! 'pos' %in% colnames(sum.raw)) sum.raw$pos = NA
  sumstats = sum.raw[sum.raw$chr == chr,c("chr", "rsid", "pos", "a0", "a1", "beta", "beta_se", "n_eff")]
  
  set.seed(2020)
  map <- obj.bigSNP$map[-c(3)]
  names(map) <- c("chr", "rsid", "pos", "a0", "a1")
  info_snp <- snp_match(sumstats, map, strand_flip = T, join_by_pos = F)
  rownames(info_snp) = info_snp$rsid
  
  POS2 <- snp_asGeneticPos(CHR, POS, dir = paste0(out_path, "/tmp/ref_files"), ncores = 3)
  ## indices in info_snp
  ind.chr <- which(info_snp$chr == chr)
  df_beta <- info_snp[ind.chr, c("beta", "beta_se", "n_eff")]
  if ( opt$verbose == 2 ) cat(paste0("*** CHR ",chr,": ", nrow(df_beta)," SNPs are included in the analysis *** \n"))
  ## indices in G
  ind.chr2 <- info_snp$`_NUM_ID_`[ind.chr]
  
  ## Compute correlation
  corr0 <- snp_cor(G, ind.col = ind.chr2, ncores = 3, infos.pos = POS2[ind.chr2], size = ldr) # default
  corr <- bigsparser::as_SFBM(as(corr0, "dgCMatrix"))
  
  
  ldsc <- snp_ldsc2(corr0, df_beta)
  h2_est <- abs(ldsc[["h2"]])
  # grid of models:
  H2_seq = signif(abs(h2_est) * H2s, 3)
  h2_seq <- signif(abs(h2_est) * H2s, 3)
  params <- expand.grid(p = p_seq, h2 = h2_seq, sparse = sparse)
  
  rm(sum.raw, tun, ref_tmp, sumstats, corr0)
  
  if ( opt$verbose >= 1 ) cat("\n** Step 2: run LDpred2 under all tuning parameter settings **\n")
  
  beta_grid <- snp_ldpred2_grid(corr, df_beta, params, ncores = NCORES)
  beta_grid = as.data.frame(beta_grid)
  beta_grid[is.na(beta_grid)] = 0; beta_grid[abs(beta_grid) > 5] = 0
  rownames(beta_grid) = info_snp$rsid
  beta_grid = cbind(info_snp$rsid, info_snp$a0, info_snp$a1, beta_grid)
  colnames(beta_grid) = c(c('marker.ID', 'a0', 'a1'),paste0('e',1:nrow(params)))
  
  beta_grid = beta_grid[,-2]
  write.table(beta_grid,file = paste0(out_path, '/tmp/beta_files/beta_in_all_settings_bychrom/ldpred2effect-chr',chr,'.txt'),
              col.names = T,row.names = F,quote=F)
  
  # The rows are in the same order as the columns in beta_grid.
  params <- expand.grid(p = p_seq, h2 = H2s, sparse = sparse)
  fwrite2(params, paste0(out_path, '/tmp/beta_files/beta_in_all_settings_bychrom/params_chr',chr,'.txt'), col.names = T, sep="\t", nThread=1)
  
  if ( opt$verbose >= 1 ) cat(paste0("\n** Completed: results saved in ", paste0(out_path, "/tmp/beta_files/beta_in_all_settings_bychrom/params_chr",chr,".txt","**\n")))
  ############
  rm(list=c("corr"))
}

  
