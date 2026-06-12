library(limma)  
library(edgeR)
library(DESeq2)
library(glue)
library(optparse)
library(BiocParallel)
library(dplyr)

##################################################### 日志辅助区 #############################################################

# 统一日志打印函数
log_msg <- function(level = "INFO", message = "") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(glue("[{timestamp}] [{level}] {message}\n"))
}

##################################################### 参数解析与校验 #############################################################

option_list <- list(
  make_option(c("--bed"), type="character", help="Path to BED file (Required)"),
  make_option(c("--cov"), type="character", help="Path to covariates file (Required)"),
  make_option(c("--normalization"), type = "character", help="Normalization method: ratio_norma or raw_count_norma (Required)"),
  make_option(c("--method"), type = "character", default="pca", help="Method to use: pca or peer [default: %default]"),
  make_option(c("--output"), type = "character", help="Output directory (Required)"),
  make_option(c("--prefix"), type = "character", help="Output prefix (Required)"),
  make_option(c("--threads"), type = "integer", default=1, help="Number of threads [default: %default]"),
  
  # PCA 特异参数
  make_option(c("--tier"), type="character", default=NULL, help="PCA specific: tier vector, e.g., '-5,0,5' [default: 0]"),
  
  # PEER 特异参数
  make_option(c("--peer_factors"), type = "integer", default=NULL, help="PEER specific: number of peer factors (Required for PEER)"),
  make_option(c("--peer_iter"), type = "integer", default=NULL, help="PEER specific: max iteration times [default: 1000]")
)

opt <- parse_args(OptionParser(option_list=option_list))

log_msg("INFO", "==== Start Parameter Validation ====")

# 检查常规必要参数是否存在
required_args <- c("bed", "cov", "normalization", "output", "prefix")
for (arg in required_args) {
  if (is.null(opt[[arg]])) {
    log_msg("ERROR", glue("Argument --{arg} is required."))
    stop(glue("Error: Argument --{arg} is required."))
  }
}

# 根据 method 严格限制和校验特异性参数
if (opt$method == "pca") {
  if (!is.null(opt$peer_factors) || !is.null(opt$peer_iter)) {
    log_msg("ERROR", "Cannot specify PEER parameters (--peer_factors, --peer_iter) when --method is 'pca'.")
    stop("Error: Cannot specify PEER parameters when --method is 'pca'.")
  }
  if (is.null(opt$tier)) {
    opt$tier <- "0"
  }
  tier <- as.numeric(strsplit(opt$tier, ",")[[1]])
  peer_factors_val <- 0
  peer_iter_val <- 0

} else if (opt$method == "peer") {
  if (!is.null(opt$tier)) {
    log_msg("ERROR", "Cannot specify PCA parameters (--tier) when --method is 'peer'.")
    stop("Error: Cannot specify PCA parameters when --method is 'peer'.")
  }
  if (is.null(opt$peer_factors)) {
    log_msg("ERROR", "--peer_factors is required when --method is 'peer'.")
    stop("Error: --peer_factors is required when --method is 'peer'.")
  }
  if (is.null(opt$peer_iter)) {
    opt$peer_iter <- 1000
  }
  peer_factors_val <- opt$peer_factors
  peer_iter_val <- opt$peer_iter
  tier <- 0 
  
} else {
  log_msg("ERROR", "--method must be either 'pca' or 'peer'.")
  stop("Error: --method must be either 'pca' or 'peer'.")
}

# 打印最终确定的核心参数
log_msg("INFO", glue("Input BED file: {opt$bed}"))
log_msg("INFO", glue("Input Covariates file: {opt$cov}"))
log_msg("INFO", glue("De-confounding Method: {opt$method}"))
log_msg("INFO", glue("Normalization Method: {opt$normalization}"))
log_msg("INFO", glue("Threads allocated: {opt$threads}"))
if (opt$method == "pca") log_msg("INFO", glue("PCA tier vector: {paste(tier, collapse=',')}"))
if (opt$method == "peer") log_msg("INFO", glue("PEER Factors: {peer_factors_val} | Max Iterations: {peer_iter_val}"))
log_msg("SUCCESS", "Parameter validation passed successfully.")

# 环境配置
register(MulticoreParam(workers = opt$threads))
Sys.setenv(
  OMP_NUM_THREADS = opt$threads,
  MKL_NUM_THREADS = opt$threads,
  OPENBLAS_NUM_THREADS = opt$threads
)

##################################################### 函数定义区 #############################################################

RINT <- function(x){
  r <- rank(x, na.last = "keep", ties.method = "random")
  qnorm((r - 3/8) / (sum(!is.na(x)) + 1/4))
}

qn_rint_matrix <- function(expr_mat){
  expr_qn <- normalizeBetweenArrays(expr_mat, method = "quantile")
  expr_rint <- t(apply(expr_qn, 1, RINT))
  rownames(expr_rint) <- rownames(expr_mat)
  colnames(expr_rint) <- colnames(expr_mat)
  return(expr_rint)
}

ratio_norma <- function(bed){
  log_msg("INFO", "Applying Ratio Quantile Normalization + RINT...")
  bed_info <- bed[, 1:4]
  expr_mat <- as.matrix(bed[, -c(1:4)])
  rownames(expr_mat) <- bed$gene_id
  expr_rint <- qn_rint_matrix(expr_mat)
  bed <- cbind(bed_info, expr_rint)
  bed <- bed[order(as.numeric(bed$chrom), bed$start), ]
  return(bed)
}

raw_count_norma <- function(bed){
  log_msg("INFO", "Applying EdgeR TMM + CPM Log2 + RINT Normalization...")
  gene_info <- bed[, 1:4]
  counts <- as.matrix(bed[, -c(1:4)])
  dge <- DGEList(counts = counts)
  dge <- calcNormFactors(dge, method = "TMM")
  logcpm <- cpm(dge, log = TRUE, prior.count = 1)
  expr_rint <- t(apply(logcpm, 1, RINT))
  bed_rint <- cbind(gene_info, expr_rint)
  return(bed_rint)
}

convert_to_dummy <- function(df) {
  log_msg("INFO", "Converting categorical covariates to dummy variables...")
  id_col <- df[, 1, drop = FALSE]
  other_cols <- df[, -1, drop = FALSE]
  result_cols <- list()
  for (col_name in colnames(other_cols)) {
    col_data <- other_cols[[col_name]]
    if (is.numeric(col_data)) {
      result_cols[[col_name]] <- col_data
    } else {
      unique_vals <- unique(col_data)
      n_levels <- length(unique_vals)
      log_msg("INFO", glue("  Categorical covariate found: '{col_name}' with {n_levels} levels."))
      if (n_levels == 2) {
        level1 <- unique_vals[1]
        level2 <- unique_vals[2]
        dummy_vec <- ifelse(col_data == level2, 1, 0)
        new_col_name <- paste0(col_name, "_", level2)
        result_cols[[new_col_name]] <- dummy_vec
      } else {
        ref_level <- unique_vals[1]
        other_levels <- unique_vals[-1]
        for (level in other_levels) {
          dummy_vec <- ifelse(col_data == level, 1, 0)
          new_col_name <- paste0(col_name, "_", level)
          result_cols[[new_col_name]] <- dummy_vec
        }
      }
    }
  }
  result <- bind_cols(id_col, as.data.frame(result_cols))
  log_msg("SUCCESS", glue("Dummy conversion complete. Total covariate features: {ncol(result) - 1}"))
  return(result)
}

getResiduals<-function(dataResponse,dataPredictors){
  vecOfOnes<-rep(1,nrow(dataResponse)) 
  X<-cbind(vecOfOnes,as.matrix(dataPredictors)) 
  betaEstimate<-solve(t(X)%*%X)%*%t(X)%*%as.matrix(dataResponse) 
  toReturn<-as.matrix(dataResponse)-X%*%betaEstimate
  return(toReturn)
}

PCA_getResiduals<-function(normalized_bed,cov,tier){
    log_msg("INFO", "Starting PCA decomposition...")
    expr<-t(normalized_bed[,-(1:4)])
    metadata<-normalized_bed[,(1:4)]
    prcompResult<-prcomp(expr,center=TRUE,scale.=TRUE)
    PCs<-prcompResult$x
    
    # 保存 PC vectors
    pc_file <- glue("{opt$output}/{opt$prefix}_PCA_vectors.txt")
    write.table(data.frame(sample = rownames(PCs), PCs, check.names = FALSE), 
                file = pc_file, sep = "\t", row.names = FALSE, quote = FALSE)
    log_msg("SUCCESS", glue("Saved PCA vectors to: {pc_file}"))

    log_msg("INFO", "Estimating optimal number of PCs via PCAForQTL...")
    resultRunElbow<-PCAForQTL::runElbow(prcompResult=prcompResult)
    RNGkind("L'Ecuyer-CMRG")
    set.seed(1)
    resultRunBE<-PCAForQTL::runBE(expr,B=20,alpha=0.05,mc.cores=opt$threads)
    K_BE<-resultRunBE$numOfPCsChosen
    log_msg("SUCCESS", glue("Optimal PC number chosen by Buja & Eyuboglu (BE) algorithm: K_pc = {K_BE}"))
    
    if(length(tier) == 1 && tier == 0) {
        log_msg("INFO", glue("Processing default tier (PC={K_BE})..."))
        PCsTop<-PCs[,1:K_BE]
        knownCovariatesFiltered<-PCAForQTL::filterKnownCovariates(cov,PCsTop,unadjustedR2_cutoff=0.9)
        PCsTop<-scale(PCsTop)
        covariatesToUse<-cbind(knownCovariatesFiltered,PCsTop)
        Residuals=getResiduals(expr,covariatesToUse)
        final_bed<- cbind(metadata, t(Residuals))
        return(setNames(list(final_bed), paste0("tier_PC", K_BE))) 
    }
    
    residual_list <- list()
    for(t in tier){
        nPC <- K_BE + t
        nPC <- max(1, min(nPC, ncol(PCs)))
        log_msg("INFO", glue("Processing customized tier: K_pc ({K_BE}) + ({t}) = {nPC} PCs..."))
        PCsTop <- scale(PCs[,1:nPC])
        knownCovariatesFiltered<-PCAForQTL::filterKnownCovariates(cov,PCsTop,unadjustedR2_cutoff=0.9)
        covariatesToUse <- cbind(knownCovariatesFiltered, PCsTop)
        Residuals <- getResiduals(expr, covariatesToUse)
        residual_list[[paste0("tier_", nPC)]] <- cbind(metadata, t(Residuals))
    }
    return(residual_list)
}

saveResiduals <- function(residual_list,output_dir,prefix){
  for(tier_name in names(residual_list)){
      nPC <- sub("tier_", "", tier_name)
      outfile <- glue("{output_dir}/{prefix}_tier_PC{nPC}.bed")
      write.table(residual_list[[tier_name]], file=outfile, sep="\t", row.names=FALSE, quote=FALSE)
      log_msg("SUCCESS", glue("Exported PCA Residual BED: {outfile}"))
  }
}

save_peer_input <- function(normalize_bed_df, cov, output_dir, prefix){
  expr_file <- glue::glue("{output_dir}/{prefix}_normalized.tab")
  cov_file  <- glue::glue("{output_dir}/{prefix}_cov.tab")
  expression_matrix <- normalize_bed_df[, 5:ncol(normalize_bed_df)]
  transposed_matrix <- t(expression_matrix)
  write.table(transposed_matrix,file = expr_file,sep = "\t",row.names = FALSE,col.names = FALSE,quote = FALSE)
  write.table(cov,file = cov_file,sep = "\t",row.names = FALSE,col.names = FALSE,quote = FALSE)
  log_msg("INFO", glue("Generated PEER internal inputs -> Expr: {expr_file} | Cov: {cov_file}"))
  return(list(expr_file = expr_file,cov_file = cov_file))
}

peer_get_save_residuals <- function(expr_normalized_matrix,cov_matrix,num_factors,out_dir,prefix,normalize_bed_df,cov,out_peer_residual_bed,iteration_time){
  create_dir <- paste('mkdir','-p',out_dir)
  system(create_dir)
  
  peer <- paste(
      "peertool",
      "-f", expr_normalized_matrix,
      "-c", cov_matrix,
      "--add_mean",
      "-n", num_factors,
      "-o", out_dir,
      '-i',iteration_time
  )
  log_msg("INFO", glue("Executing CLI command: {peer}"))
  system(peer)
  
  log_msg("INFO", "PEER computation finished. Organising and renaming output files...")
  change_residuals <- paste('mv',glue("{out_dir}/residuals.csv"),glue("{out_dir}/{prefix}_peer_residuals.csv"))
  change_covs <- paste('mv',glue("{out_dir}/X.csv"),glue("{out_dir}/{prefix}_peer_factor.csv"))
  change_w<- paste('mv',glue("{out_dir}/W.csv"),glue("{out_dir}/{prefix}_peer_W.csv"))
  change_Alpha <- paste('mv',glue("{out_dir}/Alpha.csv"),glue("{out_dir}/{prefix}_peer_Alpha.csv"))
  system(change_residuals)
  system(change_covs)
  system(change_w)
  system(change_Alpha)

  pheno_col=normalize_bed_df[1:4]
  sample_col <- data.frame(sample = rownames(cov), stringsAsFactors = FALSE)
  numeric_matrix <- read.table(glue("{out_dir}/{prefix}_peer_residuals.csv"), header = FALSE,sep = ",")
  numeric_matrix <- as.matrix(numeric_matrix)
  counts_df <- as.data.frame(numeric_matrix)
  colnames(counts_df) <- sample_col$sample
  bed_expr_df <- cbind(pheno_col, counts_df)
  write.table(bed_expr_df, file=out_peer_residual_bed, sep="\t", row.names=FALSE, quote=FALSE)
  log_msg("SUCCESS", glue("Exported PEER Residual BED: {out_peer_residual_bed}"))
}

peer_for_HFI <- function(normalize_bed_df,cov,output_dir,prefix,out_peer_residual_bed,iteration_time, num_factors){
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  peer_input=save_peer_input(normalize_bed_df,cov,output_dir,prefix)
  expr_normalized_matrix <- peer_input$expr_file
  cov_matrix <- peer_input$cov_file
  peer_get_save_residuals(expr_normalized_matrix,cov_matrix,num_factors,output_dir,prefix,normalize_bed_df,cov,out_peer_residual_bed,iteration_time)
}

pca_for_HFI <- function(normalize_bed_df,cov,out_dir,prefix,tier) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  residual_gene_count_list <- PCA_getResiduals(normalize_bed_df,cov,tier)
  saveResiduals(residual_gene_count_list,out_dir,prefix)
}

main <- function(bed_path,cov_path,normalization_method,residual_methods,tier,output_dir,prefix,num_factors,iteration_time){
  start_time <- Sys.time()
  log_msg("INFO", "==== Pipeline Core Execution Started ====")
  
  log_msg("INFO", "Loading input files into memory...")
  cov <- read.delim(cov_path, header = TRUE,sep = "\t") 
  bed_df=read.delim(bed_path,sep='\t',header = TRUE)
  
  bed_sample_cols <- names(bed_df)[5:ncol(bed_df)]
  cov_samples <- cov$sample
  common_samples <- intersect(bed_sample_cols, cov_samples)
  
  log_msg("INFO", glue("Sample Alignment Summary:"))
  log_msg("INFO", glue("  - Samples in BED: {length(bed_sample_cols)}"))
  log_msg("INFO", glue("  - Samples in Covariates: {length(cov_samples)}"))
  log_msg("INFO", glue("  - Mutually Intersected Samples: {length(common_samples)}"))
  
  if (length(common_samples) == 0) {
    log_msg("ERROR", "No overlapping samples found between BED and Covariates file! Terminal exit.")
    stop("Error: Zero overlapping samples.")
  }
  
  # 数据子集过滤与对齐
  bed_df <- bed_df[, c(names(bed_df)[1:4], common_samples)]
  cov <- cov[cov$sample %in% common_samples, ]
  cov <- cov[match(common_samples, cov$sample), ]
  
  # 转换离散变量
  cov <- convert_to_dummy(cov)
  rownames(cov) <- cov[, 1]
  cov <- cov[, -1, drop = FALSE]
  
  ### 1. 标准化步骤
  log_msg("INFO", "Step 1: Running gene expression normalization...")
  if(normalization_method=='ratio_norma'){
      normalize_bed_df=ratio_norma(bed_df)
  } else if(normalization_method=='raw_count_norma'){
      normalize_bed_df=raw_count_norma(bed_df)
  }
  log_msg("SUCCESS", "Normalization phase completed.")

  ### 2. 残差计算步骤
  log_msg("INFO", "Step 2: Calculating residuals for hidden factors...")
  if (residual_methods=='pca'){
      log_msg("INFO", "Triggering PCA mode workflow...")
      pca_for_HFI(normalize_bed_df, cov, output_dir, prefix, tier)
  }

  if (residual_methods=='peer'){
      log_msg("INFO", "Triggering PEER mode workflow...")
      out_peer_residual_bed <- glue("{output_dir}/{prefix}_peer_residuals.bed")
      peer_for_HFI(normalize_bed_df, cov, output_dir, prefix, out_peer_residual_bed, iteration_time, num_factors)
  }
  
  end_time <- Sys.time()
  run_duration <- round(as.numeric(difftime(end_time, start_time, units = "secs")), 2)
  log_msg("SUCCESS", glue("==== Pipeline Finished Successfully! Total Elapsed Time: {run_duration} seconds ===="))
}

##################################################### 运行控制区 #############################################################

main(
  bed_path = opt$bed, 
  cov_path = opt$cov, 
  normalization_method = opt$normalization,
  residual_methods = opt$method,
  tier = tier,
  output_dir = opt$output,
  prefix = opt$prefix,
  num_factors = peer_factors_val,
  iteration_time = peer_iter_val
)


### usage ###########

### 1. PCA inference hidden factors
## non-tier ###
# Rscript your_script_name.R \
#   --bed ./data/expression.bed \
#   --cov ./data/covariates.txt \
#   --normalization raw_count_norma \
#   --method pca \
#   --output ./results_pca \
#   --prefix MyStudy \
#   --threads 4

# ## tier ########
# Rscript your_script_name.R \
#   --bed ./data/expression.bed \
#   --cov ./data/covariates.txt \
#   --normalization raw_count_norma \
#   --method pca \
#   --tier "-5,0,5" \
#   --output ./results_pca_tier \
#   --prefix MyStudy \
#   --threads 4

# ### 2. peer ###########
# Rscript your_script_name.R \
#   --bed ./data/expression.bed \
#   --cov ./data/covariates.txt \
#   --normalization ratio_norma \
#   --method peer \
#   --peer_factors 20 \  ## adjust according to GTEx
#   --peer_iter 1000 \
#   --output ./results_peer \
#   --prefix MyStudy \
#   --threads 2