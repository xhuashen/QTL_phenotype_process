# Pipeline for Molecular Phenotype Processing and Hidden Factor Correction in QTL Mapping

An efficient, fully automated R pipeline designed to preprocess molecular phenotypes (e.g., gene expression, exon splicing ratios) and perform confounding factor correction prior to Quantitative Trait Loci (QTL) mapping. 

By taking a standard phenotype BED file and a known covariates file, this pipeline automatically handles sample intersection alignment, robust normalization, and calculates final phenotype residuals. It regresses out both known biological/technical covariates and inferred hidden confounding factors from the phenotype matrix, outputting clean BED files ready for downstream QTL tools like MatrixEQTL, FastQTL, or tensorQTL.

---

## 🚀 Key Features

* **Dual Normalization Tracks:** * Supports ratio-based quantile normalization (`ratio_norma`) for normalized data types.
  * Supports raw count data processing using edgeR TMM + Log2 CPM (`raw_count_norma`).
  * All processing tracks conclude with a Rank-based Inverse Normal Transformation (RINT) to ensure phenotypes strictly follow a standard normal distribution.
* **Two Hidden Factor Inference Engines:**
  * **PCA Engine:** Integrated with `PCAForQTL` to automatically estimate the optimal number of Principal Components (PCs) via the Buja & Eyuboglu (BE) permutation or Elbow algorithm. Supports **multi-tier residual calculation** (e.g., testing relative variations like $-5$ or $+5$ PCs around the optimum) to maximize your cis-QTL discovery power.
  * **PEER Engine:** Seamlessly interfaces with the classic Probabilistic Estimation of Expression Residuals (PEER) algorithm using the high-performance `peertool` command-line utility.
* **Smart Covariate Handling:** Automatically detects categorical/discrete variables in your covariate file, expands them into dummy variables, and filters out known covariates that exhibit high collinearity with inferred PCs.
* **Zero Manual Alignment:** Automatically intersects and matches sample IDs between the phenotype BED headers and the covariates matrix, eliminating formatting mismatch errors.

---

## 📦 Prerequisites & Installation

### 1. Create a Conda Environment (Recommended)
We highly recommend isolating your dependencies using a dedicated environment:
```bash
conda create -n QTL_env r-base=4.3 -y
conda activate QTL_env
2. Install Core R DependenciesRun the following commands inside your R session to install CRAN and Bioconductor packages:Rinstall.packages(c("dplyr", "glue", "optparse"))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("limma", "edgeR", "DESeq2", "BiocParallel"))
3. Install Advanced Inference ToolsPCAForQTL (Required for automated PCA selection):Rdevtools::install_github("heatherjzhou/PCAForQTL")
PEER Tool (Required only if using --method peer):Please follow the compilation guidelines on the PEER GitHub Repository. Ensure that the compiled peertool binary executable is added to your system's $PATH environment variable.📂 Input File Formats1. Phenotype BED File (--bed)A tab-delimited text file where the first 4 columns contain genomic metadata, followed by individual sample quantification columns.Plaintext#chrom  start    end      gene_id   Sample1  Sample2  Sample3
chr1    10000    10500    GENE001   12.5     14.2     9.1
chr1    20000    21000    GENE002   104.0    95.2     110.1
2. Covariates File (--cov)A tab-delimited text file featuring a mandatory header column named sample (values must match the column names in your BED file), followed by known covariates (numeric or discrete).Plaintextsample   Age  Gender  Batch
Sample1  45   M       Batch_A
Sample2  52   F       Batch_B
Sample3  38   M       Batch_A
💻 Usage ExamplesClone this repository and navigate to the directory:Bashgit clone [https://github.com/xhuashen/QTL_phenotype_process.git](https://github.com/xhuashen/QTL_phenotype_process.git)
cd QTL_phenotype_process
Scenario 1: PEER Inferred Latent Factors(We recommend configuring --peer_factors based on sample size thresholds defined in the GTEx official guidelines.)BashRscript QTL_get_residual_test.R \
  --bed ./data/expression.bed \
  --cov ./data/covariates.txt \
  --normalization raw_count_norma \
  --method peer \
  --peer_factors 15 \
  --peer_iter 1000 \
  --output ./results_peer \
  --prefix MyStudy \
  --threads 4
Scenario 2: PCA with Automated Factor Selection (BE Algorithm)The pipeline will automatically determine the optimal PC threshold and directly output the corresponding residual dataset.BashRscript QTL_get_residual_test.R \
  --bed ./data/expression.bed \
  --cov ./data/covariates.txt \
  --normalization raw_count_norma \
  --method pca \
  --output ./results_pca \
  --prefix MyStudy \
  --threads 4
Scenario 3: PCA Multi-Tier StrategyProvide a relative vector via --tier to run and export multiple custom PC steps in parallel, allowing you to benchmark downstream cis-QTL discovery rates.BashRscript QTL_get_residual_test.R \
  --bed ./data/expression.bed \
  --cov ./data/covariates.txt \
  --normalization ratio_norma \
  --method pca \
  --tier "-5,0,5" \
  --output ./results_pca_tier \
  --prefix MyStudy \
  --threads 4
🛠️ Command-Line ArgumentsArgumentData TypeRequiredDefaultDescription--bedcharacterYes-Path to the standard tab-delimited phenotype BED file.--covcharacterYes-Path to the known biological/technical covariates file.--normalizationcharacterYes-Normalization strategy: ratio_norma or raw_count_norma.--methodcharacterNopcaConfounding factor inference method: pca or peer.--outputcharacterYes-Output destination directory path.--prefixcharacterYes-Prefix tag assigned to all generated output files.--threadsintegerNo1Core count for parallel execution (BiocParallel).--tiercharacterNoNULL[PCA Specific] Relative vector string for fine-tuning PC thresholds (e.g., -5,0,5).--peer_factorsintegerNoNULL[PEER Specific] Total number of hidden factors to estimate.--peer_iterintegerNo1000[PEER Specific] Maximum iterations allowed for PEER model convergence.📤 Output FilesUpon successful execution, the pipeline creates the following data files inside your designated --output directory:For --method pca:{prefix}_PCA_vectors.txt: The full calculated principal components (PC) matrix for all samples.{prefix}_tier_PC*.bed: The final corrected phenotype residual BED files. If you specified --tier "-5,0,5", three independent files corresponding to those configurations will be written (e.g., _tier_PC10.bed, _tier_PC15.bed, _tier_PC20.bed) for parallel QTL mapping tests.For --method peer:{prefix}_peer_residuals.bed: The corrected phenotype residual BED file ready for downstream QTL analysis.{prefix}_peer_factor.csv: Inferred latent sample confounding factors matrix ($X$).{prefix}_peer_W.csv: Genomewide weights/loadings matrix ($W$).{prefix}_peer_Alpha.csv: Factor precision parameters ($\alpha$). Smaller values denote hidden factors that explain a larger proportion of phenotype variation.
