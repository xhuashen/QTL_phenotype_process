# Pipeline for Molecular Phenotype Processing and Hidden Factor Correction in QTL Mapping

An efficient, fully automated R pipeline designed to preprocess molecular phenotypes (e.g., gene expression, exon splicing ratios) and perform confounding factor correction prior to Quantitative Trait Loci (QTL) mapping. 

By taking a standard phenotype BED file and a known covariates file, this pipeline automatically handles sample intersection alignment, robust normalization, and calculates final phenotype residuals. It regresses out both known biological/technical covariates and inferred hidden confounding factors from the phenotype matrix, outputting clean BED files ready for downstream QTL tools like MatrixEQTL, FastQTL, or tensorQTL.

---

## 🚀 Key Features

* **Dual Normalization Tracks:**
  * Supports ratio-based quantile normalization (`ratio_norma`) for normalized data types.
  * Supports raw count data processing using edgeR TMM + Log2 CPM (`raw_count_norma`).
  * All processing tracks conclude with a Rank-based Inverse Normal Transformation (RINT) to ensure phenotypes strictly follow a standard normal distribution.
* **Two Hidden Factor Inference Engines:**
  * **PCA Engine:** Integrated with `PCAForQTL` to automatically estimate the optimal number of Principal Components (PCs) via the Buja & Eyuboglu (BE) permutation or Elbow algorithm. Supports **multi-tier residual calculation** (e.g., testing relative variations like $-5$ or $+5$ PCs around the optimum) to maximize your cis-QTL discovery power.
  * **PEER Engine:** Seamlessly interfaces with the classic Probabilistic Estimation of Expression Residuals (PEER) algorithm using the high-performance `peertool` command-line utility.
* **Smart Covariate Handling:** Automatically detects categorical/discrete variables in your covariate file, expands them into dummy variables, and filters out known covariates that exhibit high collinearity with inferred PCs.
* **Zero Manual Alignment:** Automatically intersects and matches sample IDs between the phenotype BED headers and the covariates matrix, eliminating formatting mismatch errors.

---

## 📦 Install and Usage

we recommend creat a new conda environment
```bash
conda create -n QTL_env

then you need to install the PCAForQTL from https://github.com/heatherjzhou/PCAForQTL and command tools of peer (peertool) from https://github.com/PMBio/peer

📂 Input File Formats
1. Phenotype BED File (--bed)
A tab-delimited text file where the first 4 columns contain genomic metadata, followed by individual sample quantification columns.
#chrom  start    end      gene_id   Sample1  Sample2  Sample3
chr1    10000    10500    GENE001   12.5     14.2     9.1
chr1    20000    21000    GENE002   104.0    95.2     110.1

2. Covariates File (--cov)
A tab-delimited text file featuring a mandatory header column named samplefollowed by known covariates (numeric or discrete,this script will convert the discrete to dummy).
sample   Age  Gender  Batch
Sample1  45   M       Batch_A
Sample2  52   F       Batch_B
Sample3  38   M       Batch_A

Note, ensuring the samples in Phenotype BED File are a subset of Covariates File.

💻 Usage
git clone [https://github.com/xhuashen/QTL_phenotype_process.git](https://github.com/xhuashen/QTL_phenotype_process.git)
cd QTL_phenotype_process

1. peer model
Rscript QTL_get_residual_test.R \
  --bed ./data/expression.bed \
  --cov ./data/covariates.txt \
  --normalization raw_count_norma \
  --method peer \
  --peer_factors 15 \
  --peer_iter 1000 \
  --output ./results_peer \
  --prefix MyStudy \
  --threads 4

2. pca BE methods define model
Rscript QTL_get_residual_test.R \
  --bed ./data/expression.bed \
  --cov ./data/covariates.txt \
  --normalization raw_count_norma \
  --method pca \
  --output ./results_pca \
  --prefix MyStudy \
  --threads 4

3. pca multi tier model
Rscript QTL_get_residual_test.R \
  --bed ./data/expression.bed \
  --cov ./data/covariates.txt \
  --normalization raw_count_norma \
  --method pca \
  --tier=-5,0,5 \
  --output ./results_pca_tier \
  --prefix MyStudy \
  --threads 4

🛠️ Parameters
--bed (character | Required)
Path to the standard tab-delimited phenotype BED file.

--cov (character | Required)
Path to the known biological/technical covariates file.

--normalization (character | Required)
Normalization strategy to use: ratio_norma or raw_count_norma.

--method (character | Default: pca)
Confounding factor inference method: pca or peer.

--output (character | Required)
Output destination directory path.

--prefix (character | Required)
Prefix tag assigned to all generated output files.

--threads (integer | Default: 1)
Core count for parallel execution (BiocParallel).

--tier (character | Default: NULL)
[PCA Specific] Relative vector string for fine-tuning PC thresholds (e.g., -5,0,5).

--peer_factors (integer | Default: NULL)
[PEER Specific] Total number of hidden factors to estimate.

--peer_iter (integer | Default: 1000)
[PEER Specific] Maximum iterations allowed for PEER model convergence.

📤 Output Files
If you use the PCA method, the pipeline generates:

prefix_PCA_vectors.txt : The full calculated principal components matrix for all samples.

prefix_tier_PC_num.bed : The final corrected phenotype residual BED files. If you used tier -5,0,5, multiple files will be generated for parallel testing.

If you use the PEER method, the pipeline generates:

prefix_peer_residuals.bed : The corrected phenotype residual BED file ready for downstream QTL analysis.

prefix_peer_factor.csv : Inferred latent sample confounding factors matrix.

prefix_peer_W.csv : Genomewide weights and loadings matrix.

prefix_peer_Alpha.csv : Factor precision parameters where smaller values denote higher importance.



