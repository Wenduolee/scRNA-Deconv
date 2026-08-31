# scRNA-seq-Deconv

A specialized bioinformatic framework for donor-reference-free deconvolution of same-tissue biological mixtures using single-cell transcriptomics.

## Pipeline Overview
1. **Pseudo-bulk Variant Pre-filtering:** `FreeBayes` aggregates reads to discover high-confidence candidate SNPs.
2. **Single-cell Genotyping:** `CellSNP-lite` calculates per-cell allele depths at candidate loci.
3. **Variational Bayesian Deconvolution:** `Vireo` separates mixed cells into individual donor clusters.
4. **Forensic Statistical Evaluation:** PLINK LD pruning, Likelihood Ratio (LR) calculation, and biogeographic ancestry inference.
