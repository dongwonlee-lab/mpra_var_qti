# MPRAnalyze Workflow for Allelic Comparison

# Allelic comparison analysis was parallelized for efficient computation.
# In this script, 10 variants are processed at a time, specificed by job_num

# get input parameters to script, assigning pool number, version number, and job number
args <- commandArgs(trailingOnly=TRUE)
dna_annot_fn <- args[1]
rna_annot_fn <- args[2]
dna_counts_fn <- args[3]
rna_counts_fn <- args[4]
pool_num <- as.numeric(args[5])
job_num <- as.numeric(args[6])

library(MPRAnalyze)
library(tidyverse)
library(ggplot2)

# read in DNA allelic annotations
dna_annot_allelic <- read.table(dna_annot_fn, header=TRUE, colClasses="factor")

dna_annot_allelic$alleletype <- factor(dna_annot_allelic$alleletype, levels = c("ref", "alt"))

# read in RNA allelic annotations
rna_annot_allelic <- read.table(rna_annot_fn, header=TRUE, colClasses="factor")

rna_annot_allelic$alleletype <- factor(rna_annot_allelic$alleletype, levels = c("ref", "alt"))

# read in DNA allelic counts
dna_counts_allelic <- read.table(dna_counts_fn, header=TRUE)

# read in RNA allelic counts
rna_counts_allelic <- read.table(rna_counts_fn, header=TRUE)

# run with every 10 lines of data
if (job_num < 51) {
    slice <- c(((job_num*10)-9):(job_num*10))
} else if (job_num == 51) {
    slice <- c(501:509)
}

# make MPRAObject
expr_obj_allelic <- MpraObject(dnaCounts = as.matrix(dna_counts_allelic[slice,]),
                               rnaCounts = as.matrix(rna_counts_allelic[slice,]), 
                               dnaAnnot = dna_annot_allelic,
                               rnaAnnot = rna_annot_allelic)

print("Created MPRAObject")

# perform library size normalization separately for DNA and RNA counts
expr_obj_allelic <- estimateDepthFactors(expr_obj_allelic, lib.factor = c("batch"),
                            which.lib = "dna", 
                            depth.estimator = "uq")
expr_obj_allelic <- estimateDepthFactors(expr_obj_allelic, lib.factor = c("batch"),
                            which.lib = "rna", 
                            depth.estimator = "uq")

# allelic comparison analysis
expr_obj_allelic <- analyzeComparative(obj = expr_obj_allelic, 
                           dnaDesign = ~ alleletype + barcode_allelic_version + version, 
                           rnaDesign = ~ alleletype + version, 
                           reducedDesign = ~ version)

print("Set up design matrix for Allelic Comparison")

# save fitted model
model_fn <- paste("../results/pool",
                pool_num,
                "/MPRAObject_allelic_", 
                slice[1], '-', tail(slice, n=1),
                ".rds",
                sep="")
saveRDS(expr_obj_allelic, file = model_fn)

# extract alpha values (transcription rates) from the fitted model
expr_alpha_allelic <- getAlpha(expr_obj_allelic)

# run likelihood ratio test
res <- testLrt(expr_obj_allelic)

# save likelihood ratio test results
lrt_res_fn <- paste("../results/pool",
                pool_num,
                "/lrt_results_", 
                slice[1], '-', tail(slice, n=1),
                ".txt",
                sep="")
write.table(res, lrt_res_fn)

# split row names into enhancer information
names_lrt <- data.frame(do.call(rbind, strsplit(rownames(res), "_")))
colnames(names_lrt) <- c("chrom", "pos", "snp", "allele", "mutation")

# create new dataframe with snp, ref/alt, and alpha values
results_stat <- data.frame('snp'= names_lrt$snp, 'logFC'=res$logFC, lrt_pval=res$pval, 'lrt_nlog10_pval'=-log10(res$pval))

# save MPRAnalyze allelic comparison results to text file
res_fn <- paste("../results/pool",
                pool_num,
                "/pool",
                pool_num,
                "_mpranalyze_allelic_results_", 
                slice[1], '-', tail(slice, n=1),
                ".txt",
                sep="")
write.table(results_stat, res_fn)
