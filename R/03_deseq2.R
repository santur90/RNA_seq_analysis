#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) stop("Usage: 03_deseq2.R <counts> <metadata> <outdir> <case> <control>")
suppressPackageStartupMessages({library(DESeq2); library(ggplot2)})
counts_file <- args[[1]]; metadata_file <- args[[2]]; outdir <- args[[3]]; case_label <- args[[4]]; control_label <- args[[5]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
counts_raw <- read.delim(counts_file, comment.char = "#", check.names = FALSE, row.names = 1)
counts <- counts_raw[, grep("Aligned.sortedByCoord.out.bam$", colnames(counts_raw)), drop = FALSE]
colnames(counts) <- sub("_Aligned.sortedByCoord.out.bam$", "", basename(colnames(counts)))
metadata <- read.delim(metadata_file, check.names = FALSE, row.names = 1)
metadata <- metadata[colnames(counts), , drop = FALSE]
stopifnot(all(rownames(metadata) == colnames(counts)))
metadata$condition <- factor(metadata$condition)
dds <- DESeqDataSetFromMatrix(round(as.matrix(counts)), metadata, design = ~ batch + condition)
dds <- dds[rowSums(counts(dds) >= 10) >= 2, ]
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", case_label, control_label), alpha = 0.05)
all_results <- as.data.frame(res); all_results$gene <- rownames(all_results)
sig <- subset(all_results, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) >= log2(1.5))
write.csv(all_results, file.path(outdir, "DESeq2_all_results.csv"), row.names = FALSE)
write.csv(sig, file.path(outdir, "DESeq2_significant_genes.csv"), row.names = FALSE)
write.csv(counts(dds, normalized = TRUE), file.path(outdir, "DESeq2_normalized_counts.csv"))
pdf(file.path(outdir, "DESeq2_QC_plots.pdf"), width = 9, height = 7)
plotMA(res, ylim = c(-5, 5)); plotDispEsts(dds); print(plotPCA(vst(dds, blind = FALSE), intgroup = c("condition", "batch"))); dev.off()
