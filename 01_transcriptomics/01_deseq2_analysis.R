#!/usr/bin/env Rscript
# 01_deseq2_analysis.R
# RNA-Seq Differential Expression Analysis for M. tuberculosis H37Rv

library(DESeq2)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(EnhancedVolcano)

# 1. Load Data
counts_mat <- read.csv("counts_mat.csv.gz", row.names = 1)
counts_mat <- as.matrix(counts_mat)

# 2. Setup Phenotype MetaData
coldata <- data.frame(
  sample_id = colnames(counts_mat),
  condition = c(rep("pstA1_ctrl", 3), rep("pstA1_RIF_30min", 3), rep("pstA1_RIF_60min", 3),
                rep("WT_ctrl", 3), rep("WT_RIF_30min", 3), rep("WT_RIF_60min", 3)),
  row.names = colnames(counts_mat)
)
coldata$condition <- as.factor(coldata$condition)

# 3. Quality Control & DESeq2 Design
dds <- DESeqDataSetFromMatrix(countData = counts_mat, colData = coldata, design = ~ condition)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# 4. Run Differential Expression
dds <- DESeq(dds)

# 5. Extract Results & Export Plots
res_ctrl <- results(dds, contrast = c("condition", "pstA1_ctrl", "WT_ctrl"))
res_30min <- results(dds, contrast = c("condition", "pstA1_RIF_30min", "WT_RIF_30min"))
res_60min <- results(dds, contrast = c("condition", "pstA1_RIF_60min", "WT_RIF_60min"))

# Save Volcano Plot Example
png("01_transcriptomics/plots/volcano_30min.png", width=800, height=800)
EnhancedVolcano(res_30min, lab = rownames(res_30min), x = 'log2FoldChange', y = 'padj', pCutoff = 0.05, FCcutoff = 0.5)
dev.off()
