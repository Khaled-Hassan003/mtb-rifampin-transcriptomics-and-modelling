#!/usr/bin/env Rscript
# 02_functional_enrichment.R
# ORA and GSEA Functional Enrichment Analysis via clusterProfiler

library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(dplyr)

# Load custom annotation mapping for M. tuberculosis H37Rv
go_map <- read.csv("go_map_annotation.csv") 
term2gene <- go_map %>% dplyr::select(GO_TERM, GENE_ID)

# Example ORA Enrichment for Signifcant DEGs
sig_genes <- read.csv("significant_degs.csv")$gene_id

go_enrich <- enricher(gene = sig_genes,
                      TERM2GENE = term2gene,
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.05)

# Plot Dotplot and save
png("01_transcriptomics/plots/go_enrichment_dotplot.png", width=900, height=600)
dotplot(go_enrich, showCategory=20) + ggtitle("Functional GO Enrichment under Rifampin Stress")
dev.off()
