## Analysis of bulk RNA-seq of E11.5 Pbx1/2 mutant tissue of lambdoidal junction
# Author: Vera Laub, adapted from Martas analysis code
# Date last modified: 2025-09-23

# Set working directory (if not already done)
setwd("/Users/veralaub/Documents/postdoc/bioinformatics/data/RNA-seq/midface/E11.5_Pbx1+2_mutants_vs_ctrl/Louis_bulk_analyses/")

# Read the RDS file
adata_deseq2 <- readRDS("adata_deseq2 (1).rds")

# Check the object type and content
class(adata_deseq2)
str(adata_deseq2)

# -------------------------------
# Violin plots for genes of interest
# -------------------------------

# Load libraries
library(DESeq2)
library(ggplot2)
library(tidyr)
library(org.Mm.eg.db)

# Set working directory (if not already done)
setwd("/Users/veralaub/Documents/postdoc/bioinformatics/data/RNA-seq/midface/E11.5_Pbx1+2_mutants_vs_ctrl/2025-09-23_Vera_analysis/")

# Extract normalized counts
norm_counts <- counts(adata_deseq2, normalized = TRUE)

# Map Ensembl IDs to gene symbols
symbols <- mapIds(org.Mm.eg.db,
                  keys = rownames(norm_counts),
                  column = "SYMBOL",
                  keytype = "ENSEMBL",
                  multiVals = "first")

# Add as a column
df <- as.data.frame(norm_counts)
df$gene_symbol <- symbols[rownames(df)]

# Reshape to long format
df_long <- pivot_longer(df, -gene_symbol, names_to = "sample", values_to = "counts")

# Add condition info from colData
df_long$condition <- colData(adata_deseq2)$condition[match(df_long$sample, colnames(adata_deseq2))]

# Genes of interest
genes_of_interest <- c("Pbx1", "Pbx2", "Zfhx3", "Bmp4", "Bmpr1a", 
                       "Smad4", "Id2", "Itga3", "Itga6")

# Subset only those genes (remove NAs)
df_plot <- subset(df_long, gene_symbol %in% genes_of_interest & !is.na(gene_symbol))

# Now plot
p <- ggplot(df_plot, aes(x = condition, y = counts, fill = condition)) +
  geom_violin(trim = FALSE) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 1) +
  facet_wrap(~ gene_symbol, scales = "free_y") +
  theme_bw() +
  labs(title = "Expression of selected genes",
       y = "Normalized counts", x = "") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5))

# Save
ggsave("violin_genes_of_interest.png", p, width = 12, height = 8, dpi = 300)
ggsave("violin_genes_of_interest.pdf", p, width = 12, height = 8)

print(p)


################################################################################
### CODE USING DATA FROM MARTA PROBABLY UNSUCCESSFUL EXPERIMENT BC PBX1/2 KO BAD


#### Analysis of bulk RNA-seq of FACS-sorted (Epcam-APC) E11.5 epithelial cells 
#### from Pbx1/2 mutant embryos
# Author: Vera Laub, adapted from Martas analysis code
# Date last modified: 2025-09-23

## Combine readcounts files in one table
##first column, gene name,
##One column per sample.
##*********all readcounts file must contain all genes and the same number of rows.

# Set wd to retrieve data
setwd("/Users/veralaub/Documents/postdoc/bioinformatics/data/RNA-seq/midface/EPI/analysis_only_Mut_vs_C")

#Read data from each replicate
#controls
c4epi <- read.table("C4EPI_read_counts.txt", header=F, row.names=1, quote="", sep="")
c6epi <- read.table("C6EPI_read_counts.txt", header=F, row.names=1, quote="", sep="")
c7epi <- read.table("C7EPI_read_counts.txt", header=F, row.names=1, quote="", sep="")
c8epi <- read.table("C8EPI_read_counts.txt", header=F, row.names=1, quote="", sep="")
c9epi <- read.table("C9EPI_read_counts.txt", header=F, row.names=1, quote="", sep="")

#mutants
m2epi <- read.table("M2EPI_read_counts.txt", header=F, row.names=1, quote="", sep="")
m6epi <- read.table("M6EPI_read_counts.txt", header=F, row.names=1, quote="", sep="")
m7epi <- read.table("M7EPI_read_counts.txt", header=F, row.names=1, quote="", sep="")
m8epi <- read.table("M8_EPI_read_counts.txt", header=F, row.names=1, quote="", sep="")
m9epi <- read.table("M9EPI_read_counts.txt", header=F, row.names=1, quote="", sep="")

#add data together
epifinal<- cbind(c4epi, c6epi, c7epi, c8epi, c9epi, m2epi, m6epi, m7epi, m8epi, m9epi)
head(epifinal)


#add label to top row
names(epifinal) <- c("c4epi", "c6epi", "c7epi", "c8epi", "c9epi", "m2epi", "m6epi", "m7epi", "m8epi", "m9epi")
head(epifinal)

# Set wd to place data
#set wd
setwd("/Users/veralaub/Documents/postdoc/bioinformatics/data/RNA-seq/midface/EPI/analysis_VL_only_Mut_vs_C")

#save file
write.table(epifinal, "epimerged_MvsC_read_counts_ENSEMBLE_IDs.txt", row.names=T, quote=F)
epifinal[555,]


### Add gene names to file

# Load library
library(biomaRt)

# Connect to the Ensembl mouse dataset
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# Clean Ensembl IDs if they have version suffixes like ".1", ".2"
rownames(epifinal) <- sub("\\..*", "", rownames(epifinal))

# Get mapping from Ensembl to gene symbol
mapping <- getBM(
  attributes = c("ensembl_gene_id", "mgi_symbol"),
  filters = "ensembl_gene_id",
  values = rownames(epifinal),
  mart = mart
)

# Create lookup vector: Ensembl ID → gene symbol
id2symbol <- setNames(mapping$mgi_symbol, mapping$ensembl_gene_id)

# Add gene symbol column to your count table
epifinal_annot <- data.frame(
  gene_symbol = id2symbol[rownames(epifinal)],
  ensembl_id  = rownames(epifinal),
  epifinal,
  row.names = NULL
)

# Save to file
write.table(
  epifinal_annot,
  "epimerged_MvsC_read_counts_with_gene_names.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

head(epifinal_annot)


###############################################################################
# RNA-seq downstream analysis (adapted for your epifinal dataset)
# - input: `epifinal` (raw counts), rownames = Ensembl IDs (may include version suffixes)
# - output: filtered counts, rlog-normalized matrix, DESeq2 results, PCA, MA plot, heatmaps
# Author: adapted for Vera Laub dataset (Sep 23 2025)
###############################################################################

# ----------------------------
# Install / load packages
# ----------------------------
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
pkgs_bioc <- c("DESeq2", "edgeR", "limma", "org.Mm.eg.db")
for (p in pkgs_bioc) if (!requireNamespace(p, quietly=TRUE)) BiocManager::install(p)

pkgs_cran <- c("ggplot2", "gdata", "gplots", "Hmisc")
for (p in pkgs_cran) if (!requireNamespace(p, quietly=TRUE)) install.packages(p)

# load libraries
library(DESeq2)
library(edgeR)
library(limma)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(ggplot2)
library(gplots)
library(ggrepel)

# ----------------------------
# (A) Prepare data
# ----------------------------
# If you already have 'epifinal' loaded in the environment (from earlier steps) skip reading;
# otherwise uncomment and read the file produced earlier:
# epifinal <- read.table("epimerged_MvsC_read_counts_ENSEMBLE_IDs.txt", header=TRUE, row.names=1, sep="\t", quote="", check.names=FALSE)

# Ensure rownames are present
if (is.null(rownames(epifinal))) stop("epifinal must have rownames with Ensembl IDs")

# strip version suffix (ENSMUSG00000012345.2 -> ENSMUSG00000012345)
rownames(epifinal) <- sub("\\..*$", "", rownames(epifinal))

# add gene symbol column using org.Mm.eg.db (keep Ensembl IDs as unique rownames)
symbols <- mapIds(org.Mm.eg.db,
                  keys = rownames(epifinal),
                  column = "SYMBOL",
                  keytype = "ENSEMBL",
                  multiVals = "first")
# symbols will have NA where no symbol found

# create a data.frame with gene_symbol and ensembl_id as columns and counts after
epifinal_annot <- data.frame(gene_symbol = symbols[rownames(epifinal)],
                             ensembl_id = rownames(epifinal),
                             epifinal,
                             row.names = NULL,
                             check.names = FALSE)

# Save annotated counts (both columns included)
write.table(epifinal_annot,
            file = "epimerged_MvsC_read_counts_with_symbols.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

# For DESeq2 we will keep raw counts as a matrix with Ensembl IDs as rownames:
counts <- as.matrix(epifinal)
rownames(counts) <- epifinal_annot$ensembl_id

# ----------------------------
# (B) Build coldata (sample metadata)
# ----------------------------
# Based on your sample variable names set earlier: c4epi,c6epi,c7epi,c8epi,c9epi,m2epi,m6epi,m7epi,m8epi,m9epi
sample_names <- colnames(counts)

# Create Condition column: controls (c*) vs mutants (m*)
Condition <- ifelse(grepl("^c", sample_names, ignore.case = TRUE), "Control",
                    ifelse(grepl("^m", sample_names, ignore.case = TRUE), "Mutant", NA))
if (any(is.na(Condition))) stop("Some sample names not matching expected ^c or ^m pattern: ", paste(sample_names[is.na(Condition)], collapse=", "))

coldata <- data.frame(SampleName = sample_names,
                      Condition = factor(Condition, levels = c("Control","Mutant")),
                      row.names = sample_names,
                      stringsAsFactors = FALSE)

# Save coldata
write.table(coldata, file = "coldata_MvsC.tsv", sep = "\t", quote = FALSE, row.names = TRUE)

# ----------------------------
# (C) Filtering (CPM-based) and save filtered raw counts
# ----------------------------
# Use edgeR's cpm to filter: keep genes with CPM > 0 in at least 3 samples (same threshold as tutorial)
cpm_mat <- cpm(counts)
keep <- rowSums(cpm_mat > 0) >= 3
counts_filt <- counts[keep, ]
cat("Genes before filtering:", nrow(counts), "after filtering:", nrow(counts_filt), "\n")
write.table(counts_filt, file = "featurecounts-genes-cpm0-filtered-raw.txt", sep = "\t", quote = FALSE, row.names = TRUE)

# ----------------------------
# (D) DESeq2 normalization + rlog (for PCA/plots)
# ----------------------------
dds <- DESeqDataSetFromMatrix(countData = counts_filt, colData = coldata, design = ~ Condition)

# set reference level to Control so fold changes are Mutant vs Control
dds$Condition <- relevel(dds$Condition, ref = "Control")

# run DESeq (this performs size factor estimation, dispersion estimation, negative binomial fitting)
dds <- DESeq(dds)

# rlog transform for visualization (useful for PCA, heatmaps)
rld <- rlog(dds, blind = FALSE)
exprsR <- assay(rld)
write.table(exprsR, file = "DESEQ-genes-cpm0-RDN.txt", sep = "\t", quote = FALSE, row.names = TRUE)

# ----------------------------
# (E) Correlation matrix and heatmap of sample correlations
# ----------------------------
data_corr <- cor(exprsR, method = "spearman")
write.table(data_corr, file = "sample_correlation_matrix.tsv", sep = "\t", quote = FALSE, row.names = TRUE)

# plot correlation heatmap (saves to PNG)
png("sample_correlation_heatmap.png", width = 1000, height = 1000, res = 150)
par(oma = c(3,4,3,4))
heatmap.2(data_corr, scale="none",
          breaks=seq(min(data_corr, na.rm=TRUE), max(data_corr, na.rm=TRUE), length.out=21),
          col = colorpanel(20, "orange", "cornsilk"),
          key=TRUE, symkey=FALSE, density.info="none", trace="none", cexCol=1.2)
dev.off()

# ----------------------------
# (F) PCA plotting
# ----------------------------
pcaData <- plotPCA(rld, intgroup = "Condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pcaData$Sample <- rownames(pcaData)

png("PCA_plot_PC1_PC2.png", width = 1000, height = 800, res = 150)
ggplot(pcaData, aes(PC1, PC2, color = Condition, label = Sample)) +
  geom_point(size = 4) +
  geom_text(vjust = -1, size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw() +
  theme(text = element_text(size = 14))
dev.off()

# Also standard prcomp alternative (if you prefer):
pca2 <- prcomp(t(exprsR), center = TRUE, scale. = FALSE)
scores <- as.data.frame(pca2$x)
scores$Sample <- rownames(scores)
# quick plot (PC1, PC2)
png("PCA_prcomp_PC1_PC2.png", width = 1000, height = 800, res = 150)
ggplot(scores, aes(PC1, PC2)) + geom_point(aes(color = coldata[Sample, "Condition"]), size = 4) +
  geom_text(aes(label = Sample), vjust = -1, size = 3) + theme_bw()
dev.off()

# ----------------------------
# (G) Differential expression with DESeq2
# ----------------------------
res <- results(dds, alpha = 0.05)   # default contrast is Mutant vs Control because of relevel above
resOrdered <- res[order(res$log2FoldChange), ]
# add gene_symbol and ensembl id columns to results for easier reading
resDF <- as.data.frame(resOrdered)
resDF$ensembl_id <- rownames(resDF)
resDF$gene_symbol <- symbols[resDF$ensembl_id]

# Write results (full) and filtered lists
write.table(resDF, file = "DESeq2_all_results.tsv", sep = "\t", quote = FALSE, row.names = TRUE)

resSig <- subset(resDF, padj < 0.05)
write.table(resSig, file = "DESeq2_padj0.05_results.tsv", sep = "\t", quote = FALSE, row.names = TRUE)

resSigFC <- subset(resSig, abs(log2FoldChange) > 1)
write.table(resSigFC, file = "DESeq2_padj0.05_log2FC_gt1_results.tsv", sep = "\t", quote = FALSE, row.names = TRUE)

# ----------------------------
# (H) MA plot (with coloring and optional labels)
# ----------------------------
# attach plotting color column
resDF$plotColor <- "black"
resDF$plotColor[resDF$log2FoldChange > 1 & resDF$padj < 0.05] <- "red"
resDF$plotColor[resDF$log2FoldChange < -1 & resDF$padj < 0.05] <- "deepskyblue"

png("MA_plot_all_genes.png", width = 1200, height = 900, res = 150)
ggplot(data = resDF, aes(x = log10(baseMean + 1), y = log2FoldChange)) +
  geom_point(alpha = 0.5, size = 1, aes(color = plotColor)) +
  scale_color_identity() +
  geom_hline(yintercept = 0, colour = "black", size = 0.6) +
  geom_hline(yintercept = 1, colour = "red", size = 0.6) +
  geom_hline(yintercept = -1, colour = "blue", size = 0.6) +
  xlab("log10 mean expression") + ylab("log2 fold change (Mutant / Control)") +
  theme_bw()
dev.off()

# Optional: label a few genes of interest by symbol (only if they exist in results)
genes_of_interest <- c("Pbx1", "Pbx2", "Zfhx3", "Bmp4", "Bmpr1a", "Smad4", "Id2", "Itga3", "Itga6")  # change to desired symbols
selRows <- subset(resDF, gene_symbol %in% genes_of_interest)
if (nrow(selRows) > 0) {
  png("MA_plot_highlight_genes.png", width = 1200, height = 900, res = 150)
  print(
    ggplot(data = resDF, aes(x = log10(baseMean + 1), y = log2FoldChange)) +
      geom_point(alpha = 0.3, size = 1, aes(color = plotColor)) +
      scale_color_identity() +
      geom_point(data = selRows, color = "green", size = 3) +
      geom_text(data = selRows, aes(label = gene_symbol), vjust = -1, size = 4) +
      geom_hline(yintercept = 0) +
      theme_bw()
  )
  dev.off()
}


if (nrow(selRows) > 0) {
  # Create a dynamic nudge: positive log2FC → nudge up, negative log2FC → nudge down
  selRows$nudge_y <- ifelse(selRows$log2FoldChange >= 0, 1, -1)
  
  png("MA_plot_highlight_genes_Bmp4-signaling.png", width = 1200, height = 900, res = 150)
  print(
    ggplot(data = resDF, aes(x = log10(baseMean + 1), y = log2FoldChange)) +
      geom_point(alpha = 0.3, size = 1, aes(color = plotColor)) +
      scale_color_identity() +
      geom_point(data = selRows, color = "green", size = 3) +
      geom_text_repel(
        data = selRows,
        aes(label = gene_symbol),
        nudge_y = selRows$nudge_y,        # dynamic vertical nudging
        direction = "y",
        segment.color = 'grey50',
        size = 4
      ) +
      geom_hline(yintercept = 0) +
      theme_bw()
  )
  dev.off()
}


# ----------------------------
# (I) Heatmaps
# ----------------------------
# prepare distance function
pcorr <- function(x) as.dist(1 - cor(t(x), method = "pearson"))
my_palette2 <- colorRampPalette(c("blue4", "white", "red4"))(n = 299)

# heatmap for the whole rlog matrix (slow for large matrices)
png("heatmap_all_genes_exprsR.png", width = 1800, height = 1400, res = 150)
heatmap.2(as.matrix(exprsR),
          dendrogram = "both",
          col = my_palette2,
          scale = "row",
          key = TRUE, density.info = "none", trace = "none",
          cexCol = 1, Colv = TRUE, Rowv = TRUE,
          labRow = FALSE,
          distfun = pcorr, margins = c(6,8))
dev.off()

# heatmap for a small selectable gene list: you can use either Ensembl IDs or gene symbols
# Example vector of symbols you want to plot
sel_symbols <- c("Pbx1", "Pbx2", "Zfhx3", "Bmp4", "Bmpr1a", "Smad4", "Id2", "Itga3", "Itga6")

# Use your mapping table to look up Ensembl IDs
sel_ensembls <- mapping$ensembl_gene_id[mapping$mgi_symbol %in% sel_symbols]

# Keep only those Ensembl IDs that actually exist in exprsR
sel_ensembls <- sel_ensembls[sel_ensembls %in% rownames(exprsR)]

# map symbols to rownames (ensembl) using our symbol vector
sel_ensembls <- names(symbols)[symbols %in% sel_symbols]
if (length(sel_ensembls) > 0) {
  sel_mat <- exprsR[sel_ensembls, , drop = FALSE]
  rownames(sel_mat) <- symbols[rownames(sel_mat)]  # use gene symbols in plotted rows
  png("heatmap_selected_genes.png", width = 1000, height = 800, res = 150)
  heatmap.2(as.matrix(sel_mat), dendrogram = "both", col = my_palette2, scale = "row",
            margins = c(6,8), key = TRUE, density.info = "none", trace = "none",
            cexRow = 1, cexCol = 1, Colv = TRUE, Rowv = TRUE)
  dev.off()
}

# ----------------------------
# Volcano plot
# ----------------------------
library(ggplot2)
library(ggrepel)

# attach a color column based on significance thresholds
resDF$plotColor <- "grey"
resDF$plotColor[resDF$log2FoldChange > 1 & resDF$padj < 0.05] <- "red"
resDF$plotColor[resDF$log2FoldChange < -1 & resDF$padj < 0.05] <- "deepskyblue"

# define genes to highlight
genes_of_interest <- c("Pbx1", "Pbx2", "Zfhx3", "Bmp4", "Bmpr1a", "Smad4", "Id2", "Itga3", "Itga6")
selRows <- subset(resDF, gene_symbol %in% genes_of_interest)

# create a dynamic vertical nudge: positive log2FC → up, negative → down
if (nrow(selRows) > 0) {
  selRows$nudge_y <- ifelse(selRows$log2FoldChange >= 0, 1, -1)
}

# open PNG device
png("volcano_plot_highlight_genes.png", width = 1200, height = 900, res = 150)

# volcano plot
print(
  ggplot(resDF, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = plotColor), alpha = 0.6, size = 1.5) +
    scale_color_identity() +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
    xlab("Log2 Fold Change") + ylab("-log10(adj p-value)") +
    theme_minimal() +
    theme(legend.position = "bottom") +
    # highlight genes of interest
    geom_point(data = selRows, color = "green", size = 3) +
    geom_text_repel(
      data = selRows,
      aes(label = gene_symbol),
      nudge_y = selRows$nudge_y,
      direction = "y",
      segment.color = "grey50",
      size = 4
    )
)

dev.off()


# ----------------------------
# (J) Final messages
# ----------------------------
cat("Analysis complete.\n")
cat("Files written:\n")
cat("- epimerged_MvsC_read_counts_with_symbols.tsv\n")
cat("- featurecounts-genes-cpm0-filtered-raw.txt\n")
cat("- DESEQ-genes-cpm0-RDN.txt\n")
cat("- sample_correlation_matrix.tsv\n")
cat("- DESeq2_all_results.tsv\n")
cat("- DESeq2_padj0.05_results.tsv\n")
cat("- DESeq2_padj0.05_log2FC_gt1_results.tsv\n")
cat("- MA_plot_all_genes.png, MA_plot_highlight_genes.png (if genes found)\n")
cat("- heatmap images\n")

###############################################################################
# End of script
###############################################################################


#### Analysis of bulk RNA-seq of FACS-sorted (Epcam-APC) E11.5 epithelial cells 
#### from Pbx1/2 mutant embryos
#### Filtered dataset: only c4epi, c6epi, c7epi, c8epi, c9epi, m2epi, m6epi, m8epi
# Author: Vera Laub, adapted from Martas analysis code
# Date last modified: 2025-09-23

# ----------------------------
# (A) Combine readcounts files (filtered dataset)
# ----------------------------
setwd("/Users/veralaub/Documents/postdoc/bioinformatics/data/RNA-seq/midface/EPI/analysis_only_Mut_vs_C")

# Read data
c4epi <- read.table("C4EPI_read_counts.txt", header=F, row.names=1, quote="", sep="")
c8epi <- read.table("C8EPI_read_counts.txt", header=F, row.names=1, quote="", sep="")
c9epi <- read.table("C9EPI_read_counts.txt", header=F, row.names=1, quote="", sep="")

m2epi <- read.table("M2EPI_read_counts.txt", header=F, row.names=1, quote="", sep="")
m6epi <- read.table("M6EPI_read_counts.txt", header=F, row.names=1, quote="", sep="")
m8epi <- read.table("M8_EPI_read_counts.txt", header=F, row.names=1, quote="", sep="")

# Combine only the selected samples
epifinal <- cbind(c4epi, c8epi, c9epi, m2epi, m6epi, m8epi)

# Add column names
names(epifinal) <- c("c4epi", "c8epi", "c9epi", "m2epi", "m6epi", "m8epi")

# Set working directory for all outputs
setwd("/Users/veralaub/Documents/postdoc/bioinformatics/data/RNA-seq/midface/EPI/analysis_VL_only_Mut_vs_C")

# ----------------------------
# (B) Annotate gene symbols
# ----------------------------
library(biomaRt)

# Connect to Ensembl
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
rownames(epifinal) <- sub("\\..*", "", rownames(epifinal))

mapping <- getBM(
  attributes = c("ensembl_gene_id", "mgi_symbol"),
  filters = "ensembl_gene_id",
  values = rownames(epifinal),
  mart = mart
)

id2symbol <- setNames(mapping$mgi_symbol, mapping$ensembl_gene_id)
epifinal_annot <- data.frame(
  gene_symbol = id2symbol[rownames(epifinal)],
  ensembl_id  = rownames(epifinal),
  epifinal,
  row.names = NULL
)

write.table(epifinal_annot,
            "epimerged_MvsC_read_counts_with_gene_names_filtered.txt",
            sep="\t", quote=FALSE, row.names=FALSE)

# ----------------------------
# (C) DESeq2 analysis
# ----------------------------
library(DESeq2)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(ggplot2)
library(ggrepel)
library(gplots)
library(edgeR)
library(limma)

counts <- as.matrix(epifinal)
rownames(counts) <- epifinal_annot$ensembl_id

# Build coldata
sample_names <- colnames(counts)
Condition <- ifelse(grepl("^c", sample_names), "Control", "Mutant")
coldata <- data.frame(SampleName = sample_names,
                      Condition = factor(Condition, levels=c("Control","Mutant")),
                      row.names = sample_names)

# Filter genes with low counts
cpm_mat <- cpm(counts)
keep <- rowSums(cpm_mat > 0) >= 3
counts_filt <- counts[keep, ]

dds <- DESeqDataSetFromMatrix(countData = counts_filt, colData = coldata, design = ~ Condition)
dds$Condition <- relevel(dds$Condition, ref = "Control")
dds <- DESeq(dds)
rld <- rlog(dds, blind = FALSE)
exprsR <- assay(rld)

# ----------------------------
# (D) Differential expression
# ----------------------------
res <- results(dds, alpha=0.05)
resOrdered <- res[order(res$log2FoldChange), ]
symbols <- mapIds(org.Mm.eg.db,
                  keys = rownames(resOrdered),
                  column = "SYMBOL",
                  keytype = "ENSEMBL",
                  multiVals = "first")
resDF <- as.data.frame(resOrdered)
resDF$ensembl_id <- rownames(resDF)
resDF$gene_symbol <- symbols[resDF$ensembl_id]

# Save results
write.table(resDF, "DESeq2_all_results_filtered.tsv", sep="\t", quote=FALSE, row.names=TRUE)
write.table(subset(resDF, padj<0.05), "DESeq2_padj0.05_results_filtered.tsv", sep="\t", quote=FALSE, row.names=TRUE)
write.table(subset(resDF, padj<0.05 & abs(log2FoldChange)>1), "DESeq2_padj0.05_log2FC_gt1_results_filtered.tsv", sep="\t", quote=FALSE, row.names=TRUE)

# ----------------------------
# (E) MA plot with highlighted genes
# ----------------------------
resDF$plotColor <- "black"
resDF$plotColor[resDF$log2FoldChange>1 & resDF$padj<0.05] <- "red"
resDF$plotColor[resDF$log2FoldChange< -1 & resDF$padj<0.05] <- "deepskyblue"

genes_of_interest <- c("Pbx1","Pbx2","Zfhx3","Bmp4","Bmpr1a","Smad4","Id2","Itga3","Itga6")
selRows <- subset(resDF, gene_symbol %in% genes_of_interest)
selRows$nudge_y <- ifelse(selRows$log2FoldChange >=0, 1, -1)

png("MA_plot_highlight_genes_filtered.png", width=1200, height=900, res=150)
print(
  ggplot(resDF, aes(x=log10(baseMean+1), y=log2FoldChange)) +
    geom_point(alpha=0.3, size=1, aes(color=plotColor)) +
    scale_color_identity() +
    geom_point(data=selRows, color="green", size=3) +
    geom_text_repel(data=selRows, aes(label=gene_symbol),
                    nudge_y=selRows$nudge_y,
                    direction="y",
                    segment.color='grey50', size=4) +
    geom_hline(yintercept=0) +
    theme_bw()
)
dev.off()

# ----------------------------
# Volcano plot for filtered dataset
# ----------------------------
library(ggplot2)
library(ggrepel)

# attach a color column based on significance thresholds
resDF$plotColor <- "grey"
resDF$plotColor[resDF$log2FoldChange > 1 & resDF$padj < 0.05] <- "red"
resDF$plotColor[resDF$log2FoldChange < -1 & resDF$padj < 0.05] <- "deepskyblue"

# define genes to highlight
genes_of_interest <- c("Pbx1", "Pbx2", "Zfhx3", "Bmp4", "Bmpr1a", "Smad4", "Id2", "Itga3", "Itga6")
selRows <- subset(resDF, gene_symbol %in% genes_of_interest)

# create dynamic vertical nudge: positive log2FC → up, negative → down
if (nrow(selRows) > 0) {
  selRows$nudge_y <- ifelse(selRows$log2FoldChange >= 0, 1, -1)
}

# open PNG device
png("volcano_plot_filtered_highlight_genes.png", width = 1200, height = 900, res = 150)

# volcano plot
print(
  ggplot(resDF, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(color = plotColor), alpha = 0.6, size = 1.5) +
    scale_color_identity() +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
    xlab("Log2 Fold Change") + ylab("-log10(adj p-value)") +
    theme_minimal() +
    theme(legend.position = "bottom") +
    # highlight genes of interest
    geom_point(data = selRows, color = "green", size = 3) +
    geom_text_repel(
      data = selRows,
      aes(label = gene_symbol),
      nudge_y = selRows$nudge_y,
      direction = "y",
      segment.color = "grey50",
      size = 4
    )
)

dev.off()

